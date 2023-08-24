! This file is part of the Kestrel software for simulations
! of sediment-laden Earth surface flows.
!
! Version 1.0
!
! Copyright 2023 Mark J. Woodhouse, Jake Langham, (University of Bristol).
!
! This program is free software: you can redistribute it and/or modify it 
! under the terms of the GNU General Public License as published by the Free 
! Software Foundation, either version 3 of the License, or (at your option) 
! any later version.
!
! This program is distributed in the hope that it will be useful, but WITHOUT 
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
! FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for 
! more details.
!
! You should have received a copy of the GNU General Public License along with 
! this program. If not, see <https://www.gnu.org/licenses/>. 


! This module contains routines for reading input files that define
! a simulation.  The input file is passed as a command line argument.
! The input file is structured in blocks consisting of keyword = value
! declarations.  The blocks are handled by separate modules.
! The conditions are stored in a type(RunSet) :: RunParams object
! (see RunSettings.f90)

module read_input_file_module

   use set_precision_module, only: wp
   use messages_module, only: InfoMessage, WarningMessage, FatalErrorMessage
   use varstring_module, only: varString, ReadFileLine
   use utilities_module, only: pair, AddToVector
   use runsettings_module, only: RunSet, SourceStr, CapStr, CubeStr
   use domain_settings_module, only: DomainSettings_Set
   use initial_conditions_module, only: Source_Set, Cap_Set, Cube_Set
   use parameters_module, only: Params_Set
   use solver_settings_module, only: Solver_Set
   use outputsettings_module, only: Output_Set
   use topog_settings_module, only: Topog_Set

   implicit none

   private
   public :: ReadInputFile

contains

   ! ReadInputFile reads the input file called at the command line and builds
   ! the RunParams structures.
   ! Inputs: none -- input file specified at the command line.
   ! Outputs: RunParams [type RunSet] - the run settings
   subroutine ReadInputFile(RunParams)

      implicit none

      type(RunSet), intent(out) :: RunParams
      logical :: FileExists
      integer :: stat ! Status of input/output
      type(varString) :: line

      character(len=4096) :: InputFile

      character(len=8) :: date
      character(len=10) :: time

      call date_and_time(date=date, time=time)

      RunParams%RunDate%s = date
      RunParams%RunTime%s = time

      call get_command_argument(1, InputFile)
      RunParams%InputFile%s = trim(InputFile)
      inquire(file=RunParams%InputFile%s, exist=FileExists)
      if (.not.FileExists) call FatalErrorMessage('The requested input file ' // RunParams%InputFile%s // ' cannot be found.')
      call InfoMessage("Reading input file " // RunParams%InputFile%s)
      open(51,file=RunParams%InputFile%s) ! Open input file

      ! Read the first line
      ! read(51,'(A)',iostat=stat) line ! Read line of input file
      call ReadFileLine(51, line, iostat=stat)
      if (stat/=0) call FatalErrorMessage('Could not read the input file' // RunParams%InputFile%s) ! If iostat returns read error, exit loop
      close(51)
      line = line%adjustl() ! Ignore leading blanks
      call parseListInputs(RunParams)

   end subroutine ReadInputFile

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function linetype(s) result(t)
      ! Classify the line type of string s
      ! Inputs: s [character string] -- the input line
      ! Output: t [integer] -- the line type flag with:
      !                        t = 0: blank line
      !                        t = -1: comment
      !                        t = 1: block name (of the form "Block:")
      !                        t = 2: keyword-value pair (of the form "keyword=value")
      !                        t = 3: keyword (of the form "keyword")

      implicit none

      type(varString), intent(in) :: s
      integer :: t

      if (s%len()==0) then ! blank line
         t = 0
         return
      else if ((s%first_char()=='%').or.(s%first_char()=='#')) then ! comment line
         t = -1
         return
      else if (s%contains(':')) then ! Block name
         t = 1
         return
      else if (s%contains('=')) then ! Keyword-value pair
         t = 2
         return
      else ! assume it is a keyword
         t = 3
         return
      end if

   end function linetype

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine parseListInputs(RunParams)

      ! parseListInputs reads the input file called at the command line and builds the
      ! RunParams structures.
      ! Inputs: RunParams [type RunSet] - the run settings at input contains the input file name
      ! Outputs: RunParams [type RunSet] - the run settings is built from the input file

      implicit none

      type(RunSet), intent(inout) :: RunParams

      logical :: FileExists
      integer :: stat ! Status of input/output

      logical :: FoundDomain ! Domain is a required block in input file
      logical :: FoundTopog ! Topog is a required block in input file

      logical :: FoundIC     ! Initial conditions are required
      logical :: FoundSource, FoundCap, FoundCube ! Either a source, a cap or a cube input is required

      logical :: FoundParams ! Parameters not required, but will need to set specified parameters from input file
      logical :: FoundSolver ! Solver not required, but will need to set specified solver settings from input file
      logical :: FoundOutput ! Output not required, but will need to set specified output settings from input file

    !   character(len=MaxLineLength) :: line
      type(varString) :: line

      ! BlockNames are the block names in input file
      integer, parameter :: BlocksNumber = 8
      type(varString), dimension(:) :: InputBlocks(BlocksNumber)
      type(varString) :: currentBlock

      type(varString) :: variableName
      type(varString) :: variableValue

      ! Block Labels
      type(varString), dimension(:), allocatable :: DomainLabels ! Domain block labels
      type(varString), dimension(:), allocatable :: ParamLabels ! Parameters block labels
      type(varString), dimension(:), allocatable :: SolverLabels ! Solver block labels
      type(varString), dimension(:), allocatable :: OutputLabels ! Output block labels
      type(varString), dimension(:), allocatable :: TopogLabels ! Topog labels
      type(SourceStr) :: SourceStrings
      type(CapStr) :: CapStrings
      type(CubeStr) :: CubeStrings
      integer :: Nsources

      ! Block values
      type(varString), dimension(:), allocatable :: DomainValues
      type(varString), dimension(:), allocatable :: ParamValues
      type(varString), dimension(:), allocatable :: SolverValues
      type(varString), dimension(:), allocatable :: OutputValues
      type(varString), dimension(:), allocatable :: TopogValues

      FoundDomain = .FALSE.
      FoundIC = .FALSE.
      FoundSource = .FALSE.
      FoundCap = .FALSE.
      FoundCube = .FALSE.
      FoundParams = .FALSE.
      Nsources = 0

      ! Blocks
      InputBlocks = [varString('Domain'), &
         varString('Source'), &
         varString('Cap'), &
         varString('Cube'), &
         varString('Parameters'), &
         varString('Solver'), &
         varString('Topog'), &
         varString('Output')]

      inquire(File=RunParams%InputFile%s, Exist=FileExists)
      if (.not.FileExists) call FatalErrorMessage('The requested input file ' // RunParams%InputFile%s // ' cannot be found.')
      open(51,file=RunParams%InputFile%s) ! Open input file

      do ! Run through input file
         call ReadFileLine(51, line, iostat=stat)
         ! read(51,'(A)',iostat=stat) line ! Read line of input file
         if (stat/=0) exit ! If iostat returns read error, exit loop
         line = line%adjustl() ! Ignore leading blanks
         line = line%drop_comments() ! Drop comments (following #)

         if (linetype(line)==1) then ! this line is a Block name
            line = line%drop_comments(comment_char_in=":")
            ! Does the line contain a Block keyword?
            if (line%in_list(InputBlocks)) then
               currentBlock=line
               if (line=='Source') call incrementSourceStr(SourceStrings)
               if (line=='Cap') call incrementCapStr(CapStrings)
               if (line=='Cube') call incrementCubeStr(CubeStrings)
            end if

         else if (linetype(line)==2) then ! this line is a keyword=value pair
            call line%split("=", variableName, remain=variableValue)

            ! If currentBlock is "Domain", read the Domain keyword=value
            if (currentBlock=='Domain') then
               FoundDomain = .TRUE.
               call addToVector(DomainLabels,variableName)
               call addToVector(DomainValues,variableValue)
            end if

            ! If currentBlock is "Source", read the Source keyword=value
            if (currentBlock=='Source') then
               FoundSource = .TRUE.
               if (variableName=='sourceX') call updateSourceStr(SourceStrings,'x',variableValue)
               if (variableName=='sourceY') call updateSourceStr(SourceStrings,'y',variableValue)
               if (variableName=='sourceLon') call updateSourceStr(SourceStrings,'Lon',variableValue)
               if (variableName=='sourceLat') call updateSourceStr(SourceStrings,'Lat',variableValue)
               if (variableName=='sourceRadius') call updateSourceStr(SourceStrings,'Radius',variableValue)
               if (variableName=='sourceFlux') call updateSourceStr(SourceStrings,'flux',variableValue)
               if (variableName=='sourceTime') call updateSourceStr(SourceStrings,'time',variableValue)
               if (variableName=='sourceConc') call updateSourceStr(SourceStrings,'Conc',variableValue)
            end if

            ! If currentBlock is "Cap", read the Cap keyword=value
            if (currentBlock=='Cap') then
               FoundCap = .TRUE.

               if (variableName=='capX')      call updateCapStr(CapStrings,'x',variableValue)
               if (variableName=='capY')      call updateCapStr(CapStrings,'y',variableValue)
               if (variableName=='capLat')    call updateCapStr(CapStrings,'Lat',variableValue)
               if (variableName=='capLon')    call updateCapStr(CapStrings,'Lon',variableValue)
               if (variableName=='capRadius') call updateCapStr(CapStrings,'Radius',variableValue)
               if (variableName=='capHeight') call updateCapStr(CapStrings,'Height',variableValue)
               if (variableName=='capVolume') call updateCapStr(CapStrings,'Volume',variableValue)
               if (variableName=='capU')      call updateCapStr(CapStrings,'U',variableValue)
               if (variableName=='capV')      call updateCapStr(CapStrings,'V',variableValue)
               if (variableName=='capConc')   call updateCapStr(CapStrings,'Conc',variableValue)
               if (variableName=='capShape')  call updateCapStr(CapStrings,'Shape',variableValue)
            end if

            ! If currentBlock is "Cube", read the Cube keyword=value
            if (currentBlock=='Cube') then
               FoundCube = .TRUE.

               if (variableName=='cubeX')      call updateCubeStr(CubeStrings,'x',variableValue)
               if (variableName=='cubeY')      call updateCubeStr(CubeStrings,'y',variableValue)
               if (variableName=='cubeLat')    call updateCubeStr(CubeStrings,'Lat',variableValue)
               if (variableName=='cubeLon')    call updateCubeStr(CubeStrings,'Lon',variableValue)
               if (variableName=='cubeLength') then
                  call updateCubeStr(CubeStrings,'Length',variableValue)
               end if
               if (variableName=='cubeHeight') call updateCubeStr(CubeStrings,'Height',variableValue)
               if (variableName=='cubeWidth')  call updateCubeStr(CubeStrings,'Width',variableValue)
               if (variableName=='cubeU')      call updateCubeStr(CubeStrings,'U',variableValue)
               if (variableName=='cubeV')      call updateCubeStr(CubeStrings,'V',variableValue)
               if (variableName=='cubeConc')   call updateCubeStr(CubeStrings,'Conc',variableValue)
               if (variableName=='cubeShape')  call updateCubeStr(CubeStrings,'Shape',variableValue)
            end if

            ! If currentBlock is "Parameters", read the Parameters keyword=value
            if (currentBlock=='Parameters') then
               FoundParams = .TRUE.
               call addToVector(ParamLabels,variableName)
               call addToVector(ParamValues,variableValue)
            end if

            ! If currentBlock is "Solver", read the Solver keyword=value
            if (currentBlock=='Solver') then
               FoundSolver = .TRUE.
               call addToVector(SolverLabels,variableName)
               call addToVector(SolverValues,variableValue)
            end if

            ! If currentBlock is "Output", read the Output keyword=value
            if (currentBlock=='Output') then
               FoundOutput = .TRUE.
               call addToVector(OutputLabels,variableName)
               call addToVector(OutputValues,variableValue)
            end if

            ! If currentBlock is "Topog", read the Topog keyword=value
            if (currentBlock=='Topog') then
               FoundTopog = .TRUE.
               call addToVector(TopogLabels,variableName)
               call addToVector(TopogValues,variableValue)
            end if

         end if
      end do

      close(51)

      if (FoundDomain) then
         call DomainSettings_Set(DomainLabels,DomainValues,RunParams) ! Use DomainLabels and DomainValues to populate the RunParams structure
         deallocate(DomainLabels,DomainValues) ! Don't need DomainLabels and DomainValues anymore
      else
         call FatalErrorMessage("Required block 'Domain' not found in input file "// RunParams%InputFile%s)
      end if

      if (FoundTopog) then
         call Topog_Set(TopogLabels,TopogValues,RunParams) ! Use TopogLabels and TopogValues to populate the RunParams structure
         deallocate(TopogLabels,TopogValues) ! Don't need TopogLabels and TopogValues anymore so deallocate
      else
         call FatalErrorMessage("Required block 'Topog' not found in input file "// RunParams%InputFile%s)
      end if

      FoundIC = FoundSource .or. FoundCap .or. FoundCube
      if (FoundIC) then
         if (FoundCap) then
            call Cap_Set(CapStrings,RunParams) ! Use CapStrings to populate the RunParams structure
         else
            RunParams%nCaps = 0
         end if
         if (FoundCube) then
            call Cube_Set(CubeStrings,RunParams) ! Use CubeStrings to populate the RunParams structure
         else
            RunParams%nCubes = 0
         end if
         if (FoundSource) then
            call Source_Set(SourceStrings,RunParams) ! Use SourceStrings to populate the RunParams structure
         else
            RunParams%nSources = 0
         end if
      else
         call WarningMessage("An 'Initial Conditions' block does not contain any block variables.")
      end if

      if (FoundParams) then
         call Params_Set(ParamLabels,ParamValues,RunParams) ! Use ParamLabels and ParamValues to populate the RunParams structure
         deallocate(ParamLabels,ParamValues) ! Don't need ParamLabels and ParamValues anymore so deallocate
      end if

      if (FoundSolver) then
         call Solver_Set(SolverLabels,SolverValues,RunParams) ! Use SolverLabels and SolverValues to populate the RunParams structure
         deallocate(SolverLabels,SolverValues) ! Don't need SolverLabels and SolverValues anymore so deallocate
      end if

      if (FoundOutput) then
         call Output_Set(OutputLabels,OutputValues,RunParams) ! Use OutputLabels and OutputValues to populate the RunParams structure
         deallocate(OutputLabels,OutputValues) ! Don't need OutputLabels and OutputValues anymore so deallocate
      end if

    !   call DomainSettings_Disp(RunParams)

      return

   end subroutine parseListInputs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine incrementCapStr(CapString)

      implicit none

      ! Increment the CapString structure by initializing a new entry to zero
      type(CapStr), intent(inout) :: CapString

      call addToVector(CapString%x,'None')
      call addToVector(CapString%y,'None')
      call addToVector(CapString%Lat,'None')
      call addToVector(CapString%Lon,'None')
      call addToVector(CapString%Radius,'None')
      call addToVector(CapString%Height,'None')
      call addToVector(CapString%Volume,'None')
      call addToVector(CapString%U,'None')
      call addToVector(CapString%V,'None')
      call addToVector(CapString%psi,'None')
      call addToVector(CapString%Shape,'None')

      return
   end subroutine incrementCapStr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine updateCapStr(CapString,var,val)

      implicit none

      ! Update the CapString structure by initializing a new entry to zero
      type(CapStr), intent(inout) :: CapString
      character(len=*), intent(in) :: var
      type(varString), intent(in) :: val
      integer :: N

      N = size(CapString%Radius)

      select case (var)
       case ('x')
         CapString%x(N) = val
       case ('y')
         CapString%y(N) = val
       case ('Lat')
         CapString%Lat(N) = val
       case ('Lon')
         CapString%Lon(N) = val
       case ('Radius')
         CapString%Radius(N) = val
       case ('Height')
         CapString%Height(N) = val
       case ('Volume')
         CapString%Volume(N) = val
       case ('U')
         CapString%U(N) = val
       case ('V')
         CapString%V(N) = val
       case ('Conc')
         CapString%psi(N) = val
       case ('Shape')
         CapString%Shape(N) = val
       case default
         call WarningMessage("in updateCapStr, variable " // var // " not recognized.")
      end select

      return
   end subroutine updateCapStr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine incrementCubeStr(CubeString)

      implicit none

      ! Increment the CubeString structure by initializing a new entry to zero
      type(CubeStr), intent(inout) :: CubeString

      call addToVector(CubeString%x,'None')
      call addToVector(CubeString%y,'None')
      call addToVector(CubeString%Lat,'None')
      call addToVector(CubeString%Lon,'None')
      call addToVector(CubeString%Length,'None')
      call addToVector(CubeString%Height,'None')
      call addToVector(CubeString%Width,'None')
      call addToVector(CubeString%U,'None')
      call addToVector(CubeString%V,'None')
      call addToVector(CubeString%psi,'None')
      call addToVector(CubeString%Shape,'None')

      return
   end subroutine incrementCubeStr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine updateCubeStr(CubeString,var,val)

      implicit none

      ! Update the CubeString structure by initializing a new entry to zero
      type(CubeStr), intent(inout) :: CubeString
      character(len=*), intent(in) :: var
      type(varString), intent(in) :: val
      integer :: N

      N = size(CubeString%Length)

      select case (var)
       case ('x')
         CubeString%x(N) = val
       case ('y')
         CubeString%y(N) = val
       case ('Lat')
         CubeString%Lat(N) = val
       case ('Lon')
         CubeString%Lon(N) = val
       case ('Length')
         CubeString%Length(N) = val
       case ('Height')
         CubeString%Height(N) = val
       case ('Width')
         CubeString%Width(N) = val
       case ('U')
         CubeString%U(N) = val
       case ('V')
         CubeString%V(N) = val
       case ('Conc')
         CubeString%psi(N) = val
       case ('Shape')
         CubeString%Shape(N) = val
       case default
         call WarningMessage("in updateCubeStr, variable " // var // " not recognized.")
      end select

      return
   end subroutine updateCubeStr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine incrementSourceStr(SourceString)

      implicit none

      ! Increment the CapString structure by initializing a new entry to zero
      type(SourceStr), intent(inout) :: SourceString

      call addToVector(SourceString%x,'None')
      call addToVector(SourceString%y,'None')
      call addToVector(SourceString%Lat,'None')
      call addToVector(SourceString%Lon,'None')
      call addToVector(SourceString%Radius,'None')
      call addToVector(SourceString%flux,'None')
      call addToVector(SourceString%time,'None')
      call addToVector(SourceString%psi,'None')

      return
   end subroutine incrementSourceStr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine updateSourceStr(SourceString,var,val)

      implicit none

      ! Update the CapString structure by initializing a new entry to zero
      type(SourceStr), intent(inout) :: SourceString
      character(len=*), intent(in) :: var
      type(varString), intent(in) :: val
      integer :: N

      N = size(SourceString%Radius)

      select case (var)
       case ('x')
         SourceString%x(N) = val
       case ('y')
         SourceString%y(N) = val
       case ('Lat')
         SourceString%Lat(N) = val
       case ('Lon')
         SourceString%Lon(N) = val
       case ('Radius')
         SourceString%Radius(N) = val
       case ('flux')
         SourceString%flux(N) = val
       case ('time')
         SourceString%time(N) = val
       case ('Conc')
         SourceString%psi(N) = val
       case default
         call WarningMessage("in updateSourceStr, variable " // var // " not recognized.")
      end select

      return
   end subroutine updateSourceStr

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module read_input_file_module

