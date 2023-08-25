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


! This module assigns values to variables in RunSet associated with the model
! source conditions, as read from "Source"/"Cap"/"Cube" blocks of an input file.
! Some input validation is performed, producing Warning or FatalError messages.
module initial_conditions_module

   use set_precision_module, only: wp, pi
   use messages_module, only: WarningMessage, FatalErrorMessage
   use varstring_module, only: varString
   use utilities_module, only: pair
   use runsettings_module, only: RunSet, SourceStr, CapStr, CubeStr

   implicit none

   private
   public :: Cap_Set, Cube_Set, Source_Set

   character(len=4), parameter :: capShape_d = 'para'
   character(len=4), parameter :: cubeShape_d = 'flat'

contains

   ! Set RunParams related to flux source initial conditions in input file.
   ! Inputs: SourceString - "Source" Conditions;
   ! Output: RunParams - The run settings information.
   subroutine Source_Set(SourceString,RunParams)
      type(SourceStr), intent(in) :: SourceString
      type(RunSet), intent(inout) :: RunParams

      integer :: nSrc
      integer :: NsrcFlux, NsrcTime, Nsrcpsi

      type(pair) :: distxy
      type(pair) :: latlong

      logical :: set_srcX
      logical :: set_srcY
      logical :: set_srcRadius
      logical :: set_srcLat
      logical :: set_srcLon
      logical :: set_srcpsi
      logical :: set_srcFlux
      logical :: set_srcTime

      integer :: K

      RunParams%set_Sources=.TRUE.

      nSrc = size(SourceString%x)

      RunParams%nSources = nSrc
      allocate(RunParams%FluxSources(nSrc))

      do K=1,nSrc
         set_srcX = .FALSE.
         set_srcY = .FALSE.
         set_srcLat = .FALSE.
         set_srcLon = .FALSE.
         set_srcRadius = .FALSE.
         set_srcpsi = .FALSE.
         set_srcFlux = .FALSE.
         set_srcTime = .FALSE.
         if (SourceString%x(K)%s.NE.'None') then
            if (set_srcLon) then
               call Warning_TwoGiven("srcLon","srcX")
            else
               call setSourceRunParams(RunParams,SourceString,K,'x',set_srcX)
            end if
         end if
         if (SourceString%y(K)%s.NE.'None') then
            if (set_srcLat) then
               call Warning_TwoGiven("srcLat","srcY")
            else
               call setSourceRunParams(RunParams,SourceString,K,'y',set_srcY)
            end if
         end if
         if (SourceString%Lat(K)%s.NE.'None') then
            if (set_srcY) then
               call Warning_TwoGiven("srcY","srcLat")
            else
               call setSourceRunParams(RunParams,SourceString,K,'Lat',set_srcLat)
            end if
         end if
         if (SourceString%Lon(K)%s.NE.'None') then
            if (set_srcX) then
               call Warning_TwoGiven("srcX","srcLon")
            else
               call setSourceRunParams(RunParams,SourceString,K,'Lon',set_srcLon)
            end if
         end if
         if (SourceString%Radius(K)%s.NE.'None') then
            call setSourceRunParams(RunParams,SourceString,K,'Radius',set_srcRadius)
         end if

         if (SourceString%psi(K)%s.NE.'None') then
            set_srcpsi=.TRUE.
            Nsrcpsi = SourceString%psi(K)%count_substring(',')+1
            if (set_srcTime) then
               if (Nsrcpsi.ne.RunParams%FluxSources(K)%nFluxSeries) call FatalError_NotEqualSets(RunParams%InputFile%s,"sourceConc","sourceTime")
            else
               RunParams%FluxSources(K)%nFluxSeries = Nsrcpsi
            end if
            call SourceString%psi(K)%read_set(RunParams%FluxSources(K)%psi)
         end if

         if (SourceString%flux(K)%s.NE.'None') then
            set_srcFlux=.TRUE.
            NsrcFlux = SourceString%flux(K)%count_substring(',')+1
            if (set_srcTime) then
               if (NsrcFlux.ne.RunParams%FluxSources(K)%nFluxSeries) call FatalError_NotEqualSets(RunParams%InputFile%s,"sourceFlux","sourceTime")
            else
               RunParams%FluxSources(K)%nFluxSeries = NsrcFlux
            end if
            call SourceString%flux(K)%read_set(RunParams%FluxSources(K)%flux)
         end if

         if (SourceString%time(K)%s.NE.'None') then
            set_srcTime=.TRUE.
            NsrcTime = SourceString%Time(K)%count_substring(',')+1
            if (set_srcFlux) then
               if (NsrcTime.ne.RunParams%FluxSources(K)%nFluxSeries) call FatalError_NotEqualSets(RunParams%InputFile%s,"sourceTime","sourceFlux")
            else
               RunParams%FluxSources(K)%nFluxSeries = NsrcTime
            end if
            call SourceString%time(K)%read_set(RunParams%FluxSources(K)%time)
         end if

         if (.not. ((set_srcX.and.set_srcY).or.(set_srcLat.and.set_srcLon))) call FatalErrorMessage("In the 'Initial Conditions' block in the input file " &
            // trim(RunParams%InputFile%s) // new_line('A') &
            // "Nsource>0 but no position data (sourceX and sourceY or sourceLat and sourceLon) is given.")

         if (.not. set_srcRadius) call FatalErrorMessage("In the 'Initial Conditions' block in the input file " &
            // trim(RunParams%InputFile%s) // new_line('A') &
            // " Nsource>0 but sourceRadius is not given.")
         if (.not. set_srcpsi) call FatalErrorMessage("In the 'Initial Conditions' block in the input file " &
            // trim(RunParams%InputFile%s) // new_line('A') &
            // " Nsource>0 but sourceConcentration is not given.")
         if (.not. set_srcTime) call FatalErrorMessage("In the 'Initial Conditions' block in the input file " &
            // trim(RunParams%InputFile%s) // new_line('A') &
            // " Nsource>0 but sourceTime is not given.")
         if (.not. set_srcFlux) call FatalErrorMessage("In the 'Initial Conditions' block in the input file " &
            // trim(RunParams%InputFile%s) // new_line('A') &
            // " Nsource>0 but sourceFlux is not given.")

         if (.not.RunParams%isOneD) then

            if (RunParams%Georeference) then
                if ((set_srcLat).and.(set_srcLon)) then
                   distxy = RunParams%projTransformer%wgs84_to_utm(RunParams%FluxSources(K)%Lat, RunParams%FluxSources(K)%Lon)
                   RunParams%FluxSources(K)%x = distxy%first - RunParams%centerUTM%first
                   RunParams%FluxSources(K)%y = distxy%second - RunParams%centerUTM%second
                end if

                if ((set_srcX).and.(set_srcY)) then
                   distxy%first = RunParams%FluxSources(K)%x + RunParams%centerUTM%first
                   distxy%second = RunParams%FluxSources(K)%y + RunParams%centerUTM%second
                   latlong = RunParams%projTransformer%utm_to_wgs84(distxy%first, distxy%second)
                   RunParams%FluxSources(K)%Lat = latlong%first
                   RunParams%FluxSources(K)%Lon = latlong%second
                end if
            else
                if ((set_srcLat).and.(set_srcLon)) then
                    call FatalErrorMessage("In the 'Initial Conditions' block in the input file " &
                       // trim(RunParams%InputFile%s) // new_line('A') &
                       // " If the domain settings give a 2D simulation on artificial topography," // new_line('A') &
                       // " the source position must be given as sourceX and sourceY.")
                 end if
                 RunParams%FluxSources(:)%Lat = 0.0_wp
                 RunParams%FluxSources(:)%Lon = 0.0_wp
            end if

            if ((RunParams%FluxSources(K)%x-RunParams%FluxSources(K)%Radius)<-0.5_wp*RunParams%xSize) call FatalErrorMessage('In input file, Source ', K, 'outside of the domain')
            if ((RunParams%FluxSources(K)%x+RunParams%FluxSources(K)%Radius)>0.5_wp*RunParams%xSize) call FatalErrorMessage('In input file, Source ', K, 'outside of the domain')
            if ((RunParams%FluxSources(K)%y-RunParams%FluxSources(K)%Radius)<-0.5_wp*RunParams%ySize) call FatalErrorMessage('In input file, Source ', K, 'outside of the domain')
            if ((RunParams%FluxSources(K)%y+RunParams%FluxSources(K)%Radius)>0.5_wp*RunParams%ySize) call FatalErrorMessage('In input file, Source ', K, 'outside of the domain')
         else
            if ((set_srcLat).and.(set_srcLon)) then
               call FatalErrorMessage("In the 'Initial Conditions' block in the input file " &
                  // trim(RunParams%InputFile%s) // new_line('A') &
                  // " If the domain settings give a 1D simulation, the source position must be given as sourceX and sourceY.")
            end if
            RunParams%FluxSources(:)%Lat = 0.0_wp
            RunParams%FluxSources(:)%Lon = 0.0_wp
         end if
      end do
      RunParams%FluxSources(:)%NumCellsInSrc = 0

   end subroutine Source_Set

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Set RunParams related to a single flux source from a sourceString.
   ! This allows for multiple sources, contained in RunParams%FluxSources
   ! In: SourceString - "Source" Conditions;
   !     K - index of a flux source
   !     var - variable to set. Can be "x", "y", "Lat", "Lon", "Radius"
   ! InOut: RunParams - The run settings information.
   ! Out: set_flag - .TRUE. if input correctly set, otherwise .FALSE.
   subroutine setSourceRunParams(RunParams,SourceString,K,var,set_flag)
      type(RunSet), intent(inout) :: RunParams
      type(SourceStr), intent(in) :: SourceString
      integer, intent(in) :: K
      character(len=*), intent(in) :: var
      logical, intent(out) :: set_flag

      select case (var)
       case ('x')
         RunParams%FluxSources(K)%x = SourceString%x(K)%to_real()
         set_flag = .TRUE.
       case ('y')
         RunParams%FluxSources(K)%y = SourceString%y(K)%to_real()
         set_flag = .TRUE.
       case ('Lat')
        RunParams%FluxSources(K)%Lat = SourceString%Lat(K)%to_real()
         set_flag = .TRUE.
       case ('Lon')
        RunParams%FluxSources(K)%Lon = SourceString%Lon(K)%to_real()
         set_flag = .TRUE.
       case ('Radius')
         RunParams%FluxSources(K)%Radius = SourceString%Radius(K)%to_real()
         set_flag = .TRUE.
       case default
         call WarningMessage("In setSrcRunParams, variable " // var // " not recognized.")
         set_flag = .FALSE.
      end select

   end subroutine setSourceRunParams

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Set RunParams related to Cap initial conditions in input file.
   ! In: CapString - Cap Initial Conditions
   ! InOut: RunParams - The run settings information.
   subroutine Cap_Set(CapString,RunParams)
      type(CapStr), intent(in) :: CapString
      type(RunSet), intent(inout) :: RunParams

      integer :: Ncap

      type(pair) :: distxy
      type(pair) :: latlong

      logical :: set_capHeight
      logical :: set_capRadius
      logical :: set_capVol
      logical :: set_capX
      logical :: set_capY
      logical :: set_capLat
      logical :: set_capLon
      logical :: set_capU
      logical :: set_capV
      logical :: set_capShape
      logical :: set_cappsi

      integer :: K

      RunParams%set_Caps=.TRUE.

      Ncap = size(CapString%x)

      RunParams%Ncaps = Ncap
      call InitiateCapSources(RunParams, Ncap)

      do K=1,Ncap
         set_capX = .FALSE.
         set_capY = .FALSE.
         set_capLat = .FALSE.
         set_capLon = .FALSE.
         set_capVol = .FALSE.
         set_capHeight = .FALSE.
         set_capRadius = .FALSE.
         set_capu = .FALSE.
         set_capv = .FALSE.
         set_cappsi = .FALSE.
         set_capShape = .FALSE.
         if (CapString%x(K)%s.NE.'None') then
            if (set_capLon) then
               call Warning_TwoGiven("capLon","capX")
            else
               call setCapRunParams(RunParams,CapString,K,'x',set_capX)
            end if
         end if
         if (CapString%y(K)%s.NE.'None') then
            if (set_capLat) then
               call Warning_TwoGiven("capLat","capY")
            else
               call setCapRunParams(RunParams,CapString,K,'y',set_capY)
            end if
         end if
         if (CapString%Lat(K)%s.NE.'None') then
            if (set_capY) then
               call Warning_TwoGiven("capY","capLat")
            else
               call setCapRunParams(RunParams,CapString,K,'Lat',set_capLat)
            end if
         end if
         if (CapString%Lon(K)%s.NE.'None') then
            if (set_capX) then
               call Warning_TwoGiven("capX","capLon")
            else
               call setCapRunParams(RunParams,CapString,K,'Lon',set_capLon)
            end if
         end if
         if (CapString%Radius(K)%s.NE.'None') then
            if ((set_capVol).and.(set_capHeight)) then
               call Warning_3Given("capVolume","capHeight","capRadius")
            else
               call setCapRunParams(RunParams,CapString,K,'Radius',set_capRadius)
            end if
         end if
         if (CapString%Height(K)%s.NE.'None') then
            if ((set_capVol).and.(set_capRadius)) then
               call Warning_3Given("capVol","capRadius","capHeight")
            else
               call setCapRunParams(RunParams,CapString,K,'Height',set_capHeight)
            end if
         end if
         if (CapString%Volume(K)%s.NE.'None') then
            if ((set_capRadius).and.(set_capHeight)) then
               call Warning_3Given("capRadius","capHeight","capVolume")
            else
               call setCapRunParams(RunParams,CapString,K,'Volume',set_capVol)
            end if
         end if
         if (CapString%U(K)%s.NE.'None') call setCapRunParams(RunParams,CapString,K,'U',set_capU)
         if (CapString%V(K)%s.NE.'None') call setCapRunParams(RunParams,CapString,K,'V',set_capV)
         if (CapString%psi(K)%s.NE.'None') call setCapRunParams(RunParams,CapString,K,'conc',set_cappsi)
         if (CapString%Shape(K)%s.NE.'None') then
            call setCapRunParams(RunParams,CapString,K,'Shape',set_capShape)
         else
            RunParams%CapSources(K)%Shape = capShape_d
         end if

         if ((set_capRadius).and.(set_capHeight)) then
            if (RunParams%CapSources(K)%Shape == 'flat') then
               RunParams%CapSources(K)%Volume = pi*RunParams%CapSources(K)%Radius*RunParams%CapSources(K)%Radius*RunParams%CapSources(K)%Height
            elseif (RunParams%CapSources(K)%Shape == 'para') then
               RunParams%CapSources(K)%Volume = 0.5_wp*pi*RunParams%CapSources(K)%Radius*RunParams%CapSources(K)%Radius*RunParams%CapSources(K)%Height
            end if
         end if

         if ((set_capRadius).and.(set_capVol)) then
            if (RunParams%CapSources(K)%Shape == 'flat') then
               RunParams%CapSources(K)%Height = RunParams%CapSources(K)%Volume/pi/RunParams%CapSources(K)%Radius/RunParams%CapSources(K)%Radius
            elseif (RunParams%CapSources(K)%Shape == 'para') then
               RunParams%CapSources(K)%Height = 2.0_wp*RunParams%CapSources(K)%Volume/pi/RunParams%CapSources(K)%Radius/RunParams%CapSources(K)%Radius
            end if
         end if

         if ((set_capHeight).and.(set_capVol)) then
            if (RunParams%CapSources(K)%Shape == 'flat') then
               RunParams%CapSources(K)%Radius = sqrt(RunParams%CapSources(K)%Volume/pi/RunParams%CapSources(K)%Height)
            elseif (RunParams%CapSources(K)%Shape == 'para') then
               RunParams%CapSources(K)%Radius = sqrt(2.0_wp*RunParams%CapSources(K)%Volume/pi/RunParams%CapSources(K)%Height)
            end if
         end if

         if (RunParams%CapSources(K)%Shape == 'level' .and.  &
            ((.not. set_capHeight) .or. (.not. set_capRadius))) then
            call FatalErrorMessage("Must set height and radius for level cap")
         end if

         if (.not.RunParams%isOneD) then

            if (RunParams%Georeference) then
                if ((set_capLat).and.(set_capLon)) then
                   distxy = RunParams%projTransformer%wgs84_to_utm(RunParams%CapSources(K)%Lat, RunParams%CapSources(K)%Lon)
                   RunParams%CapSources(K)%x = distxy%first - RunParams%centerUTM%first
                   RunParams%CapSources(K)%y = distxy%second - RunParams%centerUTM%second
                end if

                if ((set_capX).and.(set_capY)) then
                   distxy%first = RunParams%CapSources(K)%x + RunParams%centerUTM%first
                   distxy%second = RunParams%CapSources(K)%y + RunParams%centerUTM%second
                   latlong = RunParams%projTransformer%utm_to_wgs84(distxy%first, distxy%second)
                   RunParams%CapSources(K)%Lat = latlong%first
                   RunParams%CapSources(K)%Lon = latlong%second
                end if
            else
                if ((set_capLat).and.(set_capLon)) then
                    call FatalErrorMessage("In the 'Initial Conditions' block in the input file " &
                       // trim(RunParams%InputFile%s) // new_line('A') &
                       // " If the domain settings give a 2D simulation on artificial topography," // new_line('A') &
                       // " the cap position must be given as capX and capY.")
                 end if
                 RunParams%CapSources(:)%Lat = 0.0_wp
                 RunParams%CapSources(:)%Lon = 0.0_wp
            end if

            if ((RunParams%CapSources(K)%x-RunParams%CapSources(K)%Radius)<-0.5_wp*RunParams%xSize) call FatalErrorMessage('In input file, cap ', K, 'outside of the domain')
            if ((RunParams%CapSources(K)%x+RunParams%CapSources(K)%Radius)>0.5_wp*RunParams%xSize) call FatalErrorMessage('In input file, cap ', K, 'outside of the domain')
            if ((RunParams%CapSources(K)%y-RunParams%CapSources(K)%Radius)<-0.5_wp*RunParams%ySize) call FatalErrorMessage('In input file, cap ', K, 'outside of the domain')
            if ((RunParams%CapSources(K)%y+RunParams%CapSources(K)%Radius)>0.5_wp*RunParams%ySize) call FatalErrorMessage('In input file, cap ', K, 'outside of the domain')
         else
            if ((set_capLat).and.(set_capLon)) then
               call FatalErrorMessage("In the 'Initial Conditions' block in the input file " &
                  // trim(RunParams%InputFile%s) // new_line('A') &
                  // " If the domain settings give a 1D simulation, the cap position must be given as capX and capY.")
            end if
            RunParams%CapSources(:)%Lat = 0.0_wp
            RunParams%CapSources(:)%Lon = 0.0_wp
         end if
      end do

   end subroutine Cap_Set

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Make empty containers for Cap sources
   ! InOut: RunParams - Run settings to update
   ! In: Ncaps - Number of cap sources to initialize
   subroutine InitiateCapSources(RunParams, Ncaps)
      type(RunSet), intent(inout) :: RunParams
      integer, intent(in) :: Ncaps

      integer :: K

      allocate(RunParams%CapSources(Ncaps))

      do K=1, Ncaps
         RunParams%CapSources(K)%x = 0.0_wp
         RunParams%CapSources(K)%y = 0.0_wp
         RunParams%CapSources(K)%Lat = 0.0_wp
         RunParams%CapSources(K)%Lon = 0.0_wp
         RunParams%CapSources(K)%Radius = 0.0_wp
         RunParams%CapSources(K)%Volume = 0.0_wp
         RunParams%CapSources(K)%Height = 0.0_wp
         RunParams%CapSources(K)%psi = 0.0_wp
         RunParams%CapSources(K)%u = 0.0_wp
         RunParams%CapSources(K)%v = 0.0_wp
         RunParams%CapSources(K)%Shape = 'flat'
      end do

   end subroutine InitiateCapSources

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Set the RunParams element var using CapStr structure for item K
   subroutine setCapRunParams(RunParams,CapString,K,var,set_flag)
      type(RunSet), intent(inout) :: RunParams
      type(CapStr), intent(in) :: CapString
      integer, intent(in) :: K
      character(len=*), intent(in) :: var
      logical, intent(out) :: set_flag

      type(varString) :: tmpStr

      select case (var)
       case ('x')
         RunParams%CapSources(K)%x = CapString%x(K)%to_real()
         set_flag = .TRUE.
       case ('y')
         RunParams%CapSources(K)%y = CapString%y(K)%to_real()
         set_flag = .TRUE.
       case ('Lat')
         RunParams%CapSources(K)%Lat = CapString%Lat(K)%to_real()
         set_flag = .TRUE.
       case ('Lon')
         RunParams%CapSources(K)%Lon = CapString%Lon(K)%to_real()
         set_flag = .TRUE.
       case ('Radius')
         RunParams%CapSources(K)%Radius = CapString%Radius(K)%to_real()
         set_flag = .TRUE.
       case ('Height')
         RunParams%CapSources(K)%Height = CapString%Height(K)%to_real()
         set_flag = .TRUE.
       case ('Volume')
         RunParams%CapSources(K)%Volume = CapString%Volume(K)%to_real()
         set_flag = .TRUE.
       case ('VelX','U')
         RunParams%CapSources(K)%U = CapString%U(K)%to_real()
         set_flag = .TRUE.
       case ('VelY','V')
         RunParams%CapSources(K)%V = CapString%V(K)%to_real()
         set_flag = .TRUE.
       case ('conc')
         RunParams%CapSources(K)%psi = CapString%psi(K)%to_real()
         set_flag = .TRUE.
       case ('Shape')
         tmpStr = CapString%Shape(K)%to_lower()
         select case (tmpStr%s)
          case ('parabolic','para')
            RunParams%CapSources(K)%Shape = 'para'
          case ('flat')
            RunParams%CapSources(K)%Shape = 'flat'
          case ('level')
            RunParams%CapSources(K)%Shape = 'level'
          case default
            call WarningMessage("In 'Initial Conditions' block 'CapShape' value not recognised.  " &
                // " Using default CapShape = " // capShape_d)
            RunParams%CapSources(K)%Shape = capShape_d
         end select
         set_flag = .TRUE.
       case default
         call WarningMessage("in setCapRunParams, variable " // var // " not recognized.")
         set_flag = .FALSE.
      end select

   end subroutine setCapRunParams

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Set RunParams related to "Cube" initial conditions in input file.
   ! In: CubeString - Cube Initial Conditions;
   ! InOut: RunParams - Updated run settings information.
   subroutine Cube_Set(CubeString,RunParams)
      type(CubeStr), intent(in) :: CubeString
      type(RunSet), intent(inout) :: RunParams

      integer :: Ncubes

      type(pair) :: distxy
      type(pair) :: latlong

      logical :: set_cubeHeight
      logical :: set_cubeWidth
      logical :: set_cubeLength
      logical :: set_cubeShape
      logical :: set_cubeX
      logical :: set_cubeY
      logical :: set_cubeLat
      logical :: set_cubeLon
      logical :: set_cubepsi
      logical :: set_cubeu
      logical :: set_cubev

      integer :: K

      RunParams%set_Cubes=.TRUE.

      Ncubes = size(CubeString%x)

      RunParams%Ncubes = Ncubes
      allocate(RunParams%CubeSources(Ncubes))

      do K=1,Ncubes
         set_cubeX = .FALSE.
         set_cubeY = .FALSE.
         set_cubeLat = .FALSE.
         set_cubeLon = .FALSE.
         set_cubeHeight = .FALSE.
         set_cubeLength = .FALSE.
         set_cubeWidth = .FALSE.
         set_cubeShape = .FALSE.
         set_cubeu = .FALSE.
         set_cubev = .FALSE.
         set_cubepsi = .FALSE.
         if (CubeString%x(K)%s.NE.'None') then
            if (set_cubeLon) then
               call Warning_TwoGiven("cubeLon","cubeX")
            else
               call setCubeRunParams(RunParams,CubeString,K,'x',set_cubeX)
            end if
         end if
         if (CubeString%y(K)%s.NE.'None') then
            if (set_cubeLat) then
               call Warning_TwoGiven("cubeLat","cubeY")
            else
               call setCubeRunParams(RunParams,CubeString,K,'y',set_cubeY)
            end if
         end if
         if (CubeString%Lat(K)%s.NE.'None') then
            if (set_cubeY) then
               call Warning_TwoGiven("cubeY","cubeLat")
            else
               call setCubeRunParams(RunParams,CubeString,K,'Lat',set_cubeLat)
            end if
         end if
         if (CubeString%Lon(K)%s.NE.'None') then
            if (set_cubeX) then
               call Warning_TwoGiven("cubeX","cubeLon")
            else
               call setCubeRunParams(RunParams,CubeString,K,'Lon',set_cubeLon)
            end if
         end if
         if (CubeString%Length(K)%s.NE.'None') then
            call setCubeRunParams(RunParams,CubeString,K,'Length',set_cubeLength)
         end if
         if (CubeString%Height(K)%s.NE.'None') then
            call setCubeRunParams(RunParams,CubeString,K,'Height',set_cubeHeight)
         end if
         if (CubeString%Width(K)%s.NE.'None') then
            call setCubeRunParams(RunParams,CubeString,K,'Width',set_cubeWidth)
         end if
         if (CubeString%U(K)%s.NE.'None') call setCubeRunParams(RunParams,CubeString,K,'U',set_cubeu)
         if (CubeString%V(K)%s.NE.'None') call setCubeRunParams(RunParams,CubeString,K,'V',set_cubev)
         if (CubeString%psi(K)%s.NE.'None') call setCubeRunParams(RunParams,CubeString,K,'conc',set_cubepsi)

         if (CubeString%Shape(K)%s.NE.'None') then
            call setCubeRunParams(RunParams,CubeString,K,'Shape',set_cubeShape)
         else
            RunParams%CubeSources(K)%Shape = cubeShape_d
         end if

         if (.not.RunParams%isOneD) then

            if (RunParams%Georeference) then
               if ((set_cubeLat).and.(set_cubeLon)) then
                  distxy = RunParams%projTransformer%wgs84_to_utm(RunParams%CubeSources(K)%Lat, RunParams%CubeSources(K)%Lon)

                  RunParams%CubeSources(K)%x = distxy%first - RunParams%centerUTM%first
                  RunParams%CubeSources(K)%y = distxy%second - RunParams%centerUTM%second
               end if

               if ((set_cubeX).and.(set_cubeY)) then
                  distxy%first = RunParams%CubeSources(K)%x + RunParams%centerUTM%first
                  distxy%second = RunParams%CubeSources(K)%y + RunParams%centerUTM%second
                  latlong = RunParams%projTransformer%utm_to_wgs84(distxy%first, distxy%second)
                  RunParams%CubeSources(K)%Lat = latlong%first
                  RunParams%CubeSources(K)%Lon = latlong%second
               end if
            else
               if ((set_cubeLat).and.(set_cubeLon)) then
                  call FatalErrorMessage("In the 'Initial Conditions' block in the input file " &
                     // trim(RunParams%InputFile%s) // new_line('A') &
                     // " If the domain settings give a 2D simulation on artificial topography," // new_line('A') &
                     // " the cube position must be given as cubeX and cubeY.")
               end if
               RunParams%CubeSources(:)%Lat = 0.0_wp
               RunParams%CubeSources(:)%Lon = 0.0_wp
            end if

            if ((RunParams%CubeSources(K)%x-0.5_wp*RunParams%CubeSources(K)%Length)<-0.5_wp*RunParams%xSize) call FatalErrorMessage('In input file, cube ', K, 'outside of the domain')
            if ((RunParams%CubeSources(K)%x+0.5_wp*RunParams%CubeSources(K)%Length)> 0.5_wp*RunParams%xSize) call FatalErrorMessage('In input file, cube ', K, 'outside of the domain')
            if ((RunParams%CubeSources(K)%y-0.5_wp*RunParams%CubeSources(K)%Width) <-0.5_wp*RunParams%ySize) call FatalErrorMessage('In input file, cube ', K, 'outside of the domain')
            if ((RunParams%CubeSources(K)%y+0.5_wp*RunParams%CubeSources(K)%Width) > 0.5_wp*RunParams%ySize) call FatalErrorMessage('In input file, cube ', K, 'outside of the domain')
         else
            if ((set_cubeLat).and.(set_cubeLon)) then
               call FatalErrorMessage("In the 'Initial Conditions' block in the input file " &
                  // trim(RunParams%InputFile%s) // new_line('A') &
                  // " If the domain settings give a 1D simulation, the cube position must be given as cubeX and cubeY.")
            end if
            RunParams%CubeSources(:)%Lat = 0.0_wp
            RunParams%CubeSources(:)%Lon = 0.0_wp
         end if
      end do

   end subroutine Cube_Set

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine setCubeRunParams(RunParams,CubeString,K,var,set_flag)
      ! Set the Run Params element var using CubeStr structure from item K
      implicit none

      type(RunSet), intent(inout) :: RunParams
      type(CubeStr), intent(in) :: CubeString
      integer, intent(in) :: K
      character(len=*), intent(in) :: var
      logical, intent(out) :: set_flag

      type(varString) :: tmpStr

      select case (var)
       case ('x')
         RunParams%CubeSources(K)%x  = CubeString%x(K)%to_real()
         set_flag = .TRUE.
       case ('y')
         RunParams%CubeSources(K)%y  = CubeString%y(K)%to_real()
         set_flag = .TRUE.
       case ('Lat')
         RunParams%CubeSources(K)%Lat = CubeString%Lat(K)%to_real()
         set_flag = .TRUE.
       case ('Lon')
         RunParams%CubeSources(K)%Lon = CubeString%Lon(K)%to_real()
         set_flag = .TRUE.
       case ('Length')
         RunParams%CubeSources(K)%Length  = CubeString%Length(K)%to_real()
         set_flag = .TRUE.
       case ('Height')
         RunParams%CubeSources(K)%Height  = CubeString%Height(K)%to_real()
         set_flag = .TRUE.
       case ('Width')
         RunParams%CubeSources(K)%Width  = CubeString%Width(K)%to_real()
         set_flag = .TRUE.
       case ('VelX','U')
         RunParams%CubeSources(K)%U  = CubeString%U(K)%to_real()
         set_flag = .TRUE.
       case ('VelY','V')
         RunParams%CubeSources(K)%V  = CubeString%V(K)%to_real()
         set_flag = .TRUE.
       case ('conc')
         RunParams%CubeSources(K)%psi  = CubeString%psi(K)%to_real()
         set_flag = .TRUE.
       case ('Shape')
         tmpStr = CubeString%Shape(K)%to_lower()
         select case (tmpStr%s)
          case ('flat')
            RunParams%CubeSources(K)%Shape = 'flat'
          case ('level')
            RunParams%CubeSources(K)%Shape = 'level'
          case default
            call WarningMessage("In 'Initial Conditions' block 'CubeShape' value not recognised.  " &
               // "Using default CubeShape = " // cubeShape_d)
            RunParams%CubeSources(K)%Shape = cubeShape_d
         end select
         set_flag = .TRUE.

       case default
         call WarningMessage("in setCubeRunParams, variable " // var // " not recognized.")
         set_flag = .FALSE.
      end select

   end subroutine setCubeRunParams

   subroutine FatalError_NotEqualSets(InputFile,var1,var2)

      implicit none

      character(len=*), intent(in) :: InputFile
      character(len=*), intent(in) :: var1
      character(len=*), intent(in) :: var2

      call FatalErrorMessage("In the 'Initial Conditions' block in the input file " &
         // trim(InputFile) // new_line('A') &
         // "the number of values in the sets of '" &
         // trim(var1) // "' and '" // trim(var2) // "' must be equal.")

   end subroutine FatalError_NotEqualSets

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine Warning_3Given(var1,var2,varIgnore)

      implicit none

      character(len=*), intent(in) :: var1
      character(len=*), intent(in) :: var2
      character(len=*), intent(in) :: varIgnore

      call WarningMessage("In the 'Initial Conditions' block the block variables " // trim(var1) // " and " // trim(var2) // " are set before " // trim(varIgnore) // "; " // trim(varIgnore) // " will be ignored.")

   end subroutine Warning_3Given

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine Warning_TwoGiven(var1,var2)

      implicit none

      character(len=*), intent(in) :: var1
      character(len=*), intent(in) :: var2

      call WarningMessage("The 'Initial Conditions' block" &
         // " contains both " // trim(var1) // " and " // trim(var2) // ".  Only one should be set. " &
         // "Using " // trim(var1) // " and ignoring " // trim(var2) // ".")

      return
   end subroutine Warning_TwoGiven

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module initial_conditions_module
