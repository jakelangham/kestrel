! This file is part of the Kestrel software for simulations
! of sediment-laden Earth surface flows.
!
! Version v1.1.1
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


! This module contains the primary routines for integrating general governing
! equations for shallow layer flow models with bed evolution (morphodynamics) 
! and a suspended sediment load. Our method employs the method of lines, i.e. 
! we discretise the governing equations in space (using finite volumes) and 
! time step the resulting ODE system, which has the form
!
! dq/dt = H(q) + M(q)
!
! where q = (w, Hnpsi, rhoHnu, rhoHnv, b) is a vector of the primary flow 
! variables [for definitions, see main.f90], while H(q) and M(q) are (spatially 
! discretised) PDE operators respectively associated with the hydraulic and 
! morphodynamic physics. Their contributions are computed in 
! CalculateHydraulicRHS [see HydraulicRHS.f90] and CalculateMorphodynamicRHS 
! [see MorphodynamicRHS.f90], respectively.
! 
! The time stepper used is a second-order semi-implicit Runge-Kutta scheme, 
! [given in Chertock, Cui, Kurganov & Wu, Int. J. Numer. Meth. Fluids (2015)]
! in which only the stiff basal drag term is treated implicitly. This allows
! updates at each substep to be rearranged into explicit formulae. 
! Additionally, we employ Strang operator splitting to time step the hydraulic
! and morphodynamic operators separately.
! 
! More details on our approach, including the model equations and the numerical
! scheme are available in the following paper: https://arxiv.org/abs/2306.16185
module timestepper_module

   use set_precision_module, only: wp
   use messages_module, only: InfoMessage, TimestepMessage, RevisedTimestepMessage, WarningMessage
   use grid_module, only: GridType, TileList, TileType
   use redistribute_module, only: AddToRedistList, ExcessDeposition, RedistributeGrid, RedistList
   use runsettings_module, only: RunSet
   use closures_module, only : ComputeHn, FlowSquaredSpeedSlopeAligned, GeometricCorrectionFactor
   use hydraulic_rhs_module, only: CalculateHydraulicRHS, ComputeDesingularisedVariables
   use morphodynamic_rhs_module, only: CalculateMorphodynamicRHS, ComputeCellCentredTopographicData, ComputeInterfacialTopographicData
   use update_tiles_module, only: AddTile
   use varstring_module, only: varString
   use output_module, only: OutputSolutionData, OutputAggregateData, OutputInfo, CalculateVolume
   use utilities_module, only: Int2String

   implicit none

   private
   public :: Run

contains

   ! Run the simulation, assuming RunParams has been loaded and the numerical
   ! grid has been set up correctly by the caller in main.f90. This routine
   ! mainly handles dividing the simulation time into separate integrations, 
   ! between which the appropriate output routine are called to record the
   ! simulation data.
   subroutine Run(RunParams, grid)

      implicit none

      type(RunSet), intent(inout) :: RunParams
      type(GridType), target, intent(inout) :: grid

      type(varString) :: RunInfoFilename

      real(kind=wp) :: timestep

      integer :: i
      
      RunInfoFilename = RunParams%OutDir + RunParams%InfoFilename

      RunParams%CurrentOut = RunParams%FirstOut

      if (.not. RunParams%Restart) then
         call OutputSolutionData(RunParams, RunParams%CurrentOut, grid)
         call OutputInfo(RunParams)
         call CalculateVolume(RunParams, RunParams%tstart, grid, create = .true.)
      end if

      ! Main loop over Nout output intervals of length RunParams%DeltaT. 
      ! We integrate from t = RunParams%tstart to RunParams%tend. If this is 
      ! a restarted simulation, then RunParams%FirstOutput > 0.
      do i = 1, RunParams%Nout - RunParams%FirstOut

         timestep = RunParams%tstart + (i + RunParams%FirstOut) * RunParams%DeltaT
         call IntegrateTo(RunParams, timestep, grid)

         RunParams%CurrentOut = RunParams%CurrentOut + 1

         call OutputSolutionData(RunParams, RunParams%CurrentOut, grid)
         call CalculateVolume(RunParams, timestep, grid)
         call OutputAggregateData(RunParams, grid)
         call OutputInfo(RunParams)

      end do

   end subroutine Run

   ! Integrate the solution data in the numerical grid from t = grid%t to tend.
   subroutine IntegrateTo(RunParams, tend, grid)

      implicit none

      type(RunSet), intent(inout) :: RunParams
      real(kind=wp), intent(in) :: tend
      type(GridType), target, intent(inout) :: grid

      real(kind=wp) :: advisedTimeStep, nextT, tmax
      real(kind=wp) :: dt_hydro, dt_morpho

      integer :: nTiles, nActiveTiles, nDimensions, nFlux

      logical :: integrating, firstStep, refineTimeStep

      nTiles = grid%nTiles
      nActiveTiles = grid%activeTiles%size
      nDimensions = RunParams%nDimensions
      nFlux = RunParams%nFlux

      if (tend <= grid%t) then
         call WarningMessage("End-time Tend is less than initial time")
         integrating = .false.
         return
      else
         integrating = .True.
      end if

      firstStep = .true.

      ! We integrate using Strang operator splitting method to separate the
      ! hydraulic and morphodynamic parts. If morphodynamics is not set
      ! we just take single time steps of the hydraulic operator H.
      ! Otherwise, we do H(dt) * M(2dt) * H(dt), where M is the morphodynamic
      ! operator.
      do while (integrating)
         call CheckIfNearBoundaries(RunParams, grid, nActiveTiles)
         call InitialiseTimeSteppingArrays(nActiveTiles, grid)

         call CalculateHydraulicRHS(RunParams, grid, grid%tileContainer, &
                                    grid%t, 1, advisedTimeStep)

         ! This ensures that we take time steps that exactly align with any flux
         ! source time series, which is needed to exactly deliver the total flux
         ! specified (and hence needed for conservativity).
         tmax = min(tend, NextFluxSeriesTime(RunParams, grid%t))

         if (.not. RunParams%MorphodynamicsOn) then
            advisedTimeStep = min(advisedTimeStep, tmax - grid%t)
         else
            advisedTimeStep = min(advisedTimeStep, 0.5_wp * (tmax - grid%t))
         end if
         grid%dt = advisedTimeStep
         dt_hydro = advisedTimeStep

         ! This loop tries to complete a single H(dt) * M(2dt) * H(dt) step.
         grid%t0 = grid%t ! t0 is the t to fall back to if dt gets revised.
         do while (.true.)
            ! Copy necessary data over to the intermed{0-3} arrays, which 
            ! store the solution vector at each substep of the method.
            call CopySolutionData(RunParams, grid, grid%tileContainer, grid%intermed0)
            call CopyMutableTopographicData(RunParams, grid, grid%tileContainer, grid%intermed1)
            call CopyMutableTopographicData(RunParams, grid, grid%tileContainer, grid%intermed2)
            call CopySolutionData(RunParams, grid, grid%tileContainer, grid%intermed3)

            call HydraulicTimeStepper(RunParams, grid, dt_hydro, nextT, refineTimeStep)

            if (refineTimeStep) then
#if DEBUG_TIMESTEP==1 || DEBUG_TIMESTEP==2
               call InfoMessage('Refining at first hydraulic substep')
#if DEBUG_TIMESTEP==2
               call exit
#endif
#endif
               cycle
            end if

            ! If no erosion, we're already done
            if (.not.RunParams%MorphodynamicsOn) exit
            ! else copy solution data back into I0 for the next operator.
            call CopySolutionData(RunParams, grid, grid%intermed3, grid%intermed0)
            call CopySolutionData(RunParams, grid, grid%intermed3, grid%intermed1)
            call CopySolutionData(RunParams, grid, grid%intermed3, grid%intermed2)
            dt_morpho = 2.0_wp * dt_hydro

            call MorphodynamicTimeStepper(RunParams, grid, dt_morpho, refineTimeStep)

            if (refineTimeStep) then
#if DEBUG_TIMESTEP==1 || DEBUG_TIMESTEP==2
               call InfoMessage('Refining at morphodynamic substep')
#if DEBUG_TIMESTEP==2
               call exit
#endif
#endif
               ! Set hydraulic dt to new suggested step.
               dt_hydro = 0.5_wp * dt_morpho
               grid%dt = dt_hydro
               cycle
            end if

            ! N.B. Updated topographic data needed for intermediate computations.
            call CopySolutionData(RunParams, grid, grid%intermed3, grid%intermed0)
            call CopyMutableTopographicData(RunParams, grid, grid%intermed3, grid%intermed1)
            call CopyMutableTopographicData(RunParams, grid, grid%intermed3, grid%intermed2)
            grid%t = grid%t + dt_hydro

            call CalculateHydraulicRHS(RunParams, grid, grid%intermed0, &
                                       grid%t, 1, advisedTimeStep)
            if (advisedTimeStep < dt_hydro) then
               ! Make sure dt changes by at least 10%
               if ((dt_hydro - advisedTimeStep) / dt_hydro < 0.1_wp) then
                  dt_hydro = 0.9_wp * dt_hydro
               else
                  dt_hydro = advisedTimeStep
               end if
               grid%dt = dt_hydro
               grid%t = grid%t0 ! otherwise need to reset t
#if DEBUG_TIMESTEP==1 || DEBUG_TIMESTEP==2
               call InfoMessage('Refining at second hydraulic substep')
#if DEBUG_TIMESTEP==2
               call exit
#endif
#endif
               call ReportRevisedTimeStep(RunParams, grid%t0, dt_hydro)
               cycle
            end if

            call HydraulicTimeStepper(RunParams, grid, dt_hydro, nextT, refineTimeStep)

            if (.not. refineTimeStep) then
               exit ! We're done
            end if
            grid%t = grid%t0 ! ...otherwise need to reset t.
#if DEBUG_TIMESTEP==1 || DEBUG_TIMESTEP==2
            call InfoMessage('Refining at second hydraulic substep')
#if DEBUG_TIMESTEP==2
            call exit
#endif
#endif
         end do

         ! Copy final time step back into main solution vector.
         call CopySolutionData(RunParams, grid, grid%intermed3, grid%tileContainer)

         ! Record previous simulation time and set current value.
         if (.not.RunParams%MorphodynamicsOn) then
            grid%t = grid%t0 + dt_hydro
         else
            grid%t = grid%t0 + 2.0_wp * dt_hydro
         end if

         ! Output simulation time.
         call TimestepMessage(grid%t, dt_hydro)
         
         if (grid%t >= tend) then
            integrating = .false.
         end if

         firstStep = .false.
      end do ! end integration

   end subroutine IntegrateTo

   ! Given t, look through the time series for flux sources and determine which
   ! is the next point in the future. Return in nextT.
   function NextFluxSeriesTime(RunParams, t) result(nextT)
      implicit none

      type(RunSet), intent(in) :: RunParams
      real(kind=wp), intent(in) :: t
      real(kind=wp) :: nextT

      integer :: i, j, nSrc, nFluxSeries
      real(kind=wp) :: tdiff, tmp

      nSrc = RunParams%nSources

      nextT = huge(1.0_wp)
      tdiff = huge(1.0_wp)
      do i = 1, nSrc
         nFluxSeries = RunParams%FluxSources(i)%nFluxSeries
         do j = 1, nFluxSeries
            tmp = RunParams%FluxSources(i)%time(j) - t
            if (tmp > 0.0_wp .and. tmp < tdiff) then
               tdiff = tmp
               nextT = RunParams%FluxSources(i)%time(j)
            end if
         end do
      end do
   end function NextFluxSeriesTime

   ! Output to let user know dt is going to be refined.
   subroutine ReportRevisedTimeStep(RunParams, t, dt)
      implicit none

      type(RunSet), intent(in) :: RunParams
      real(kind=wp), intent(in) :: t, dt

      real(kind=wp) :: nextT

      ! Next time depends on whether a morphodynamic step will be taken.
      nextT = t + dt
      if (RunParams%MorphodynamicsOn) nextT = nextT + dt

      call RevisedTimestepMessage(nextT, dt)
      
   end subroutine ReportRevisedTimeStep

   ! Try to compute a single time step of interval thisdt for the hydraulic
   ! subproblem. If thisdt is deemed to be too large for a stable time step to
   ! be taken (i.e. because the CFL condition is violated) 
   ! return refineTimeStep = .true.
   !
   ! N.B. CalculateHydraulicRHS should be called once prior to this routine in
   ! order to compute dt and gather the first RHS data, stored in 
   ! grid%{Explicit,Implicit}0_init or grid%{Explicit,Implicit}0 depending on if
   ! it's the first or second splitting step respectively.
   subroutine HydraulicTimeStepper(RunParams, grid, thisdt, &
                                   nextT, refineTimeStep)
      implicit none

      type(RunSet), intent(inout) :: RunParams
      type(GridType), target, intent(inout) :: grid
      real(kind=wp), intent(inout) :: thisdt
      real(kind=wp), intent(inout) :: nextT
      logical, intent(inout) :: refineTimeStep

      type(TileType), dimension(:), pointer :: tileContainer
      type(TileList), pointer :: activeTiles
      type(tileType), dimension(:), pointer :: intermed0, intermed1, intermed2, intermed3

      real(kind=wp), pointer :: t
      real(kind=wp) :: dt1, dt2, dt3
      integer :: tt, ttk, k, nFlux, var, i, j

      real(kind=wp), dimension(:,:), allocatable :: hp_old, hp_new, w_update

      ! These pointers simplify notation.
      tileContainer => grid%tileContainer
      activeTiles => grid%activeTiles
      intermed0 => grid%intermed0; intermed1 => grid%intermed1
      intermed2 => grid%intermed2; intermed3 => grid%intermed3

      ! Temp arrays needed for partial calculations.
      allocate(hp_old(RunParams%nXpertile, RunParams%nYpertile), &
               hp_new(RunParams%nXpertile, RunParams%nYpertile), &
               w_update(RunParams%nXpertile, RunParams%nYpertile))

      t => grid%t
      nextT = grid%t + thisdt
      grid%dt = thisdt
      nFlux = RunParams%nFlux
      refineTimeStep = .false.

      ! Compute first substep.
      do tt = 1, ActiveTiles%Size
         ttk = ActiveTiles%List(tt)

         call ApplySpongeLayer(RunParams, ttk, tileContainer, intermed0, .false.)

         do k = 1, nFlux
            var = RunParams%iFlux(k)
            if (RunParams%ImplicitStep(var) .and. (.not. RunParams%StoppedMaterialHandling)) then
               intermed1(ttk)%u(var,:,:) = (intermed0(ttk)%u(var,:,:) + &
                  thisdt * intermed0(ttk)%ddtExplicit(var,:,:)) / (1.0_wp - thisdt * intermed0(ttk)%ddtImplicit(var,:,:))
            else
               intermed1(ttk)%u(var, :, :) = intermed0(ttk)%u(var, :, :) + &
                  thisdt * intermed0(ttk)%ddtExplicit(var, :, :)
            end if
         end do
      end do ! end first substep

      ! Compute new RHS from first substep.
      call CalculateHydraulicRHS(RunParams, grid, intermed1, nextT, 2, dt1)

      ! Check that the first substep was ok. If too big, we recommend a new 
      ! dt and return.
      if (dt1 < thisdt) then
         thisdt = 0.9_wp * dt1
         call ReportRevisedTimeStep(RunParams, grid%t0, thisdt)
         refineTimeStep = .true.
#if DEBUG_TIMESTEP==1 || DEBUG_TIMESTEP==2
         call InfoMessage('HydraulicTimeStepper refining at first correction')
#if DEBUG_TIMESTEP==2
         call exit
#endif
#endif
         return
      end if

      ! Compute second substep.
      do tt = 1, ActiveTiles%Size
         ttk = ActiveTiles%List(tt)

         call ApplySpongeLayer(RunParams, ttk, tileContainer, intermed1, .false.)

         do k = 1, nFlux
            var = RunParams%iFlux(k)
            ! If using stopped material handling, our current method relies 
            ! overriding the default time stepper with (2nd order) Heun's 
            ! method. So that's what's happening here.
            if (RunParams%StoppedMaterialHandling) then
               do i = 1, RunParams%nXpertile
                  do j = 1, RunParams%nYpertile
                     if (abs(intermed1(ttk)%u(var,i,j)) < 1d-14) then
                        intermed3(ttk)%u(var,i,j) = intermed1(ttk)%u(var,i,j)
                     else
                        intermed3(ttk)%u(var,i,j) = 0.5_wp * &
                           (intermed0(ttk)%u(var,i,j) + &
                           intermed1(ttk)%u(var,i,j) + &
                           thisdt * intermed1(ttk)%ddtExplicit(var,i,j))
                     end if
                  end do
               end do
            else if (RunParams%ImplicitStep(var)) then
               intermed2(ttk)%u(var,:,:) = &
                  0.75_wp * intermed0(ttk)%u(var,:,:) + &
                  0.25_wp * (intermed1(ttk)%u(var,:,:) + &
                  thisdt * intermed1(ttk)%ddtExplicit(var,:,:)) / (1.0_wp - thisdt * intermed1(ttk)%ddtImplicit(var,:,:))
            else if (var == RunParams%Vars%w) then
               ! This way of doing the explicit update for w makes sure that
               ! if hp = 0 (for intermed{0,1}) and ddtExplicit1 = 0 then
               ! there is no update. Otherwise, finite precision errors when
               ! splitting up the expression can lead to (small) negative depths
               ! when flow is very shallow. Same for the intermed3 update below.
               hp_old = -intermed0(ttk)%u(RunParams%Vars%bt, :, :)
               hp_old = hp_old + (intermed0(ttk)%u(RunParams%Vars%w, :, :) - &
                                  intermed0(ttk)%u(RunParams%Vars%b0, :, :))
               hp_new = -intermed1(ttk)%u(RunParams%Vars%bt, :, :)
               hp_new = hp_new + (intermed1(ttk)%u(RunParams%Vars%w, :, :) - &
                                  intermed1(ttk)%u(RunParams%Vars%b0, :, :))
               w_update = intermed0(ttk)%u(RunParams%Vars%bt, :, :)
               w_update = w_update + 0.25_wp * hp_new
               w_update = w_update + 0.75_wp * hp_old
               w_update = w_update + 0.25_wp * thisdt * intermed1(ttk)%ddtExplicit(var, :, :)
               w_update = w_update + intermed0(ttk)%u(RunParams%Vars%b0, :, :)
               intermed2(ttk)%u(var, :, :) = w_update
            else
               intermed2(ttk)%u(var,:,:) = &
                  0.75_wp * intermed0(ttk)%u(var,:,:) + &
                  0.25_wp * (intermed1(ttk)%u(var,:,:) + &
                  thisdt * intermed1(ttk)%ddtExplicit(var,:,:))
            end if
         end do
      end do ! end second substep

      if (.not. RunParams%StoppedMaterialHandling) then
         ! Compute new RHS from second substep.
         call CalculateHydraulicRHS(RunParams, grid, intermed2, &
                                    t + 0.5_wp * thisdt, 3, dt2)
   
         ! Check that the second substep was ok. If too big, we recommend a new 
         ! dt and return.
         if (dt2 < thisdt) then
            thisdt = 0.9_wp * dt2
            call ReportRevisedTimeStep(RunParams, grid%t0, thisdt)
            refineTimeStep = .true.
#if DEBUG_TIMESTEP==1 || DEBUG_TIMESTEP==2
            call InfoMessage('HydraulicTimeStepper refining at second correction')
#if DEBUG_TIMESTEP==2
            call exit
#endif
#endif
            return
         end if

         ! Compute third substep.
         do tt = 1, ActiveTiles%Size
            ttk = ActiveTiles%List(tt)

            call ApplySpongeLayer(RunParams, ttk, tileContainer, intermed2, .false.)

            do k = 1, nFlux
               var = RunParams%iFlux(k)
               if (RunParams%ImplicitStep(var)) then
                  intermed3(ttk)%u(var,:,:) = &
                     (1.0_wp/3.0_wp) * intermed0(ttk)%u(var,:,:) + &
                     (2.0_wp/3.0_wp) * (intermed2(ttk)%u(var,:,:) + &
                     thisdt * intermed2(ttk)%ddtExplicit(var,:,:)) / (1.0_wp - thisdt * intermed2(ttk)%ddtImplicit(var,:,:))
               else if (var == RunParams%Vars%w) then
                  hp_old = -intermed0(ttk)%u(RunParams%Vars%bt, :, :)
                  hp_old = hp_old + (intermed0(ttk)%u(RunParams%Vars%w, :, :) - &
                                     intermed0(ttk)%u(RunParams%Vars%b0, :, :))
                  hp_new = -intermed2(ttk)%u(RunParams%Vars%bt, :, :)
                  hp_new = hp_new + (intermed2(ttk)%u(RunParams%Vars%w, :, :) - &
                                     intermed2(ttk)%u(RunParams%Vars%b0, :, :))
                  w_update = intermed0(ttk)%u(RunParams%Vars%bt, :, :)
                  w_update = w_update + (1.0_wp/3.0_wp) * hp_old
                  w_update = w_update + (2.0_wp/3.0_wp) * hp_new
                  w_update = w_update + (2.0_wp/3.0_wp) * thisdt * intermed2(ttk)%ddtExplicit(var, :, :)
                  w_update = w_update + intermed0(ttk)%u(RunParams%Vars%b0, :, :)
                  intermed3(ttk)%u(var, :, :) = w_update
               else
                  intermed3(ttk)%u(var,:,:) = &
                     (1.0_wp/3.0_wp) * intermed0(ttk)%u(var,:,:) + &
                     (2.0_wp/3.0_wp) * (intermed2(ttk)%u(var,:,:) + &
                     thisdt * intermed2(ttk)%ddtExplicit(var,:,:))
               end if
            end do
         end do ! end third substep

         ! Compute new RHS from third substep.
         call CalculateHydraulicRHS(RunParams, grid, intermed3, nextT, 4, dt3)

      end if

      ! Final implicit substep.
      do tt = 1, ActiveTiles%Size
         ttk = ActiveTiles%List(tt)

         if (.not. RunParams%StoppedMaterialHandling) then
            call ApplySpongeLayer(RunParams, ttk, tileContainer, intermed3, .false.)

            do k = 1, nFlux
               var = RunParams%iFlux(k)
               if (RunParams%ImplicitStep(var)) then
                  intermed3(ttk)%u(var,:,:) = &
                     (intermed3(ttk)%u(var,:,:) - &
                     thisdt * thisdt * intermed3(ttk)%ddtExplicit(var,:,:) * intermed3(ttk)%ddtImplicit(var,:,:)) / &
                     (1.0_wp + thisdt * thisdt * intermed3(ttk)%ddtImplicit(var,:,:) * intermed3(ttk)%ddtImplicit(var,:,:))
               end if
            end do
         end if

         call ComputeDesingularisedVariables(RunParams, tileContainer, ttk, .true.)
         call UpdateMaximumHeights(RunParams, tileContainer(ttk), nextT)
         call UpdateMaximumSpeeds(RunParams, tileContainer(ttk), nextT)
         call UpdateMaximumErosion(RunParams, tileContainer(ttk), nextT)
         call UpdateMaximumDeposit(RunParams, tileContainer(ttk), nextT)
         call UpdateMaximumSolidsFraction(RunParams, tileContainer(ttk), nextT)
      end do ! end final implicit substep

   end subroutine HydraulicTimeStepper

   ! Try to compute a single time step of interval thisdt for the morphodynamic
   ! subproblem. If thisdt is deemed to be too large for a stable time step to 
   ! be taken return refineTimeStep = .true.
   subroutine MorphodynamicTimeStepper(RunParams, grid, thisdt, refineTimeStep)
      implicit none

      type(RunSet), intent(inout) :: RunParams
      type(GridType), target, intent(inout) :: grid
      real(kind=wp), intent(inout) :: thisdt
      logical, intent(inout) :: refineTimeStep

      type(TileType), dimension(:), pointer :: tileContainer
      type(TileList), pointer :: activeTiles
      type(TileType), dimension(:), pointer :: intermed0, intermed1, intermed2, intermed3

      integer :: tt, ttk, iw, ib0, ibt, iHnpsi, ipsi, i, j
      real(kind=wp) :: newdt, Hn_old, Hn_new, Hnpsi_old, relhdiff, relhdiffmax
      real(kind=wp) :: deltaBt, psiold, excess_dep, gamold, gamnew, db, w, Hnpsi

      type(RedistList) :: redistcells

      tileContainer => grid%tileContainer
      activeTiles => grid%activeTiles
      intermed0 => grid%intermed0; intermed1 => grid%intermed1
      intermed2 => grid%intermed2; intermed3 => grid%intermed3

      iw = RunParams%Vars%w
      ib0 = RunParams%Vars%b0
      ibt = RunParams%Vars%bt
      ipsi = RunParams%Vars%psi
      iHnpsi = RunParams%Vars%Hnpsi
      refineTimeStep = .false.
      relhdiffmax = 0.0_wp

      ! We compute the initial RHS term in separate ddtExplicit0/ddtExplicitBt0 
      ! arrays. This allows us to roll back to the start of the hydraulic time 
      ! step if necessary.
      call CalculateMorphodynamicRHS(RunParams, grid, intermed0)

      ! Compute first substep for the bed evolution equation.
      do tt = 1, ActiveTiles%Size
         ttk = ActiveTiles%List(tt)

         call ApplySpongeLayer(RunParams, ttk, tileContainer, intermed0, .true.)

         do i = 1, RunParams%nXpertile + 1
            do j = 1, RunParams%nYpertile + 1
               ! We limit bt by the maximum erosion depth here.
               intermed1(ttk)%bt(i, j) = &
                  max(-RunParams%EroDepth, &
                     intermed0(ttk)%bt(i, j) + thisdt * intermed0(ttk)%ddtExplicitBt(i, j))
            end do
         end do
      end do

      ! Update the other variables affected by morphodynamics (w, Hnpsi).
      ! These updates are linearly dependent on the update to bt, so they may
      ! be explicitly determined in a conservative way.
      do tt = 1, ActiveTiles%Size
         ttk = ActiveTiles%List(tt)
         ! First we have to interpolate bt onto cell centres.
         call ComputeCellCentredTopographicData(RunParams, grid, intermed1, ttk)
         do i = 1, RunParams%nXpertile
            do j = 1, RunParams%nYpertile
               gamold = GeometricCorrectionFactor(RunParams, intermed0(ttk)%u(:, i, j))
               gamnew = GeometricCorrectionFactor(RunParams, intermed1(ttk)%u(:, i, j))
               Hn_old = ComputeHn(intermed0(ttk)%u(iw, i, j), &
                  intermed0(ttk)%u(ib0, i, j), intermed0(ttk)%u(ibt, i, j), gamold)
               db = intermed1(ttk)%u(ibt, i, j) - intermed0(ttk)%u(ibt, i, j)
               w = -db / gamnew / gamnew
               w = w  + intermed1(ttk)%u(ibt, i, j)
               w = w + Hn_old * gamold / gamnew / gamnew
               w = w + intermed1(ttk)%u(ib0, i, j)
               intermed1(ttk)%u(iw, i, j) = w
               Hnpsi_old = intermed0(ttk)%u(iHnpsi, i, j)
               Hnpsi = -(1.0_wp - RunParams%BedPorosity) * db / gamnew
               Hnpsi = Hnpsi + Hnpsi_old * gamold / gamnew
               intermed1(ttk)%u(iHnpsi, i, j) = Hnpsi
            end do
         end do
         call ComputeDesingularisedVariables(RunParams, intermed1, ttk, .false.)
      end do

      call CalculateMorphodynamicRHS(RunParams, grid, intermed1)

      ! Compute second substep for the bed evolution equation.
      do tt = 1, ActiveTiles%Size
         ttk = ActiveTiles%List(tt)

         call ApplySpongeLayer(RunParams, ttk, tileContainer, intermed1, .true.)

         ! mom update
         do i = 1, RunParams%nXpertile + 1
            do j = 1, RunParams%nYpertile + 1
               intermed2(ttk)%bt(i, j) = &
                  max(-RunParams%EroDepth, &
                      0.75_wp * intermed0(ttk)%bt(i, j) + &
                      0.25_wp * (intermed1(ttk)%bt(i, j) + &
                      thisdt * intermed1(ttk)%ddtExplicitBt(i, j)))
            end do
         end do
      end do

      ! Update the other variables affected by morphodynamics (w, Hnpsi).
      do tt = 1, ActiveTiles%Size
         ttk = ActiveTiles%List(tt)
         call ComputeCellCentredTopographicData(RunParams, grid, intermed2, ttk)
         do i = 1, RunParams%nXpertile
            do j = 1, RunParams%nYpertile
               gamold = GeometricCorrectionFactor(RunParams, intermed0(ttk)%u(:, i, j))
               gamnew = GeometricCorrectionFactor(RunParams, intermed2(ttk)%u(:, i, j))
               Hn_old = ComputeHn(intermed0(ttk)%u(iw, i, j), &
                  intermed0(ttk)%u(ib0, i, j), intermed0(ttk)%u(ibt, i, j), gamold)
               db = intermed2(ttk)%u(ibt, i, j) - intermed0(ttk)%u(ibt, i, j)
               w = -db / gamnew / gamnew
               w = w + intermed2(ttk)%u(ibt, i, j)
               w = w + Hn_old * gamold / gamnew / gamnew
               w = w + intermed2(ttk)%u(ib0, i, j)
               intermed2(ttk)%u(iw, i, j) = w
               Hnpsi_old = intermed0(ttk)%u(iHnpsi, i, j)
               Hnpsi = -(1.0_wp - RunParams%BedPorosity) * db / gamnew
               Hnpsi = Hnpsi + Hnpsi_old * gamold / gamnew
               intermed2(ttk)%u(iHnpsi, i, j) = Hnpsi
            end do
         end do
         call ComputeDesingularisedVariables(RunParams, intermed2, ttk, .false.)
      end do

      call CalculateMorphodynamicRHS(RunParams, grid, intermed2)

      ! Compute third substep for the bed evolution equation.
      do tt = 1, ActiveTiles%Size
         ttk = ActiveTiles%List(tt)

         call ApplySpongeLayer(RunParams, ttk, tileContainer, intermed2, .true.)

         do i = 1, RunParams%nXpertile + 1
            do j = 1, RunParams%nYpertile + 1
               intermed3(ttk)%bt(i, j) = &
                  max(-RunParams%EroDepth, &
                     (1.0_wp/3.0_wp) * intermed0(ttk)%bt(i,j) + &
                     (2.0_wp/3.0_wp) * (intermed2(ttk)%bt(i,j) + &
                     thisdt * intermed2(ttk)%ddtExplicitBt(i,j)))
            end do
         end do
      end do

      ! Update the other variables affected by morphodynamics (w, Hnpsi).
      do tt = 1, ActiveTiles%Size
         ttk = ActiveTiles%List(tt)
         call ComputeCellCentredTopographicData(RunParams, grid, intermed3, ttk)
         do i = 1, RunParams%nXpertile
            do j = 1, RunParams%nYpertile
               gamold = GeometricCorrectionFactor(RunParams, intermed0(ttk)%u(:, i, j))
               gamnew = GeometricCorrectionFactor(RunParams, intermed3(ttk)%u(:, i, j))
               Hn_old = ComputeHn(intermed0(ttk)%u(iw, i, j), &
                  intermed0(ttk)%u(ib0, i, j), intermed0(ttk)%u(ibt, i, j), gamold)
               db = intermed3(ttk)%u(ibt, i, j) - intermed0(ttk)%u(ibt, i, j)
               w = -db / gamnew / gamnew
               w = w + intermed3(ttk)%u(ibt, i, j)
               w = w + Hn_old * gamold / gamnew / gamnew
               w = w + intermed3(ttk)%u(ib0, i, j)
               intermed3(ttk)%u(iw, i, j) = w
               Hnpsi_old = intermed0(ttk)%u(iHnpsi, i, j)
               Hnpsi = -(1.0_wp - RunParams%BedPorosity) * db / gamnew
               Hnpsi = Hnpsi + Hnpsi_old * gamold / gamnew
               intermed3(ttk)%u(iHnpsi, i, j) = Hnpsi
            end do
         end do
         call ComputeDesingularisedVariables(RunParams, intermed3, ttk, .false.)
      end do

      ! This loop performs two checks on the morphodynamic update.
      ! * Is a cell depositing more solids than it possesses? If so, we
      !   refine the time step or make a correction.
      ! * Is the flow height changing by a large amount? (e.g. 10%)
      !   If so, then time step is likely too large and we refine it.
      !   (N.B. This process can get bogged down if erosionCriticalHeight is
      !   small. In general it should probably be set so it's at least greater
      !   than the solid diameter.)
      do tt = 1, ActiveTiles%Size
         ttk = ActiveTiles%List(tt)
         do j = 1, RunParams%nYpertile
            do i = 1, RunParams%nXpertile
               gamold = GeometricCorrectionFactor(RunParams, intermed0(ttk)%u(:,i,j))
               gamnew = GeometricCorrectionFactor(RunParams, intermed3(ttk)%u(:,i,j))
               Hn_old = ComputeHn(intermed0(ttk)%u(iw, i, j), &
                  intermed0(ttk)%u(ib0, i, j), intermed0(ttk)%u(ibt, i, j), gamold)
               Hn_new = ComputeHn(intermed3(ttk)%u(iw, i, j), &
                  intermed3(ttk)%u(ib0, i, j), intermed3(ttk)%u(ibt, i, j), gamnew)

               ! Return value excess_dep is > 0 if cell deposited more solids 
               ! than it had.
               call ExcessDeposition(RunParams, intermed0, intermed3, &
                                     ttk, i, j, excess_dep, deltaBt, psiold)

               ! n.b. The check for deltaBtn > 0 is in case psi is computed to be slightly
               ! negative, which can be the case for very dilute flows.
               if ((.not. refineTimeStep) .and. excess_dep > epsilon(excess_dep) .and. &
                   deltaBt > epsilon(deltaBt) .and. psiold > -epsilon(psiold)) then
                  ! If the cell is shallow, or only experienced a small 
                  ! amount of deposition, we shall apply a small correction by
                  ! calling RedistributeGrid (below).
                  ! Otherwise, we must refine the time step.
                  if (Hn_old < RunParams%EroCriticalHeight .or. &
                      abs(Hn_new - Hn_old) < RunParams%EroCriticalHeight) then
                     call AddToRedistList(redistcells, ttk, i, j, excess_dep)
                  else
                     refineTimeStep = .true.
                     exit
                  end if
                  cycle
               end if
               ! If the original flow depth is significant, check whether it 
               ! changed by too much (and if so, refine the time step).
               if (Hn_old < RunParams%EroCriticalHeight) cycle
               relhdiff = abs(Hn_new - Hn_old) / abs(Hn_old)
               if (relhdiff > 0.1_wp .and. relhdiff > relhdiffmax) then
                  relhdiffmax = relhdiff
                  refineTimeStep = .true.
                  exit
               end if
            end do
         end do
      end do

      ! Some cells may have been marked for redistribution in the loop above.
      ! This is handled by RedistributeGrid, which makes a small correction to
      ! the morphodynamic update for each cell vertex and its neighbours.
      if (associated(redistcells%head) .and. (.not. refineTimeStep)) then
         call RedistributeGrid(RunParams, grid, redistcells, refineTimeStep)
      end if

      if (refineTimeStep) then
         newdt = 0.5_wp * thisdt
         thisdt = newdt
         call ReportRevisedTimeStep(RunParams, grid%t0, newdt * 0.5_wp)
#if DEBUG_TIMESTEP==1 || DEBUG_TIMESTEP==2
         call InfoMessage('HydraulicTimeStepper refining at morphodynamic step')
#if DEBUG_TIMESTEP==2
         call exit
#endif
#endif
         return
      end if

      ! Update u{Plus,Minus}{X,Y}, b and grad(b) data since the bed has changed.
      do tt = 1, ActiveTiles%Size
         ttk = ActiveTiles%List(tt)
         call ComputeInterfacialTopographicData(RunParams, grid, grid%intermed3, ttk)
      end do

   end subroutine MorphodynamicTimeStepper

   ! Initialise each of the intermediate substep arrays by copying over the
   ! current solution data.
   subroutine InitialiseTimeSteppingArrays(nActiveTiles, grid)
      implicit none

      integer, intent(in) :: nActiveTiles
      type(GridType), target, intent(inout) :: grid

      integer :: tt, ttk

      ! Copy data over to the intermediate time stepping arrays
      do tt = 1, nActiveTiles
         ttk = grid%activeTiles%List(tt)
         grid%intermed0(ttk) = grid%tileContainer(ttk)
         grid%intermed1(ttk) = grid%tileContainer(ttk)
         grid%intermed2(ttk) = grid%tileContainer(ttk)
         grid%intermed3(ttk) = grid%tileContainer(ttk)
      end do
      do tt = 1, grid%ghostTiles%size
         ttk = grid%ghostTiles%List(tt)
         grid%intermed0(ttk) = grid%tileContainer(ttk)
         grid%intermed1(ttk) = grid%tileContainer(ttk)
         grid%intermed2(ttk) = grid%tileContainer(ttk)
         grid%intermed3(ttk) = grid%tileContainer(ttk)
      end do
   end subroutine InitialiseTimeSteppingArrays

   ! Copy the full solution vector from one tile container to another.
   subroutine CopySolutionData(RunParams, grid, tilesfrom, tilesto)
      implicit none

      type(RunSet), intent(in) :: RunParams
      type(GridType), target, intent(inout) :: grid
      type(TileType), dimension(:), intent(inout) :: tilesfrom, tilesto

      integer :: tt, ttk

      type(TileList), pointer :: activeTiles

      activeTiles => grid%activeTiles

      call CopyMutableTopographicData(RunParams, grid, tilesfrom, tilesto)

      do tt = 1, ActiveTiles%Size
         ttk = ActiveTiles%List(tt)
         tilesto(ttk)%u(:, :, :) = tilesfrom(ttk)%u(:, :, :)
         tilesto(ttk)%ddtExplicit(:, :, :) = tilesfrom(ttk)%ddtExplicit(:, :, :)
         tilesto(ttk)%ddtImplicit(:, :, :) = tilesfrom(ttk)%ddtImplicit(:, :, :)
         if (RunParams%MorphodynamicsOn) then
            tilesto(ttk)%ddtExplicitBt(:, :) = tilesfrom(ttk)%ddtExplicitBt(:, :)
         end if
      end do
   end subroutine CopySolutionData

   ! Copy the time dependent topographic data from one tile container to another.
   subroutine CopyMutableTopographicData(RunParams, grid, tilesfrom, tilesto)
      implicit none

      type(RunSet), intent(in) :: RunParams
      type(GridType), target, intent(inout) :: grid
      type(TileType), dimension(:), intent(inout) :: tilesfrom, tilesto

      integer :: tt, ttk

      type(TileList), pointer :: activeTiles

      activeTiles => grid%activeTiles

      do tt = 1, ActiveTiles%Size
         ttk = ActiveTiles%List(tt)
         tilesto(ttk)%u(RunParams%Vars%bt, :, :) = tilesfrom(ttk)%u(RunParams%Vars%bt, :, :)
         tilesto(ttk)%u(RunParams%Vars%dbdx, :, :) = tilesfrom(ttk)%u(RunParams%Vars%dbdx, :, :)
         tilesto(ttk)%bt(:, :) = tilesfrom(ttk)%bt(:, :)

         tilesto(ttk)%uPlusX(RunParams%Vars%bt, :, :) = tilesfrom(ttk)%uPlusX(RunParams%Vars%bt, :, :)
         tilesto(ttk)%uMinusX(RunParams%Vars%bt, :, :) = tilesfrom(ttk)%uMinusX(RunParams%Vars%bt, :, :)
         tilesto(ttk)%uPlusX(RunParams%Vars%dbdx, :, :) = tilesfrom(ttk)%uPlusX(RunParams%Vars%dbdx, :, :)
         tilesto(ttk)%uMinusX(RunParams%Vars%dbdx, :, :) = tilesfrom(ttk)%uMinusX(RunParams%Vars%dbdx, :, :)
         if (.not. RunParams%isOneD) then
            tilesto(ttk)%u(RunParams%Vars%dbdy, :, :) = tilesfrom(ttk)%u(RunParams%Vars%dbdy, :, :)

            tilesto(ttk)%uPlusX(RunParams%Vars%dbdy, :, :) = tilesfrom(ttk)%uPlusX(RunParams%Vars%dbdy, :, :)
            tilesto(ttk)%uMinusX(RunParams%Vars%dbdy, :, :) = tilesfrom(ttk)%uMinusX(RunParams%Vars%dbdy, :, :)

            tilesto(ttk)%uPlusY(RunParams%Vars%bt, :, :) = tilesfrom(ttk)%uPlusY(RunParams%Vars%bt, :, :)
            tilesto(ttk)%uMinusY(RunParams%Vars%bt, :, :) = tilesfrom(ttk)%uMinusY(RunParams%Vars%bt, :, :)
            tilesto(ttk)%uPlusY(RunParams%Vars%dbdx, :, :) = tilesfrom(ttk)%uPlusY(RunParams%Vars%dbdx, :, :)
            tilesto(ttk)%uMinusY(RunParams%Vars%dbdx, :, :) = tilesfrom(ttk)%uMinusY(RunParams%Vars%dbdx, :, :)
            tilesto(ttk)%uPlusY(RunParams%Vars%dbdy, :, :) = tilesfrom(ttk)%uPlusY(RunParams%Vars%dbdy, :, :)
            tilesto(ttk)%uMinusY(RunParams%Vars%dbdy, :, :) = tilesfrom(ttk)%uMinusY(RunParams%Vars%dbdy, :, :)
         end if
      end do
   end subroutine CopyMutableTopographicData

   ! Sponge layer modifications to explicit time stepping terms.
   ! morpho_flag = are we in the morphodynamic operator or not?
   subroutine ApplySpongeLayer(RunParams, ttk, tileContainer, intermed, morpho_flag)
      implicit none

      type(RunSet), intent(in) :: RunParams
      integer, intent(in) :: ttk
      type(tileType), dimension(:), intent(in) :: tileContainer
      type(tileType), dimension(:), intent(inout) :: intermed
      logical, intent(in) :: morpho_flag

      integer :: iw, irhoHnu, irhoHnv, iHnpsi, ibt

      iw = RunParams%Vars%w
      irhoHnu = RunParams%Vars%rhoHnu
      irhoHnv = RunParams%Vars%rhoHnv
      iHnpsi = RunParams%Vars%Hnpsi
      ibt = RunParams%Vars%bt

      if (RunParams%SpongeLayer) then
         if (((.not. RunParams%isOneD) .and. &
              ((tileContainer(ttk)%North == 0) .or. (tileContainer(ttk)%South == 0))) &
             .or. (tileContainer(ttk)%East == 0) .or. (tileContainer(ttk)%West == 0)) then
            if (morpho_flag) then
               intermed(ttk)%ddtExplicit(ibt, :, :) = &
                  intermed(ttk)%ddtExplicit(ibt, :, :) - &
                  RunParams%SpongeStrength * intermed(ttk)%u(ibt, :, :)
            else
               intermed(ttk)%ddtExplicit(iw, :, :) = &
                  intermed(ttk)%ddtExplicit(iw, :, :) - &
                  RunParams%SpongeStrength * intermed(ttk)%u(iw, :, :)
               intermed(ttk)%ddtExplicit(irhoHnu, :, :) = &
                  intermed(ttk)%ddtExplicit(irhoHnu, :, :) - &
                  RunParams%SpongeStrength * intermed(ttk)%u(irhoHnu, :, :)
               intermed(ttk)%ddtExplicit(irhoHnv, :, :) = &
                  intermed(ttk)%ddtExplicit(irhoHnv, :, :) - &
                  RunParams%SpongeStrength * intermed(ttk)%u(irhoHnv, :, :)
               intermed(ttk)%ddtExplicit(iHnpsi, :, :) = &
                  intermed(ttk)%ddtExplicit(iHnpsi, :, :) - &
                  RunParams%SpongeStrength * intermed(ttk)%u(iHnpsi, :, :)
            end if
         end if
      end if
   end subroutine ApplySpongeLayer

   ! Check if the flow as reached any edge of the active region. This is broken
   ! into four separate routines (below) that check each boundary.
   subroutine CheckIfNearBoundaries(RunParams, grid, nActiveTiles)
      implicit none

      type(RunSet), intent(in) :: RunParams
      type(GridType), target, intent(inout) :: grid
      integer, intent(inout) :: nActiveTiles

      type(TileList), pointer :: activeTiles

      activeTiles => grid%activeTiles

      if (.not. RunParams%isOneD) then
         call NearNorthBoundary(RunParams, grid)
         call NearSouthBoundary(RunParams, grid)
      end if
      call NearEastBoundary(RunParams, grid)
      call NearWestBoundary(RunParams, grid)

      ! Return updated value of nActiveTiles.
      if (ActiveTiles%size > nActiveTiles) then
         nActiveTiles = ActiveTiles%size
      end if
   end subroutine CheckIfNearBoundaries

   subroutine NearNorthBoundary(RunParams, grid)

      implicit none

      type(RunSet), intent(in) :: RunParams
      type(GridType), target, intent(inout) :: grid

      integer :: nXtiles, nYtiles
      integer :: nXpertile, nYpertile
      integer :: ii, jj
      integer :: tt, ttk, ttN
      real(kind=wp), dimension(:) :: Hn(RunParams%nYpertile)
      integer(kind=wp) :: ydist

      type(tileType), dimension(:), pointer :: tileContainer
      type(TileList), pointer :: ActiveTiles

      nXtiles = grid%nXtiles
      nYtiles = grid%nYtiles

      nXpertile = RunParams%nXpertile
      nYpertile = RunParams%nYpertile

      tileContainer => grid%tileContainer
      ActiveTiles => grid%ActiveTiles

      do tt = 1, ActiveTiles%Size
         ttk = ActiveTiles%List(tt)
         ttN = tileContainer(ttk)%North 

         ! Check if there exists an inactive neighbour that needs activating.
         if (ttN > 0 .and. (.not. tileContainer(ttN)%TileOn)) then
            ydist = 0
            do ii = 1, nXpertile
               Hn = tileContainer(ttk)%u(RunParams%Vars%Hn, ii, :)
               do jj = nYpertile, 1, -1
                  if (Hn(jj) > RunParams%heightThreshold) then
                     ydist = jj
                     exit
                  end if
               end do
               if (ydist > (RunParams%nYpertile - RunParams%TileBuffer)) then
                  call AddTile(grid, ttN, RunParams)
                  exit
               end if
            end do
         end if
      end do

   end subroutine NearNorthBoundary

   subroutine NearEastBoundary(RunParams, grid)

      implicit none

      type(RunSet), intent(in) :: RunParams
      type(GridType), target, intent(inout) :: grid

      integer :: nXtiles, nYtiles
      integer :: nXpertile, nYpertile
      integer :: ii, jj
      integer :: tt, ttk, ttE
      real(kind=wp), dimension(:) :: Hn(RunParams%nXpertile)
      integer(kind=wp) :: xdist

      type(tileType), dimension(:), pointer :: tileContainer
      type(TileList), pointer :: ActiveTiles

      nXtiles = grid%nXtiles
      nYtiles = grid%nYtiles

      nXpertile = RunParams%nXpertile
      nYpertile = RunParams%nYpertile

      tileContainer => grid%tileContainer
      ActiveTiles => grid%ActiveTiles

      do tt = 1, ActiveTiles%Size
         ttk = ActiveTiles%List(tt)
         ttE = tileContainer(ttk)%East

         ! Check if there exists an inactive neighbour that needs activating.
         if (ttE > 0 .and. (.not. tileContainer(ttE)%TileOn)) then
            xdist = 0
            do jj = 1, nYpertile
               Hn = tileContainer(ttk)%u(RunParams%Vars%Hn, :, jj)
               do ii = nXpertile, 1, -1
                  if (Hn(ii) > RunParams%heightThreshold) then
                     xdist = ii
                     exit
                  end if
               end do
               if (xdist > (RunParams%nXpertile - RunParams%TileBuffer)) then
                  call AddTile(grid, ttE, RunParams)
                  exit
               end if
            end do
         end if
      end do

   end subroutine NearEastBoundary

   subroutine NearSouthBoundary(RunParams, grid)

      implicit none

      type(RunSet), intent(in) :: RunParams
      type(GridType), target, intent(inout) :: grid

      integer :: nXtiles, nYtiles
      integer :: nXpertile, nYpertile
      integer :: ii, jj
      integer :: tt, ttk, ttS
      real(kind=wp), dimension(:) :: Hn(RunParams%nYpertile)
      integer(kind=wp) :: ydist

      type(tileType), dimension(:), pointer :: tileContainer
      type(TileList), pointer :: ActiveTiles

      nXtiles = grid%nXtiles
      nYtiles = grid%nYtiles

      nXpertile = RunParams%nXpertile
      nYpertile = RunParams%nYpertile

      tileContainer => grid%tileContainer
      ActiveTiles => grid%ActiveTiles

      do tt = 1, ActiveTiles%Size
         ttk = ActiveTiles%List(tt)
         ttS = tileContainer(ttk)%South

         ! Check if there exists an inactive neighbour that needs activating.
         if (ttS > 0 .and. (.not. tileContainer(ttS)%TileOn)) then
            ydist = nYpertile
            do ii = 1, nXpertile
               Hn = tileContainer(ttk)%u(RunParams%Vars%Hn, ii, :)
               do jj = 1, nYpertile
                  if (Hn(jj) > RunParams%heightThreshold) then
                     ydist = jj
                     exit
                  end if
               end do
               if (ydist <= RunParams%TileBuffer) then
                  call AddTile(grid, ttS, RunParams)
                  exit
               end if
            end do
         end if
      end do

   end subroutine NearSouthBoundary

   subroutine NearWestBoundary(RunParams, grid)

      implicit none

      type(RunSet), intent(in) :: RunParams
      type(GridType), target, intent(inout) :: grid

      integer :: nXtiles, nYtiles
      integer :: nXpertile, nYpertile
      integer :: ii, jj
      integer :: tt, ttk, ttW
      real(kind=wp), dimension(:) :: Hn(RunParams%nXpertile)
      integer(kind=wp) :: xdist

      type(TileType), dimension(:), pointer :: tileContainer
      type(TileList), pointer :: ActiveTiles

      nXtiles = grid%nXtiles
      nYtiles = grid%nYtiles

      nXpertile = RunParams%nXpertile
      nYpertile = RunParams%nYpertile

      tileContainer => grid%tileContainer
      ActiveTiles => grid%ActiveTiles

      do tt = 1, ActiveTiles%Size
         ttk = ActiveTiles%List(tt)
         ttW = tileContainer(ttk)%West

         ! Check if there exists an inactive neighbour that needs activating.
         if (ttW > 0 .and. (.not. tileContainer(ttW)%TileOn)) then
            xdist = nXpertile
            do jj = 1, nYpertile
               Hn = tileContainer(ttk)%u(RunParams%Vars%Hn, :, jj)
               do ii = 1, nXpertile
                  if (Hn(ii) > RunParams%heightThreshold) then
                     xdist = ii
                     exit
                  end if
               end do
               if (xdist <= RunParams%TileBuffer) then
                  call AddTile(grid, ttW, RunParams)
                  exit
               end if
            end do
         end if
      end do

   end subroutine NearWestBoundary

! The remaining routines below are used for updating data the tracts the maximum
! values of various observables over the length of the simulation. 

   subroutine UpdateMaximumHeights(RunParams, tile, t)

      implicit none

      type(RunSet), intent(in) :: RunParams
      type(tileType), intent(inout) :: tile
      real(kind=wp), intent(in) :: t

      integer :: iHn, ii, jj

      real(kind=wp) :: Hn

      iHn = RunParams%Vars%Hn

      do jj = 1, RunParams%nYpertile
         do ii = 1, RunParams%nXpertile
            Hn = tile%u(iHn, ii, jj)
            if (Hn > RunParams%heightThreshold) then
               if (tile%tfirst(ii, jj) == -1) tile%tfirst(ii, jj) = t
               if (Hn > tile%Hnmax(ii, jj, 1)) then
                  tile%Hnmax(ii, jj, 1) = Hn
                  tile%Hnmax(ii, jj, 2) = t
               end if
            end if
         end do
      end do

   end subroutine UpdateMaximumHeights

   subroutine UpdateMaximumSpeeds(RunParams, tile, t)

      implicit none

      type(RunSet), intent(in) :: RunParams
      type(tileType), intent(inout) :: tile
      real(kind=wp), intent(in) :: t

      real(kind=wp) :: spd

      integer :: iHn, ii, jj

      iHn = RunParams%Vars%Hn

      do jj = 1, RunParams%nYpertile
         do ii = 1, RunParams%nXpertile
            spd = sqrt(FlowSquaredSpeedSlopeAligned(RunParams, tile%u(:, ii, jj)))
#if DEBUG_SPD==1 || DEBUG_SPD==2
            if (spd > 50.0_wp) then
               call InfoMessage('High speed found in UpdateMaximumSpeeds at cell ' // &
                   Int2String(ii) // ',' // Int2String(jj) // ' : spd = ', spd)
#if DEBUG_SPD==2
               call exit(1)
#endif
            end if
#endif
            if ((spd > tile%umax(ii, jj, 1)) .and. &
                (tile%u(iHn, ii, jj) > RunParams%heightThreshold)) then
               tile%umax(ii, jj, 1) = spd
               tile%umax(ii, jj, 2) = t
            end if
         end do
      end do

   end subroutine UpdateMaximumSpeeds

   subroutine UpdateMaximumErosion(RunParams, tile, t)

      implicit none

      type(RunSet), intent(in) :: RunParams
      type(tileType), intent(inout) :: tile
      real(kind=wp), intent(in) :: t

      integer :: ibt, ii, jj

      real(kind=wp) :: bt

      ibt = RunParams%Vars%bt

      do jj = 1, RunParams%nYpertile
         do ii = 1, RunParams%nXpertile
            bt = tile%u(ibt, ii, jj)
            if (bt < 0) then
               if (-bt > tile%emax(ii, jj, 1)) then
                  tile%emax(ii, jj, 1) = -bt
                  tile%emax(ii, jj, 2) = t
               end if
            end if
         end do
      end do

   end subroutine UpdateMaximumErosion

   subroutine UpdateMaximumDeposit(RunParams, tile, t)

      implicit none

      type(RunSet), intent(in) :: RunParams
      type(tileType), intent(inout) :: tile
      real(kind=wp), intent(in) :: t

      integer :: ibt, ii, jj

      real(kind=wp) :: bt

      ibt = RunParams%Vars%bt

      do jj = 1, RunParams%nYpertile
         do ii = 1, RunParams%nXpertile
            bt = tile%u(ibt, ii, jj)
            if (bt > 0) then
               if (bt > tile%dmax(ii, jj, 1)) then
                  tile%dmax(ii, jj, 1) = bt
                  tile%dmax(ii, jj, 2) = t
               end if
            end if
         end do
      end do

   end subroutine UpdateMaximumDeposit

   subroutine UpdateMaximumSolidsFraction(RunParams, tile, t)

      implicit none

      type(RunSet), intent(in) :: RunParams
      type(tileType), intent(inout) :: tile
      real(kind=wp), intent(in) :: t

      integer :: iHn, ipsi, ii, jj

      real(kind=wp) :: psi

      iHn = RunParams%Vars%Hn
      ipsi = RunParams%Vars%psi

      do jj = 1, RunParams%nYpertile
         do ii = 1, RunParams%nXpertile
            psi = tile%u(ipsi, ii, jj)
            if (tile%u(iHn, ii, jj) > RunParams%heightThreshold) then
               if (psi > tile%psimax(ii, jj, 1)) then
                  tile%psimax(ii, jj, 1) = psi
                  tile%psimax(ii, jj, 2) = t
               end if
            end if
         end do
      end do

   end subroutine UpdateMaximumSolidsFraction

end module timestepper_module
