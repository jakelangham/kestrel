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


! Kestrel is a program for simulating shallow Earth-surface flows on complex 
! terrains. It has a modular design that enables it to be specialised to model
! many kinds of phenomena, including floods, debris flows, granular flow. Its 
! underlying governing equations are a formulation of the shallow water
! equations that accommodates topographic variation, a general rheological law
! (via a basal drag term) and morphodynamic transfer between the flow and bed
! (erosion and deposition), documented in https://arxiv.org/abs/2306.16185
!
! Its features include:
!
! * Fully conservative, positivity preserving, well-balanced numerical scheme.
! * Dynamic erosion and deposition of sediment (morphodynamics).
! * Variety of available model closures for basal drag and morphodynamics.
! * Simulation on user-specified topographies, via digital elevation maps (DEMs).
! * Output via geo-referenced NetCDF or plain-text.
!
! Kestrel is free software, provided for research purposes and is intended 
! neither for commercial use, nor for sensitive application such as hazards
! assessment. The terms of its license confer absolutely no warranty.
program Kestrel

   use, intrinsic :: iso_fortran_env, only : stdout=>output_unit, stderr=>error_unit
   use grid_module, only: GridType, InitialiseGrid
   use runsettings_module, only: RunSet, SetImmutableRunSettings
   use read_input_file_module, only: ReadInputFile
   use set_sources_module, only: LoadSourceConditions
   use restart_module, only: LoadInitialCondition
   use dem_module, only: LoadDEM
   use timestepper_module, only: Run
   use output_module, only: CreateOutDir
   use varstring_module, only: varString
   use version_module, only: GetVersion

   implicit none

   ! This data structure stores all the simulation paramters 
   ! (see RunSettings.f90).
   type(RunSet) :: RunParams

   ! Likewise, this contains all the data related to the numerical grid, 
   ! including the solution fields (see Grid.f90).
   type(GridType) :: grid

   character(len=:), allocatable :: version

   ! Set version
   version = GetVersion()
   RunParams%version = varString(trim(version))

   ! The following block indexes the fields that are explicitly stored during 
   ! the simulation. There are four primary flow variables (1-4). 
   !
   ! w = hp + b, where hp is the flow depth, projected on to the vertical and
   ! b is the bed elevation measured vertically
   RunParams%Vars%w = 1
   ! rhoHnu = rho * Hn * u, where rho is the flow density, Hn is the flow depth
   ! measured along the bed normal and u is the velocity in the x direction.
   RunParams%Vars%rhoHnu = 2
   ! rhoHnv = rho * Hn * v, where v is the velocity in the x direction.
   RunParams%Vars%rhoHnv = 3
   ! Hnpsi = Hn * psi, where psi is the volume fraction of the flow occupied by
   ! sediment.
   RunParams%Vars%Hnpsi = 4

   ! Variables (5-9) are derived from the primary variables (1-4), but are 
   ! sufficiently useful that we keep track of them.
   RunParams%Vars%Hn = 5
   RunParams%Vars%u = 6
   RunParams%Vars%v = 7
   RunParams%Vars%psi = 8
   RunParams%Vars%rho = 9

   ! Topographic fields. The bed elevation b is decomposed as b = b0 + bt,
   ! where b0 is the initial elevation and bt records the change in 
   ! elevation as the simulation progresses. We also store the bed gradients
   ! dbdx and dbdy.
   RunParams%Vars%b0 = 10
   RunParams%Vars%bt = 11
   RunParams%Vars%dbdx = 12
   RunParams%Vars%dbdy = 13
   RunParams%Vars%d2bdxx = 14
   RunParams%Vars%d2bdyy = 15
   RunParams%Vars%d2bdxy = 16
   RunParams%Vars%d2bdtx = 17
   RunParams%Vars%d2bdty = 18
   
   RunParams%nDimensions = 18

   write (stdout, *)
   write (stdout, *) "Starting Kestrel (" // RunParams%version%s // ")"
   write (stdout, *)

   ! Read input file and populate the type(RunSet) :: RunParams object
   ! that contains model settings and parameters
   call ReadInputFile(RunParams)
   call SetImmutableRunSettings(RunParams)

   write (stdout, *) "Completed reading of input file.  Continuing to solver."

   ! Create output directory, if needed
   call CreateOutDir(RunParams)

   ! Initialise computational domain
   call InitialiseGrid(RunParams, grid)
   ! Load topographic data, either from DEM raster files or from chosen analytic surface
   call LoadDEM(RunParams, grid)
   ! Load initial condition, if needed
   call LoadInitialCondition(RunParams, grid)
   ! Load source conditions, specified in input file
   call LoadSourceConditions(RunParams, grid)

   ! Commence simulation
   write (stdout, *) "Starting flow simulation."
   write (stdout, *)
   call Run(RunParams, grid)

end program Kestrel
