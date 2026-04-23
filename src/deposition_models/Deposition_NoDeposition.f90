submodule (deposition_types) no_deposition

   use set_precision_module, only: wp
   use dict_module, only: Dict
   use runsettings_module, only: RunSet

   implicit none

contains

   module function NoDeposition_constructor(params, RunParams) result(this)
      type(Dict), intent(in) :: params
      type(RunSet), intent(in) :: RunParams
      type(NoDeposition) :: this
      
      call this%initialize(params, RunParams)
      this%name = varString('No deposition')
   end function NoDeposition_constructor

   module subroutine NoDeposition_validator(this)
      class(NoDeposition), intent(in) :: this
      call this%base_validate()
   end subroutine NoDeposition_validator

   module pure function NoDeposition_hindering(this, psi) result(flux)
      class(NoDeposition), intent(in) :: this
      real(kind=wp), intent(in) :: psi
      real(kind=wp) :: flux

      flux = 0.0_wp
   end function NoDeposition_hindering

end submodule no_deposition