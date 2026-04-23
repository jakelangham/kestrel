submodule (deposition_types) simple_hindered_settling

! Simple quadratic hindered settling function,
! with zeros at psi = 0 and psi = maximum packing

   use set_precision_module, only: wp
   use dict_module, only: Dict
   use runsettings_module, only: RunSet

   implicit none

contains

   module function SimpleHinderedSettling_constructor(params, RunParams) result(this)
      type(Dict), intent(in) :: params
      type(RunSet), intent(in) :: RunParams
      type(SimpleHinderedSettling) :: this
      call this%initialize(params, RunParams)
      this%name = varString('simple hindered settling')
   end function SimpleHinderedSettling_constructor

   module subroutine SimpleHinderedSettling_validator(this)
      class(SimpleHinderedSettling), intent(in) :: this
      call this%base_validate()
      return
   end subroutine SimpleHinderedSettling_validator

   module pure function SimpleHinderedSettling_hindering(this, psi) result(deposition)
      ! Simple quadratic hindered settling function, with zeros at 
      ! psi = 0 and psi = maximum packing
      implicit none
      class(SimpleHinderedSettling), intent(in) :: this
      real(kind=wp), intent(in) :: psi
      real(kind=wp) :: deposition

      deposition = psi * (1.0_wp - psi / this%maxPack)
   end function SimpleHinderedSettling_hindering

end submodule simple_hindered_settling