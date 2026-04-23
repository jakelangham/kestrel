submodule (deposition_types) spearman_manning_hindered_settling

! Hindered settling function from Spearman & Manning (2017)
! doi:10.1007/s10236-017-1034-7.

   use set_precision_module, only: wp
   use dict_module, only: Dict
   use runsettings_module, only: RunSet

   implicit none

contains

   module function SpearmanManningHinderedSettling_constructor(params, RunParams) result(this)
      type(Dict), intent(in) :: params
      type(RunSet), intent(in) :: RunParams
      type(SpearmanManningHinderedSettling) :: this
      real(kind=wp) :: Rep

      call this%initialize(params, RunParams)
      this%name = varString('Spearman Manning hindered settling')
      
      Rep = sqrt(this%g*this%solid_diameter)*this%solid_diameter/this%visc_w ! Particle Reynolds number

      ! Exponent in hindered settling form
      ! (Rowe 1987; A convenient empirical equation for estimation of the Richardson-Zaki exponent)
      this%nsettling = (4.7_wp+0.41_wp*Rep**(0.75_wp))/(1.0_wp + 0.175_wp*Rep**(0.75_wp))

      this%a = 2.7_wp - 0.15_wp * this%nsettling
      this%b = 0.62_wp * this%nsettling - 1.46_wp
   end function SpearmanManningHinderedSettling_constructor

   module subroutine SpearmanManningHinderedSettling_validator(this)
      class(SpearmanManningHinderedSettling), intent(in) :: this
      call this%base_validate()
      return
   end subroutine SpearmanManningHinderedSettling_validator

   module pure function SpearmanManningHinderedSettling_hindering(this, psi) result(deposition)
      ! Hindered settling function from Spearman & Manning (2017)
      ! doi:10.1007/s10236-017-1034-7.
      implicit none
      class(SpearmanManningHinderedSettling), intent(in) :: this
      real(kind=wp), intent(in) :: psi
      real(kind=wp) :: deposition

      deposition = psi * (1.0_wp - psi)**this%a * (1.0_wp - psi / this%maxPack)**this%b
   end function SpearmanManningHinderedSettling_hindering

end submodule spearman_manning_hindered_settling