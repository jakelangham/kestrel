submodule (drag_types) drag_edwards2019

   use set_precision_module, only: wp
   use varstring_module, only: varString
   use dict_module, only: Dict
   use runsettings_module, only: RunSet
   use closures_module, only: FlowSquaredSpeedSlopeAligned, GeometricCorrectionFactor, Reciprocal

   implicit none

   ! -- Drag closure: Edwards 2019 drag --
   !
   !    Extended 'Pouliquen'-like law featuring velocity-weakening behaviour at 
   !    low Froude number.
   !    See Edwards, Russell, Johnson, Gray, J. Fluid Mech. (2019).
   !    N.B. The Fr = 0 is omitted, since stopped regions require special 
   !    attention and this is difficult to implement in Kestrel currently.
   !
   !    
   !    Parameters:
   !        pouliquen_min > 0
   !        pouliquen_max > 0
   !        pouliquen_intermediate > 0
   !        pouliquen_beta > 0
   !        edwards2019_betastar > 0
   !        edwards2019_kappa > 0
   !        edwards2019_gamma > 0
   !
contains

   ! Constructor
   module function edwards2019_constructor(params) result(this)
      type(Dict), intent(in) :: params
      type(Edwards2019Drag) :: this
      this%name = varString("Edwards2019")
      if (params%has_key('pouliquen min')) then
         call params%get('pouliquen min', this%pouliquen_min)
      else
         call WarningMessage('Pouliquen Min not given.  Using default value ', this%pouliquen_min)
      end if
      if (params%has_key('pouliquen max')) then
         call params%get('pouliquen max', this%pouliquen_max)
      else
         call WarningMessage('Pouliquen Max not given.  Using default value ', this%pouliquen_max)
      end if
      if (params%has_key('pouliquen intermediate')) then
         call params%get('pouliquen intermediate', this%pouliquen_intermediate)
      else
         call WarningMessage('Pouliquen Intermediate not given.  Using default value ', this%pouliquen_intermediate)
      end if
      if (params%has_key('pouliquen beta')) then
         call params%get('pouliquen beta', this%pouliquen_beta)
      else
         call WarningMessage('Pouliquen Beta not given.  Using default value ', this%pouliquen_beta)
      end if
      if (params%has_key('edwards2019 betastar')) then
         call params%get('edwards2019 betastar', this%edwards2019_betastar)
      else
         call WarningMessage('Edwards2019 betastar not given.  Using default value ', this%edwards2019_betastar)
      end if
      if (params%has_key('edwards2019 kappa')) then
         call params%get('edwards2019 kappa', this%edwards2019_kappa)
      else
         call WarningMessage('Edwards2019 kappa not given.  Using default value ', this%edwards2019_kappa)
      end if
      if (params%has_key('edwards2019 gamma')) then
         call params%get('edwards2019 gamma', this%edwards2019_gamma)
      else
         call WarningMessage('Edwards2019 gamma not given.  Using default value ', this%edwards2019_gamma)
      end if
      if (params%has_key('solid diameter')) then
         call params%get('solid diameter', this%solid_diameter)
      else
         call FatalErrorMessage('Required parameters Solid Diameter not given.')
      end if

      ! Parameters
      call this%parameters%append("pouliquen min", this%pouliquen_min)
      call this%parameters%append("pouliquen max", this%pouliquen_max)
      call this%parameters%append("pouliquen intermediate", this%pouliquen_intermediate)
      call this%parameters%append("pouliquen beta", this%pouliquen_beta)
      call this%parameters%append("edwards2019 betastar", this%edwards2019_betastar)
      call this%parameters%append("edwards2019 kappa", this%edwards2019_kappa)
      call this%parameters%append("edwards2019 gamma", this%edwards2019_gamma)
   end function edwards2019_constructor

   ! Validator
   module subroutine edwards2019_validate(this)
      class(Edwards2019Drag), intent(in) :: this

      if (this%pouliquen_min<=0) then
         call FatalErrorMessage('Pouliquen Min must be positive.  Received', this%pouliquen_min)
      end if
      if (this%pouliquen_max<=0) then
         call FatalErrorMessage('Pouliquen max must be positive.  Received', this%pouliquen_max)
      end if
      if (this%pouliquen_beta<=0) then
         call FatalErrorMessage('Pouliquen beta must be positive.  Received', this%pouliquen_beta)
      end if
      if (this%edwards2019_betastar<=0) then
         call FatalErrorMessage('Edwards2019 betastar must be positive.  Received', this%edwards2019_betastar)
      end if
      if (this%edwards2019_kappa<=0) then
         call FatalErrorMessage('Edwards2019 kappa must be positive.  Received', this%edwards2019_kappa)
      end if
      if (this%edwards2019_gamma<0) then
         call FatalErrorMessage('Edwards2019 gamma must be non-negative.  Received', this%edwards2019_gamma)
      end if
      if (this%solid_diameter<=0) then
         call FatalErrorMessage('Solid Diameter must be positive.  Received', this%solid_diameter)
      end if
   end subroutine edwards2019_validate

   ! Friction
   module pure function edwards2019_friction(this, RunParams, uvect) result(friction)
      class(Edwards2019Drag), intent(in) :: this
      type(RunSet), intent(in) :: RunParams
      real(kind=wp), dimension(:), intent(in) :: uvect
      real(kind=wp) :: friction

      real(kind=wp) :: modu2
      real(kind=wp) :: Hn, gam, gperp, mu

      Hn = uvect(RunParams%Vars%Hn)
      gam = GeometricCorrectionFactor(RunParams, uvect)
      gperp = RunParams%g / gam

      modu2 = FlowSquaredSpeedSlopeAligned(RunParams, uvect)

      mu = this%friction_coefficient(RunParams%heightThreshold, gperp, Hn, sqrt(modu2))

      friction = mu * gperp * Hn

   end function edwards2019_friction

   ! Friction coefficient
   module pure function edwards2019_friction_coefficient(this, heightThreshold, gcostheta, Hn, modu) result(mu)
      class(Edwards2019Drag), intent(in) :: this
      real(kind=wp), intent(in) :: heightThreshold
      real(kind=wp), intent(in) :: gcostheta
      real(kind=wp), intent(in) :: Hn
      real(kind=wp), intent(in) :: modu
      real(kind=wp) :: mu

      real(kind=wp) :: betastar, beta, mu1, mu2, mu3, capgam, L, kappa
      real(kind=wp) :: Fr, Hn_recip

      mu1 = this%pouliquen_min
      mu2 = this%pouliquen_max
      mu3 = this%pouliquen_intermediate
      beta = this%pouliquen_beta
      betastar = this%edwards2019_betastar
      kappa = this%edwards2019_kappa
      capgam = this%edwards2019_gamma
      L = this%solid_diameter

      Hn_recip = Reciprocal(Hn, heightThreshold)
      Fr = modu * sqrt(Hn_recip / gcostheta)

      if (Fr > betastar) then
         mu = mu1 + (mu2 - mu1) / (1.0_wp + Hn * beta / (L * (Fr + capgam)))
      else
         mu = ((Fr / betastar)**kappa) * &
             (mu1 + (mu2 - mu1) / (1.0_wp + Hn * beta / (L * (betastar + capgam))) - &
             mu3 - (mu2 - mu1) / (1.0_wp + Hn / L)) + &
             mu3 + (mu2 - mu1) / (1.0_wp + Hn / L)
      end if
   end function edwards2019_friction_coefficient

end submodule drag_edwards2019