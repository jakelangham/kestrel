module blur_module

   use set_precision_module, only: wp, pi

   implicit none

   public :: GaussBlur
   public :: BoxBlur
   public :: Binom3Blur
   public :: MultiBoxBlur
   public :: MultiBinom3Blur
   private :: BoxBlur_1d, BoxBlur_2d
   private :: MultiBoxBlur_1d, MultiBoxBlur_2d
   private :: MultiBinom3Blur_1d, MultiBinom3Blur_2d
   private :: Binom3BlurH

   interface BoxBlur
      module procedure :: BoxBlur_1d, BoxBlur_2d
   end interface

   interface Binom3Blur
      module procedure :: Binom3Blur_1d, Binom3Blur_2d
   end interface

   interface MultiBoxBlur
      module procedure :: MultiBoxBlur_1d, MultiBoxBlur_2d
   end interface

   interface MultiBinom3Blur
      module procedure :: MultiBinom3Blur_1d, MultiBinom3Blur_2d
   end interface

contains

   function GaussBlur(source, r) result(filtered)
      implicit none
      real(kind=wp), dimension(:,:), intent(in) :: source
      integer, intent(in) :: r
      real(kind=wp), dimension(:,:) :: filtered(size(source,1),size(source,2))

      integer :: i, j, k, l
      integer :: ii, jj      ! inside the array
      integer :: nx, ny
      integer :: dim(2)
      integer :: r2
      real(kind=wp) :: val, wsum
      real(kind=wp) :: dsq         ! distance squared
      real(kind=wp) :: wght        ! weight

      dim = shape(source)
      nx = dim(1)
      ny = dim(2)

      r2 = r*r

      do i = 1, nx
         do j = 1, ny
            val = 0
            wsum = 0
            do k = i-r, i+r     ! inner loop over kernel
               do l = j-r, j+r  ! inner loop over kernel
                  ii = min(nx, max(1,k))   ! make sure i is always inside the grid (this implies values are extendet (stretched at the boundaries))
                  jj = min(ny, max(1,l))   ! make sure j is always inside the grid (this implies values are extendet (stretched at the boundaries))
                  dsq = (j-l)**2 + (i-k)**2
                  wght = exp(-dsq / (2.0_wp*r2)) / (2.0_wp*pi*r2)
                  val = val + source(ii,jj) * wght
                  wsum = wsum + wght
               end do
            end do
            filtered(i,j) = val / wsum
         end do
      end do

      return
   end function GaussBlur

   pure function BoxBlurH(source, r) result(filtered)
      ! computes horizontal blur
      implicit none
      real(kind=wp), dimension(:,:), intent(in) :: source
      integer, intent(in) :: r
      real(kind=wp), dimension(:,:) :: filtered(size(source,1),size(source,2))
  
      integer :: nx, ny
      integer :: dim(2)
  
      real(kind=wp) :: wght  ! weight
      real(kind=wp) :: sum  
      integer :: i,j
      integer :: il    ! leftmost  pixel which should be removed from the accumulator
      integer :: ir    ! rightmost pixel which should be added   to   the accumulator
  
      dim = shape(source)
      nx = dim(1)
      ny = dim(2)
  
      wght = 1.0_wp / (2*r + 1)
  
      do j = 1, ny   ! loop over all rows
         ! compute sum at first pixel
         sum = source(1,j)
         do i = 1, r  ! loop over box kernel
            sum = sum + source(1,j) + source(i+1,j) ! always take 1 as left pixel, to not get out of grid
         end do
  
         ! generate output pixel, then update running sum
         do i = 1, nx
            filtered(i,j) = sum * wght
            il = max(i-r, 1)     ! make sure we dont get off the grid
            ir = min(i+r+1, nx)  ! make sure we dont get off the grid
            sum = sum + source(ir,j) - source(il,j)
         end do
      end do
      return         
   end function BoxBlurH

   pure function BoxBlur_1d(source, r) result(filtered)
      ! computes 1d box blur
      implicit none
      real(kind=wp), dimension(:), intent(in) :: source
      integer, intent(in) :: r
      real(kind=wp), dimension(:) :: filtered(size(source))
      
      integer :: nx
      real(kind=wp) :: wght  ! weight
      real(kind=wp) :: sum
      integer :: i
      integer :: il    ! leftmost  pixel which should be removed from the accumulator
      integer :: ir    ! rightmost pixel which should be added   to   the accumulator

      nx = size(source)

      wght = 1.0_wp / (2*r + 1)

      sum = source(1)
      do i = 1, r  ! loop over box kernel
         sum = sum + source(1) + source(i+1) ! always take 1 as left pixel, to not get out of grid
      end do
  
      ! generate output pixel, then update running sum
      do i = 1, nx
         filtered(i) = sum * wght
         il = max(i-r, 1)     ! make sure we dont get off the grid
         ir = min(i+r+1, nx)  ! make sure we dont get off the grid
         sum = sum + source(ir) - source(il)
      end do

      return
   end function BoxBlur_1d

   pure function BoxBlur_2d(source, r) result(filtered)
      ! computes 2d box blur, by calling horizontal boxblur twice, the second time with tansposed matrix
      implicit none
      real(kind=wp), dimension(:,:), intent(in) :: source
      integer, intent(in) :: r
      real(kind=wp), dimension(:,:) :: filtered(size(source,1),size(source,2))
      
      integer :: nx, ny
      integer :: dim(2)
      
      real(kind=wp), dimension(:,:), allocatable :: source_t  ! source transposed
      real(kind=wp), dimension(:,:), allocatable :: filtered_t  ! filtered transposed
      
      dim = shape(source)
      nx = dim(1)
      ny = dim(2)
      allocate(source_t(ny,nx))
      allocate(filtered_t(ny,nx))
  
      ! first horizontal blur
      filtered = BoxBlurH(source, r)
      if (ny>1) then
         ! then transpose result and call horizontal blur again, now really doing vertical blur
         source_t = transpose(filtered)
         filtered_t = BoxBlurH(source_t, r)
         ! transpose again, to get back initial shape
         filtered = transpose(filtered_t)
      end if

      return
   end function BoxBlur_2d

   pure function Binom3Blur_1d(source) result(filtered)
      ! computes 1d 3-tap binomial blur
      implicit none
      real(kind=wp), dimension(:), intent(in) :: source
      real(kind=wp), dimension(:) :: filtered(size(source))

      integer :: N

      N = size(source)

      filtered(1) = source(1)
      filtered(2:N-1) = (2*source(2:N-1) + source(1:N-2) + source(3:N)) / 4
      filtered(N) = source(N)

      return

   end function Binom3Blur_1d

   pure function Binom3BlurH(source) result(filtered)
      ! computes 2d 3-tap binomial blur
      implicit none
      real(kind=wp), dimension(:,:), intent(in) :: source
      real(kind=wp), dimension(:,:) :: filtered(size(source,1), size(source,2))

      integer :: M

      M = size(source,1)

      filtered(1,:) = source(1,:)
      filtered(2:M-1,:) = (2*source(2:M-1,:) + source(1:M-2,:) + source(3:M,:)) / 4
      filtered(M,:) = source(M,:)

      return
   end function Binom3BlurH

   pure function Binom3Blur_2d(source) result(filtered)
      ! computes 2d 3-tap binomial blur
      implicit none
      real(kind=wp), dimension(:,:), intent(in) :: source
      real(kind=wp), dimension(:,:) :: filtered(size(source,1), size(source,2))
      real(kind=wp), dimension(:,:) :: filtered_t(size(source,2), size(source,1))

      filtered = Binom3BlurH(source)
      
      filtered_t = transpose(filtered)
      filtered_t = Binom3BlurH(filtered_t)
      filtered = transpose(filtered_t)

      return

   end function Binom3Blur_2d


   pure function MultiBoxBlur_1d(source, r, n) result(filtered)
      ! performs multiple box blurs to mimic a Guassian blur
      implicit none
      real(kind=wp), dimension(:), intent(in) :: source
      integer, intent(in) :: r
      integer, intent(in) :: n
      real(kind=wp), dimension(:) :: filtered(size(source))

      integer :: i

      filtered = source
      do i=1,n
         filtered = BoxBlur_1d(filtered, r)
      end do

      return

   end function MultiBoxBlur_1d

   pure function MultiBoxBlur_2d(source, r, n) result(filtered)
      ! performs multiple box blurs to mimic a Guassian blur
      implicit none
      real(kind=wp), dimension(:,:), intent(in) :: source
      integer, intent(in) :: r
      integer, intent(in) :: n
      real(kind=wp), dimension(:,:) :: filtered(size(source,1),size(source,2))

      integer :: i

      filtered = source
      do i=1,n
         filtered = BoxBlur_2d(filtered, r)
      end do

      return

   end function MultiBoxBlur_2d

   pure function MultiBinom3Blur_1d(source, n) result(filtered)
      implicit none
      real(kind=wp), dimension(:), intent(in) :: source
      integer, intent(in) :: n
      real(kind=wp), dimension(:) :: filtered(size(source))

      integer :: i

      filtered = source
      do i=1,n
         filtered = Binom3Blur_1d(filtered)
      end do

      return

   end function MultiBinom3Blur_1d

   pure function MultiBinom3Blur_2d(source, n) result(filtered)
      implicit none
      real(kind=wp), dimension(:,:), intent(in) :: source
      integer, intent(in) :: n
      real(kind=wp), dimension(:,:) :: filtered(size(source,1),size(source,2))

      integer :: i

      filtered = source
      do i=1,n
         filtered = Binom3Blur_2d(filtered)
      end do

      return

   end function MultiBinom3Blur_2d

end module blur_module