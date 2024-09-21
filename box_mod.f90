#include "error.inc"
!----------------------------------------------------------------------------
!>
!> @author G.S.Fanourgakis
!>
!> @brief  
!>
!> @details 
!>
!> @todo 
!> 
!----------------------------------------------------------------------------
module box_mod
use errors_mod,    only: handle_error
implicit none

type t_box
    double precision, dimension(3) :: l    ! box-length
    double precision, dimension(3) :: a    ! angles
    double precision, dimension(3) :: cc   ! 3 dimensions based on which nx, ny nz are determined
    double precision, dimension(3, 3) :: rprimd, gprimd
    double precision               :: volume    ! 
    double precision, dimension(3) :: inv    ! 
    logical                        :: pbc
end type t_box

private
public :: init_box, update_box, calcRij, cart2frac, frac2cart, getvolume, t_box, rprim2la, la2rprim
public :: calcRijORTHO

interface cart2frac
   module procedure cart2frac1
   module procedure cart2fracN
end interface cart2frac

interface frac2cart
   module procedure frac2cart1
   module procedure frac2cartN
end interface frac2cart

contains
   !
   !
   !   subroutine init_box( box )
   !   subroutine la2rprim( box )
   !   subroutine rprim2la( box )
   !   subroutine update_box( box )
   !   subroutine calcRij(box, rij, drsq)
   !   subroutine cart2frac1(box, r1, r2)
   !   subroutine frac2cart1(box, r1, r2)
   !   subroutine cart2fracN(box, r1, r2)
   !   subroutine frac2cartN(box, r1, r2)
   !   subroutine calcRijORTHO(box, r, drsq)
   !
   subroutine init_box( box )
   type (t_box),    intent( inout ) :: box

   if (box% l(1) > 0.) then
      box% pbc =.true.
      if (box% l(2) <=0.d0 ) box% l(2) = box% l(1) 
      if (box% l(3) <=0.d0 ) box% l(3) = box% l(2) 

      call la2rprim (  box  )
      call update_box (  box  )

   else
      box% pbc=.false.
   endif
   end subroutine init_box
   !
   !
   !
   subroutine la2rprim( box )
   ! having rprimd  get boxlength (l) and angles (a)
   use consts,   only: PI
   type (t_box),    intent( inout ) :: box
   !... local variables
   double precision :: cos_a, cos_b, cos_c, sin_c
   double precision :: term1, term2, term3

   cos_a = cos(box% a(1))
   cos_b = cos(box% a(2))
   cos_c = cos(box% a(3))
   sin_c = sin(box% a(3))

   box% rprimd(1, 1:3) = [ box% l(1), 0.d0, 0.d0 ]
   box% rprimd(2, 1:3) = [ box% l(2)*cos_c, box% l(2)*sin_c, 0.d0 ]
   term1 = box% l(3)*cos_b
   term2 = box% l(3) * (cos_a - cos_c*cos_b) / sin_c
   term3 = box% l(3) * sqrt(1.d0 -cos_a**2 - cos_b**2 - cos_c**2 + 2.d0*cos_a*cos_b*cos_c) / sin_c
   box% rprimd(3, 1:3) = [ term1, term2, term3 ]
   box% rprimd = transpose(box% rprimd)
   end subroutine la2rprim
   !
   !
   !
   subroutine rprim2la( box )
   ! having rprimd  get boxlength (l) and angles (a)
   use consts,   only: PI
   type (t_box),    intent( inout ) :: box
   !... local variables
   double precision :: cos_a, cos_b, cos_c, sin_c

   box% l(1) = box% rprimd(1,1)
   box% l(2) = sqrt(box% rprimd(1,2)**2 + box% rprimd(2,2)**2)
   box% l(3) = sqrt(box% rprimd(1,3)**2 + box% rprimd(2,3)**2 + box% rprimd(3, 3)**2)
   cos_c = box% rprimd(1,2)  /  box% l(2)
   sin_c = box% rprimd(2,2)  /  box% l(2)
   cos_b = box% rprimd(1,3)  /  box% l(3)
   cos_a = box% rprimd(2, 3)*sin_c/box% l(3) + cos_b*cos_c
   box% a(3) = acos(cos_c)
   box% a(2) = acos(cos_b)
   box% a(1) = acos(cos_a)
   end subroutine rprim2la
   !
   !
   ! starts from box% rprimd
   subroutine update_box( box )
   ! @brief  
   !  on input:
   !  * pbc
   !  * rprimd  (use la2rprim to compute it)
   ! on output
   ! * gprimd
   ! * volume
   use math_mod, only: matr3inv
   implicit none
   type (t_box),    intent( inout ) :: box

   if (box% pbc) then
      call matr3inv(box% rprimd, box% gprimd)
      !SOS
      box% gprimd = transpose(box% gprimd)
      box% volume = getvolume( box )
   endif
   end subroutine update_box
   !
   !
   subroutine calcRij(box, rij, drsq, rf)
   !... on input rij is at fractional coordinates
   type (t_box),                   intent( in    ) :: box
   double precision, dimension(3), intent( inout ) :: rij
   double precision, optional,     intent(   out ) :: drsq
   double precision, dimension(3), optional,  intent(   out ) :: rf
   !... local variables
   integer :: ii
   double precision, dimension(3) :: r1

   if (box% pbc) then
      do ii=1, 3
         r1(ii) = Rij(ii) + 1.d0
         do while (r1(ii) .gt. 0.5d0)
            r1(ii) = r1(ii) - 1.d0
         end do
         do while (r1(ii) .lt. -0.5d0)
            r1(ii) = r1(ii) + 1.d0
         end do
      enddo
      if (present(rf)) rf=r1
      call frac2cart(box, r1, rij)
   endif
!   rij = matmul(box% rprimd, rij)

   if (present(drsq)) drsq=rij(1)*rij(1) + rij(2)*rij(2) + rij(3)*rij(3)
   end subroutine calcRij

   !--------------------------------------------------------------------------
   !>  @author  GSF
   !!  @brief if (pbc) transforms the cartesian (r1) to fractional (r2) coordinates
   !!  @param[in] r1
   !!  @param[out] r2
   !--------------------------------------------------------------------------
   subroutine cart2frac1(box, r1, r2)
   use consts, only: PI
   type (t_box),                   intent( in    ) :: box
   double precision, dimension(3), intent( in    ) :: r1
   double precision, dimension(3), intent(   out ) :: r2

   if (box% pbc) then
      r2 = matmul(box% gprimd, r1)
   else
      r2 = r1
   endif
   end subroutine cart2frac1
   !--------------------------------------------------------------------------
   !>  @author  GSF
   !!  @brief if (pbc) transforms the fractional (r1) to cartesian (r2) coordinates
   !!  @param[in] r1
   !!  @param[out] r2
   !--------------------------------------------------------------------------
   subroutine frac2cart1(box, r1, r2)
   implicit none
   type (t_box),                   intent( in    ) :: box
   double precision, dimension(3), intent( in    ) :: r1
   double precision, dimension(3), intent(   out ) :: r2

   if (box% pbc) then
!      r2(1) = box% rprimd(1,1)*r1(1) + box% rprimd(1,2)*r1(2) + box% rprimd(1,3)*r1(3)
!      r2(2) =                          box% rprimd(2,2)*r1(2) + box% rprimd(2,3)*r1(3)
!     r2(3) =                                                   box% rprimd(3,3)*r1(3)
      r2 = matmul(box% rprimd, r1)
   else
        r2 = r1
   endif
   end subroutine frac2cart1

   !--------------------------------------------------------------------------
   !>  @author  GSF
   !!  @brief if (pbc) transforms the cartesian (r1) to fractional (r2) coordinates
   !!  @param[in] r1(3,*)
   !!  @param[out] r2(3,*)
   !--------------------------------------------------------------------------
   subroutine cart2fracN(box, r1, r2)
   use consts, only: PI
   type (t_box),                     intent( in    ) :: box
   double precision, dimension(:,:), intent( in    ) :: r1
   double precision, dimension(:,:), intent(   out ) :: r2
   !.... local variables
   integer :: i, N

   n = size(r1, dim=2)
   if (n /= size(r2, dim=2)) call handle_error('wrong dimensions ', __FILE__, __LINE__)
   if (box% pbc) then
      do i=1, N
         r2(:,i) = matmul(box% gprimd, r1(:,i))
      enddo
   else
      r2 = r1
   endif
   end subroutine cart2fracN
   !--------------------------------------------------------------------------
   !>  @author  GSF
   !!  @brief if (pbc) transforms the fractional (r1) to cartesian (r2) coordinates
   !!  @param[in] r1
   !!  @param[out] r2
   !--------------------------------------------------------------------------
   subroutine frac2cartN(box, r1, r2)
   implicit none
   type (t_box),                   intent( in    ) :: box
   double precision, dimension(:, :), intent( in    ) :: r1
   double precision, dimension(:, :), intent(   out ) :: r2
   !.... local variables
   integer :: i, N

   n = size(r1, dim=2)
   if (n /= size(r2, dim=2)) call handle_error('wrong dimensions ', __FILE__, __LINE__)
   if (box% pbc) then
      do i=1, n
         r2(:, i) = matmul(box% rprimd, r1(:, i))
      enddo
   else
        r2 = r1
   endif
   end subroutine frac2cartN
   !--------------------------------------------------------------------------
   !>  @author  GSF
   !!  @brief if (pbc) gets box volume
   !--------------------------------------------------------------------------
   function getvolume(box) result (vol)
   implicit none
   type (t_box),                   intent( in    ) :: box
   double precision                                :: vol

   if (box% pbc) then
!      vol = box% l(1) * box% l(2) * box% l(3)
      vol = box% rprimd(1,1) * box% rprimd(2,2) * box% rprimd(3,3)
   endif
   end function getvolume
   !
   !
   !
   subroutine calcRijORTHO(box, r, drsq)
   implicit none
   type(t_box), intent( in    ) :: box
   double precision, dimension(3), intent( inout) :: r
   double precision, intent(   out), optional :: drsq
   !...
   r = r - box% l*anint(r/box% l)
   if (present(drsq)) drsq = r(1)*r(1) + r(2)*r(2) + r(3)*r(3)

   end subroutine calcRijORTHO


end module box_mod
