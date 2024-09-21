!------------------------------------------------------------------------------
!> MODULE bond_mod
!------------------------------------------------------------------------------
!>
!> @author G.S.Fanourgakis
!>
!> @brief handles intramolecular bonded interactions
!>
!> @todo add MORSE and other pair interactions
!> 
!>-----------------------------------------------------------------------------
module bond_mod
use errors_mod,         only: handle_error
use t_intramolecular_mod,   only: t_bondtype, t_bond
implicit none

PRIVATE
PUBLIC :: calc_bonds, bond_harmonic

Contains
   !----------------------------------------------------------------------------
   !>
   !> @author G.S.Fanourgakis
   !>
   !> @brief 
   !>
   !> @todo  correct harmonic/virial
   !> 
   !----------------------------------------------------------------------------
   subroutine calc_bonds(bonds, box, R, e, d, v, ift, ilt)
   use time_mod,    only: timer_start, timer_stop, itm_bonds
   use box_mod,       only: calcrij, t_box
   implicit none
   !...
   type( t_box ),                       intent( in    ) :: box
   type (t_bond),                       intent( in    ) :: bonds
   double precision, dimension(:, :),   intent( in    ) :: R
   double precision,                    intent(   out ) :: e
   double precision, dimension(:, :),   intent( inout ) :: d
   double precision, dimension(:, :),   intent(   out ) :: v
   integer,                             intent( in    ), optional :: ift, ilt
   !... local variables
   integer :: it, ifirst, ilast
   double precision, dimension(3, 3) :: v1
   double precision :: e1

   call timer_start(itm_bonds)
   ifirst = 1;
   ilast = bonds%ntypes
   if (present(ift)) ifirst=ift
   if (present(ilt)) ilast =ilt

   e=0.d0
   v=0.d0

   !SOS
   do it=ifirst,ilast
!   do it=1,1

      select case ( bonds%type(it)%type )
         ! 
         case ( 'NONE' ) 
         case ( 'HARMONIC', 'harmonic' ) 
            call bond_harmonic(bonds%type(it), box, r, e1, d, v1)
         !
!         case ( 'MORSE' ) 
         case default
            call handle_error('Unknown type of BOND potential', __FILE__, __LINE__)
      end select
      e = e + e1
      v = v + v1

   enddo   !  do it=ifirst,ilast
!   print*,'e_bonds=', e
!   print*,'vir_bonds=',char(10), v


   call timer_stop(itm_bonds)

   end subroutine calc_bonds
   !----------------------------------------------------------------------------
   !>
   !> @author G.S.Fanourgakis
   !>
   !> @brief   Follows several forms of the potential function
   !>  * harmonic: par=2  (k, r0)
   !>                  V=k*(r-r0)**2
   !>  * morse:    par=3  (De, a, r0)
   !>                  V = De*(1 - exp(-a*(r-r0))) 
   !> @todo  implement, morse, morseaprx (Manolopoulos), unharmonic
   !>             (Have a look also on lammps)
   !> 
   !----------------------------------------------------------------------------
   !
   !    HARMONIC
   !
   !----------------------------------------------------------------------------
   subroutine bond_harmonic(bondI, box, r, e, d, v, dedp, dvdp )
   use box_mod,       only: calcrij, t_box
   type (t_box),                             intent( in    ) :: box
   type (t_bondtype),                        intent( in    ) :: bondI
   double precision, dimension(:,:),         intent( in    ) :: r
   double precision,                         intent(   out ) :: e
   double precision, dimension(:,:),         intent( inout ) :: d
   double precision, dimension(:,:),         intent(   out ) :: v
   double precision, dimension(2), optional, intent( inout ) :: dedp
   double precision, dimension(2), optional, intent( inout ) :: dvdp
   !... local variables
   integer :: ia, ja, ibnd
   double precision :: k, r0, dr, drsq, tmp
   double precision, dimension(3) :: Rij, di

   e=0.d0
   v=0.d0

   k  = bondI%params(1) 
   r0 = bondI%params(2) 
   !... NO derivatives wrt potential parameters 
   if (.not. (present(dedp) .and. present(dvdp))) then
      do ibnd=1, bondI%N
         ia = bondI%ia(ibnd)
         ja = bondI%ja(ibnd)
         Rij=R(1:3,ia) - R(1:3,ja)
         call calcRij(box, rij, drsq)
         dr = sqrt(drsq)
         e = e + k*(dr-r0)**2
         tmp = 2. * k * (dr-r0)/dr
         di = tmp*rij
         d(1:3,ia) = d(1:3,ia) + di
         d(1:3,ja) = d(1:3,ja) - di
         v = v + reshape( (/rij*di(1), rij*di(2), rij*di(3)/), (/3,3/) )
         !SOS
!         print*,'iiii=', ia, ja, dr !k*(dr-r0)**2

      enddo
   !... derivatives wrt potential parameters 
   else
      dedp = 0.d0
      dvdp = 0.d0

      do ibnd=1, bondI%N
         ia = bondI%ia(ibnd)
         ja = bondI%ja(ibnd)
         Rij=R(1:3,ia) - R(1:3,ja)
         call calcRij(box, rij, drsq)
         dr = sqrt(drsq)
         e = e + k*(dr-r0)**2
         tmp = 2. * k * (dr-r0)/dr
         di = tmp*rij
         d(1:3,ia) = d(1:3,ia) + di
         d(1:3,ja) = d(1:3,ja) - di
         v = v - reshape( (/rij*di(1), rij*di(2), rij*di(3)/), (/3,3/) )
         ! derivatives of energy and virial wrt k and r0
         dedp(1) = dedp(1) + (dr-r0)**2               !< der of e wrt k
         dedp(2) = dedp(2) - 2.d0*k*(dr-r0)           !< der of e wrt r0
!         dvdp(:,:,1) = dvdp(:,:,1)
      enddo
   endif
   end subroutine bond_harmonic
end module bond_mod

