!--------------------------------------------------------------------------
!>  @author  GSF
!!  @brief constants, physical constants and units conversion factors
!!  @par Description
!!  @details  *
!<
!--------------------------------------------------------------------------
module consts
   !.... mathematical constants
   double precision, parameter ::     PI= 3.14159265358979323844d0
   double precision, parameter ::  TWOPI= 2.d0*PI
   double precision, parameter :: FOURPI= 4.d0*PI
   double precision, parameter ::    TSP= 2.d0 / sqrt(PI)
   double precision, parameter ::   SQPI= sqrt(PI) !  sqrt(pi)
   double precision, parameter :: DEG2RAD=  PI / 180.d0  !  degrees --> Radian 
   !... physical constants
   double precision, parameter :: BOLTZMANN = 1.9872041d-3 ! kcal/mol/Kelvin
   double precision, parameter :: COULOMB = 332.0522173d0
   double precision, parameter :: HBAR =0.0151704885d0 ! kcal/mol * ps
   double precision, parameter :: SPDLGHT = 2.997924580d6  ! Ang * ps-1
   double precision, parameter :: NAVOG = 6.022140857d23

   double precision, parameter :: BARATM = 0.986923d0   ! bar --> to Atm
   double precision, parameter :: KCALKLVN= 1.d0 / BOLTZMANN ! 503.21 
   double precision, parameter :: PRESCON=6.85695d+4/BARATM  ! kcal/mole/Ang**3 to bar !Atm
   double precision, parameter :: MASSCON=418.4d0!at.weig->(kcal/mol)*(ps/Ang)^2
   double precision, parameter :: DENSCON = 1.66053872801494672163d0*MASSCON
!   double precision, parameter :: CHARGECON = dsqrt(COULOMB)
   double precision, parameter :: CHARGECON = 18.22261534447435150793d0 !tinker
   double precision, parameter :: STP_NUMDENS =  1.d0/PRESCON/(273.15/KCALKLVN)/BARATM ! [1/A3] standard num dens

   double precision, parameter :: DEBYE  = 4.8033324d0/CHARGECON
   double precision, parameter :: FS2PS=0.001d0
   double precision, parameter :: CMIPSI=.18836515673088532773d0! cm-1 --> ps-1

   double precision, parameter :: CALJOULE=4.184d0! kcal ro kjoule
   double precision, parameter :: EVKCAL = 23.060541945329334d0

!   double precision, parameter :: KCALKJOULE = 4.184d0! from kcal to kjoule
!   double precision, parameter :: KCALKLVN=503.21659237841666799453d0
!   double precision, parameter :: CHARGECON = 18.22261720426243437986d0
!   double precision, parameter :: CHARGECON = 18.22261534447435150793d0 !tinker


   integer, parameter :: anamesize=8
end module consts

