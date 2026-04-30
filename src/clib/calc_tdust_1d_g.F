#include "phys_const.def"

!=======================================================================
!////////////////////  SUBROUTINE CALC_GR_BALANCE  \\\\\\\\\\\\\\\\\\\\\

      subroutine calc_gr_balance_g(
     &     tdust, tgas, kgr, trad4, gasgr, gamma_isrf, nh,
     &     itmask, sol, in, is, ie)

!  CALCULATE GRAIN HEAT BALANCE
!
!  written by: Britton Smith
!  date:       September, 2019
!  modified1:
!
!  PURPOSE:
!    Calculate grain heating/cooling balance

!  TODO:
!    this docstring should EXPLICITLY document how the action of this
!    function changes based on whether we are using the single-field
!    dust model or the multi-species dust model.
!
!    In the classic single single-field model, the returned values
!    are emission/absorption rate per unit grain mass with units of
!    erg/s/g
!
!    In the multi-species dust model, I **THINK** that the returned
!    values are emission/absorption rate per unit gas mass with units
!    of erg/s/g (But I'm not entirely unsure).
!
!  INPUTS:
!     in       - i dimension of 3D fields
!
!     tdust    - dust temperature
!     tgas     - gas temperature
!     kgr      - grain opacity
!     trad4    - CMB temperature to 4th power
!     gasgr    - gas/grain heat transfer rate
!     gamma_isrf - heating from interstellar radiation field
!     nh       - hydrogen number density
!
!     is,ie    - start and end indices of active region (zero based)
!
!     itmask   - iteration mask
!
!
!  OUTPUTS:
!     sol      - heating/cooling balance (heating - cooling)
!
!  PARAMETERS:
!
!-----------------------------------------------------------------------

      implicit NONE
#include "grackle_fortran_types.def"

!  Arguments

      integer in, is, ie
      real*8 tdust(in), tgas(in), kgr(in), gasgr(in), nh(in), trad4,
     &     gamma_isrf(in)

!  Iteration mask

      MASK_TYPE itmask(in)

!  Parameters

      real*8, parameter :: radf = 4._DKIND * sigma_sb

!  Locals

      integer i

!  Slice Locals

      real*8 sol(in)

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
!=======================================================================

      do i = is+1, ie+1
         if ( itmask(i) .ne. MASK_FALSE ) then

            sol(i) = gamma_isrf(i) + radf * kgr(i) *
     &           (trad4 - tdust(i)**4) +
     &           (gasgr(i) * nh(i) *
     &           (tgas(i) - tdust(i)))
       ! Historically, the following comment was present here:
       !     emission/absorption rate per unit grain mass [erg/s/g]
       !     for Z = Zsun (default)
       ! This comment is **ONLY** correct when the function used as
       ! part of the single-field dust model. See the docstring for
       ! more details.

         endif
      enddo

      return
      end
