#include "phys_const.def"

!=======================================================================
!////////////////////  SUBROUTINE CALC_KAPPA_GR_G  \\\\\\\\\\\\\\\\\\\\\

      subroutine calc_kappa_gr_g(
     &     tdust, kgr, itmask, in, is, ie, t_subl
     &   , gr_N, gr_Size, gr_dT, gr_Td, logalsp, idspecies)

!  CALCULATE GRAIN PLANK MEAN OPACITY
!
!  written by: Britton Smith
!  date:       September, 2011
!
!  PURPOSE:
!    Calculate grain plank mean opacity
!
!  TODO:
!    this docstring should EXPLICITLY document how the action of this
!    function changes based on whether we are using the single-field
!    dust model or the multi-species dust model.
!
!    In the classic single single-field model, the returned opacities
!    have units of cm^2/g and they are measured "per unit grain mass"
!
!    In the multi-species dust model, I THINK that the returned
!    opacities have units of cm^2/g and they are measured "per unit
!    gas mass."
!
!  INPUTS:
!     in       - i dimension of 3D fields
!
!     tdust    - dust temperature
!
!     is,ie    - start and end indices of active region (zero based)
!
!     itmask   - iteration mask
!
!     t_subl   - grain sublimation temperature
!
!  OUTPUTS:
!     kgr      - opacities
!
!  PARAMETERS:
!
!-----------------------------------------------------------------------

      implicit NONE
#include "grackle_fortran_types.def"

!  Arguments

      integer in, is, ie
      integer idspecies
      real*8 t_subl
      real*8 tdust(in)

!  Iteration mask

      MASK_TYPE itmask(in)

!  Parameters

! grain opacity from Omukai (2000, equation 17) normalized by
! the local dust-to-gas ratio, which in this work is 0.934e-2.
      real*8, parameter :: kgr1 = 4.0e-4_DKIND / 0.00934d0

! We then apply a Td^2 scaling up to 200 K following Dopcke et al. (2011).
! However, note that Dopcke et al. (2011) adopt a different
! normalization (i.e., for kgr1).
! Also note that Omukai (2000) says the Td^2 proportionality is only valid
! for Td < 50 K. We go with this because Omukai (2000) does not suggest
! what should be done for Td > 50 K.
      real*8, parameter :: kgr200 = 16.0_DKIND / 0.00934d0

!  Opacity table

      integer gr_N(1), gr_Size
      real*8  gr_dT, gr_Td(gr_N(1))
      R_PREC  logalsp(gr_N(1),in)

!  Locals

      integer i
      R_PREC  logalsp1(gr_Size), logkgr
      integer*8 gr_N_i64(1)
      real*8 log10tdust

!  Slice Locals

      real*8 kgr(in)

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
!=======================================================================
      if(idspecies.ne.0) then
         gr_N_i64(1) = int(gr_N(1), 8)
      endif


      do i = is+1, ie+1
      if ( itmask(i) .ne. MASK_FALSE ) then

         if(idspecies.eq.0) then

!           Temperature dependence from Dopcke et al. (2011).
!           Normalized to Omukai (2000).
!           See comment above for note about Td dependence for kgr.
            if (tdust(i) .lt. 200._DKIND) then
               kgr(i) = kgr1 * tdust(i)**2
            else if (tdust(i) .lt. t_subl) then
               kgr(i) = kgr200
            else
               kgr(i) = max(tiny, 
     &              (kgr200 * (tdust(i) / 1.5e3_DKIND)**(-12)))
            endif

         else

            logalsp1(:) = logalsp(:,i)
            log10tdust = log10(tdust(i))
       
            call interpolate_1D_g(
     &           log10tdust, gr_N_i64, gr_Td, gr_dT,
     &           gr_N_i64(1), logalsp1, logkgr)
            kgr(i) = 10._DKIND**logkgr

         endif

      endif
      enddo

      return
      end

! TODO: to be removed when transcription is done
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
