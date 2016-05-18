!=======================================================================
!//////////////  SUBROUTINE COOL1D_CLOUDY_OLD_TABLES_G  \\\\\\\\\\\\\\\\

      subroutine cool1D_cloudy_old_tables_g(d, de, rhoH, metallicity,
     &                in, jn, kn, is, ie, j, k,
     &                logtem, edot, comp2, ispecies, dom, zr,
     &                icmbTfloor, iClHeat, 
     &                clEleFra, clGridRank, clGridDim,
     &                clPar1, clPar2, clPar3, clPar4, clPar5,
     &                clDataSize, clCooling, clHeating, 
     &                itmask)

!
!  SOLVE CLOUDY METAL COOLING FOR OLD CLOUDY TABLES
!
!  written by: Britton Smith
!  date: September, 2009
!
!  PURPOSE:
!    Solve cloudy cooling by interpolating from the data.
!    This version uses tables formatted for the original
!    Cloudy cooling functionality in Enzo.
!
!  INPUTS:
!    in,jn,kn - dimensions of 3D fields
!
!    d        - total density field
!    de       - electron density field
!
!    rhoH     - total H mass density
!    metallicity - metallicity
!
!    is,ie    - start and end indices of active region (zero based)
!    ispecies - chemistry module (1 - H/He only, 2 - molecular H, 3 - D) 
!    logtem   - natural log of temperature values
!
!    dom      - unit conversion to proper number density in code units
!    zr       - current redshift
!
!    icmbTfloor - flag to include temperature floor from cmb
!    iClHeat    - flag to include cloudy heating
!    clEleFra   - parameter to account for additional electrons from metals 
!    clGridRank - rank of cloudy cooling data grid
!    clGridDim  - array containing dimensions of cloudy data
!    clPar1, clPar2, clPar3, clPar4, clPar5 - arrays containing cloudy grid parameter values
!    clDataSize - total size of flattened 1D cooling data array
!    clCooling  - cloudy cooling data
!    clHeating  - cloudy heating data
!
!    itmask     - iteration mask
!
!  OUTPUTS:
!    update edot with heating/cooling contributions from metals
!
!  PARAMETERS:
!
!-----------------------------------------------------------------------

      implicit NONE
#include "grackle_fortran_types.def"

!  General Arguments

      integer in, jn, kn, is, ie, j, k, ispecies

      real*8 comp2, dom, zr
      R_PREC    d(in,jn,kn), de(in,jn,kn)
      real*8 rhoH(in), metallicity(in), 
     &     logtem(in), edot(in)

!  Cloudy parameters and data

      integer icmbTfloor, iClHeat
      integer*8 clGridRank, clDataSize, clGridDim(3)
      real*8 clEleFra
      real*8 clPar1(clGridDim(1)), clPar2(clGridDim(2)),
     &     clPar3(clGridDim(3)), clPar4(clGridDim(4)),
     &     clPar5(clGridDim(5))
      real*8 clCooling(clDataSize), clHeating(clDataSize)

!  Iteration mask

      logical itmask(in)

!  Parameters

!  Locals

      integer i, q
      real*8 dclPar(clGridRank), inv_log10, log10_tCMB

!  Slice locals

      real*8 log_Z(in), e_frac(in), log_e_frac(in), 
     &     cl_e_frac(in), fh(in), log_n_h(in),
     &     log_cool(in), log_cool_cmb(in), log_heat(in),
     &     edot_met(in), log10tem(in)

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
!=======================================================================

      inv_log10 = 1._DKIND / log(10._DKIND)
      log10_tCMB = log10(comp2)

!     Calculate parameter value slopes

      dclPar(1) = (clPar1(clGridDim(1)) - clPar1(1)) / 
     &     real(clGridDim(1) - 1, DKIND)
      if (clGridRank > 1) then
         dclPar(2) = (clPar2(clGridDim(2)) - clPar2(1)) / 
     &        real(clGridDim(2) - 1, DKIND)
      endif
      if (clGridRank > 2) then
         dclPar(3) = (clPar3(clGridDim(3)) - clPar3(1)) / 
     &        real(clGridDim(3) - 1, DKIND)
      endif
      if (clGridRank > 3) then
         dclPar(4) = (clPar4(clGridDim(4)) - clPar4(1)) /
     &        real(clGridDim(4) - 1, DKIND)
      endif
      if (clGridRank > 4) then
         dclPar(5) = (clPar5(clGridDim(5)) - clPar5(1)) /
     &        real(clGridDim(5) - 1, DKIND)
      endif

      do i=is+1, ie+1
         if ( itmask(i) ) then

            log10tem(i) = logtem(i) * inv_log10

!           Calcualte H mass fraction

            fh(i) = rhoH(i) / d(i,j,k)

!           Calculate proper log(n_H)

            if (clGridRank > 1) then

               log_n_h(i) = log10(rhoH(i) * dom)

            endif

!           Calculate metallicity

            if (clGridRank > 2) then

               log_Z(i) = log10(metallicity(i))

            endif

!           Calculate electron fraction
            
            if (clGridRank > 3) then

               e_frac(i) = 2._DKIND * de(i,j,k) /
     &              (d(i,j,k) * (1._DKIND + fh(i)))
!           Make sure electron fraction is never above 1 
!           which can give bad cooling/heating values when 
!           extrapolating in the Cloudy data.
               log_e_frac(i) = min(log10(e_frac(i)), 0._DKIND)

!           Get extra electrons contributed by metals

               cl_e_frac(i) = e_frac(i) *
     &              (1._DKIND + (2._DKIND * clEleFra * metallicity(i) *
     &              fh(i)) / (1._DKIND + fh(i)))

            endif

!           Call interpolation functions to get heating/cooling

!           Interpolate over temperature.
            if (clGridRank == 1) then
               call interpolate_1D_g(
     &              log10tem(i), clGridDim, clPar1,
     &              dclPar(1), clDataSize, clCooling, log_cool(i))
               edot_met(i) = -10._DKIND**log_cool(i)

!     Ignore CMB term if T >> T_CMB
               if ((icmbTfloor == 1) .and. 
     &              ((log10tem(i) - log10_tCMB) < 2._DKIND)) then
                  call interpolate_1D_g(
     &                 log10_tCMB, clGridDim, clPar1, 
     &                 dclPar(1), clDataSize, clCooling, 
     &                 log_cool_cmb(i))
                  edot_met(i) = edot_met(i) + 10._DKIND**log_cool_cmb(i)
               endif

               if (iClHeat == 1) then
                  call interpolate_1D_g(
     &                 log10tem(i), clGridDim, clPar1, 
     &                 dclPar(1), clDataSize, clHeating, 
     &                 log_heat(i))
                  edot_met(i) = edot_met(i) + 10._DKIND**log_heat(i)
               endif

!           Interpolate over density and temperature.
            else if (clGridRank == 2) then
               call interpolate_2D_g(
     &              log_n_h(i), log10tem(i), clGridDim,
     &              clPar1, dclPar(1), clPar2, dclPar(2),
     &              clDataSize, clCooling, log_cool(i))
               edot_met(i) = -10._DKIND**log_cool(i)

!     Ignore CMB term if T >> T_CMB
               if ((icmbTfloor == 1) .and. 
     &              ((log10tem(i) - log10_tCMB) < 2.0)) then
                  call interpolate_2D_g(
     &                 log_n_h(i), log10_tCMB, clGridDim,
     &                 clPar1, dclPar(1), clPar2, dclPar(2),
     &                 clDataSize, clCooling, log_cool_cmb(i))
                  edot_met(i) = edot_met(i) + 10._DKIND**log_cool_cmb(i)
               endif

               if (iClHeat == 1) then
               call interpolate_2D_g(
     &                 log_n_h(i), log10tem(i), clGridDim,
     &              clPar1, dclPar(1), clPar2, dclPar(2),
     &              clDataSize, clHeating, log_heat(i))
                  edot_met(i) = edot_met(i) + 10._DKIND**log_heat(i)
               endif

!           Interpolate over density, metallicity, and temperature.
            else if (clGridRank == 3) then
               call interpolate_3D_g(
     &              log_n_h(i), log_Z(i), log10tem(i),
     &              clGridDim,
     &              clPar1, dclPar(1), clPar2, dclPar(2),
     &              clPar3, dclPar(3),
     &              clDataSize, clCooling, log_cool(i))
               edot_met(i) = -10._DKIND**log_cool(i)

!     Ignore CMB term if T >> T_CMB
               if ((icmbTfloor == 1) .and. 
     &              ((log10tem(i) - log10_tCMB) < 2._DKIND)) then
                  call interpolate_3D_g(
     &                 log_n_h(i), log_Z(i), log10_tCMB,
     &                 clGridDim,
     &                 clPar1, dclPar(1), clPar2, dclPar(2),
     &                 clPar3, dclPar(3),
     &                 clDataSize, clCooling, log_cool_cmb(i))
                  edot_met(i) = edot_met(i) + 10._DKIND**log_cool_cmb(i)
               endif

               if (iClHeat == 1) then
                  call interpolate_3D_g(
     &                 log_n_h(i), log_Z(i), log10tem(i),
     &                 clGridDim,
     &                 clPar1, dclPar(1), clPar2, dclPar(2),
     &                 clPar3, dclPar(3),
     &                 clDataSize, clHeating, log_heat(i))
                  edot_met(i) = edot_met(i) + 10._DKIND**log_heat(i)
               endif

!           Interpolate over density, metallicity, electron fraction, and temperature.
            else if (clGridRank == 4) then
               call interpolate_4D_g(
     &              log_n_h(i), log_Z(i),
     &              log_e_frac(i), log10tem(i), clGridDim,
     &              clPar1, dclPar(1), clPar2, dclPar(2),
     &              clPar3, dclPar(3), clPar4, dclPar(4),
     &              clDataSize, clCooling, log_cool(i))
               edot_met(i) = -10._DKIND**log_cool(i)

!     Ignore CMB term if T >> T_CMB
               if ((icmbTfloor == 1) .and.
     &              ((log10tem(i) - log10_tCMB) < 2._DKIND)) then
                  call interpolate_4D_g(
     &                 log_n_h(i), log_Z(i),
     &                 log_e_frac(i), log10_tCMB, clGridDim,
     &                 clPar1, dclPar(1), clPar2, dclPar(2),
     &                 clPar3, dclPar(3), clPar4, dclPar(4),
     &                 clDataSize, clCooling, log_cool_cmb(i))
                  edot_met(i) = edot_met(i) + 10._DKIND**log_cool_cmb(i)
               endif

               if (iClHeat == 1) then
                  call interpolate_4D_g(
     &                 log_n_h(i), log_Z(i),
     &                 log_e_frac(i), log10tem(i), clGridDim,
     &                 clPar1, dclPar(1), clPar2, dclPar(2),
     &                 clPar3, dclPar(3), clPar4, dclPar(4),
     &                 clDataSize, clHeating, log_heat(i))
                  edot_met(i) = edot_met(i) + 10._DKIND**log_heat(i)
               endif

!           Interpolate over density, metallicity, electron fraction, redshift,
!           and temperature.
            else
               call interpolate_5D_g(
     &              log_n_h(i), log_Z(i),
     &              log_e_frac(i), zr, log10tem(i), clGridDim,
     &              clPar1, dclPar(1), clPar2, dclPar(2),
     &              clPar3, dclPar(3), clPar4, dclPar(4),
     &              clPar5, dclPar(5),
     &              clDataSize, clCooling, log_cool(i))
               edot_met(i) = -10._DKIND**log_cool(i)

!     Ignore CMB term if T >> T_CMB
               if ((icmbTfloor == 1) .and.
     &              ((log10tem(i) - log10_tCMB) < 2._DKIND)) then
                  call interpolate_5D_g(
     &                 log_n_h(i), log_Z(i),
     &                 log_e_frac(i), zr, log10_tCMB, clGridDim,
     &                 clPar1, dclPar(1), clPar2, dclPar(2),
     &                 clPar3, dclPar(3), clPar4, dclPar(4),
     &                 clPar5, dclPar(5),
     &                 clDataSize, clCooling, log_cool_cmb(i))
                  edot_met(i) = edot_met(i) + 10._DKIND**log_cool_cmb(i)
               endif

               if (iClHeat == 1) then
                  call interpolate_5D_g(
     &                 log_n_h(i), log_Z(i),
     &                 log_e_frac(i), zr, log10tem(i), clGridDim,
     &                 clPar1, dclPar(1), clPar2, dclPar(2),
     &                 clPar3, dclPar(3), clPar4, dclPar(4),
     &                 clPar5, dclPar(5),
     &                 clDataSize, clHeating, log_heat(i))
                  edot_met(i) = edot_met(i) + 10._DKIND**log_heat(i)
               endif

            endif

            if (clGridRank > 3) then
               edot_met(i) = edot_met(i) * cl_e_frac(i)
            endif

            edot(i) = edot(i) + (edot_met(i) * rhoH(i) * d(i,j,k))

         end if
      enddo

      return
      end
