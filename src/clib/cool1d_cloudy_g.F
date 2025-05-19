!=======================================================================
!////////////////////  SUBROUTINE COOL1D_CLOUDY_G  \\\\\\\\\\\\\\\\\\\\\

      subroutine cool1d_cloudy_g(d, rhoH, metallicity,
     &                in, jn, kn, is, ie, j, k,
     &                logtem, edot, comp2, dom, zr,
     &                icmbTfloor, iClHeat, iZscale,
     &                clGridRank, clGridDim,
     &                clPar1, clPar2, clPar3,
     &                clDataSize, clCooling, clHeating, 
     &                itmask)

!
!  SOLVE CLOUDY METAL COOLING
!
!  written by: Britton Smith
!  date: September, 2009
!
!  PURPOSE:
!    Solve cloudy cooling by interpolating from the data.
!
!  INPUTS:
!    in,jn,kn - dimensions of 3D fields
!
!    d        - total density field
!
!    rhoH     - total H mass density
!    metallicity - metallicity
!
!    is,ie    - start and end indices of active region (zero based)
!    logtem   - natural log of temperature values
!
!    dom      - unit conversion to proper number density in code units
!    zr       - current redshift
!
!    icmbTfloor - flag to include temperature floor from cmb
!    iClHeat    - flag to include cloudy heating
!    iZscale    - flag to scale cooling by metallicity
!    clGridRank - rank of cloudy cooling data grid
!    clGridDim  - array containing dimensions of cloudy data
!    clPar1, clPar2, clPar3 - arrays containing cloudy grid parameter values
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

      integer in, jn, kn, is, ie, j, k

      real*8 comp2, dom, zr
      R_PREC d(in,jn,kn)
      real*8 rhoH(in), metallicity(in), logtem(in),
     &       edot(in)

!  Cloudy parameters and data

      integer icmbTfloor, iClHeat, iZscale
      integer*8 clGridRank, clDataSize,
     &     clGridDim(clGridRank)
      real*8 clPar1(clGridDim(1)), clPar2(clGridDim(2)),
     &     clPar3(clGridDim(3)),
     &     clCooling(clDataSize), clHeating(clDataSize)

!  Iteration mask

      MASK_TYPE itmask(in)

!  Parameters

!  Locals

      integer i, get_heat
      integer*8 zindex, zmidpt, zhighpt
      real*8 dclPar(clGridRank), inv_log10, log10_tCMB
      logical end_int

!  Slice locals

      real*8 log_n_h(in),
     &       log_cool(in), log_cool_cmb(in), log_heat(in),
     &       edot_met(in), log10tem(in)

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
!=======================================================================

      end_int = .false.
      get_heat = iClHeat

      inv_log10 = 1._DKIND / log(10._DKIND)
      log10_tCMB = log10(comp2)

!     Calculate parameter value slopes

      dclPar(1) = (clPar1(clGridDim(1)) - clPar1(1)) / 
     &     real(clGridDim(1) - 1, DKIND)
      if (clGridRank .gt. 1) then
         dclPar(2) = (clPar2(clGridDim(2)) - clPar2(1)) / 
     &        real(clGridDim(2) - 1, DKIND)
      endif
      if (clGridRank .gt. 2) then
         dclPar(3) = (clPar3(clGridDim(3)) - clPar3(1)) / 
     &        real(clGridDim(3) - 1, DKIND)
      endif

      do i=is+1, ie+1
         if (itmask(i) .ne. MASK_FALSE) then

            log10tem(i) = logtem(i) * inv_log10

!           Calculate proper log(n_H)

            log_n_h(i) = log10(rhoH(i) * dom)

!           Calculate index for redshift dimension

            if (clGridRank .gt. 2) then

!           Get index for redshift dimension via bisection

               if (zr .le. clPar2(1)) then
                  zindex = 1
               else if (zr .ge. clPar2(clGridDim(2)-1)) then
                  zindex = clGridDim(2)
                  end_int = .true.
                  get_heat = 0
               else if (zr .ge. clPar2(clGridDim(2)-2)) then
                  zindex = clGridDim(2) - 2
               else
                  zindex = 1
                  zhighpt = clGridDim(2) - 2
                  do while ((zhighpt - zindex) .gt. 1)
                     zmidpt = int((zhighpt + zindex) / 2)
                     if (zr .ge. clPar2(zmidpt)) then
                        zindex = zmidpt
                     else
                        zhighpt = zmidpt
                     endif
                  enddo
               endif

            endif

!           Call interpolation functions to get heating/cooling

!           Interpolate over temperature.
            if (clGridRank .eq. 1) then
               call interpolate_1D_g(log10tem(i), clGridDim, clPar1,
     &              dclPar(1), clDataSize, clCooling, log_cool(i))
               edot_met(i) = -10._DKIND**log_cool(i)

!     Ignore CMB term if T >> T_CMB
               if ((icmbTfloor .eq. 1) .and. 
     &              ((log10tem(i) - log10_tCMB) .lt. 2._DKIND)) then
                  call interpolate_1D_g(log10_tCMB, clGridDim, clPar1, 
     &                 dclPar(1), clDataSize, clCooling, 
     &                 log_cool_cmb(i))
                  edot_met(i) = edot_met(i) + 10._DKIND**log_cool_cmb(i)
               endif

               if (get_heat .eq. 1) then
                  call interpolate_1D_g(log10tem(i), clGridDim, clPar1, 
     &                 dclPar(1), clDataSize, clHeating, 
     &                 log_heat(i))
                  edot_met(i) = edot_met(i) + 10._DKIND**log_heat(i)
               endif

!           Interpolate over density and temperature.
            else if (clGridRank .eq. 2) then
               call interpolate_2D_g(log_n_h(i), log10tem(i), clGridDim,
     &              clPar1, dclPar(1), clPar2, dclPar(2),
     &              clDataSize, clCooling, log_cool(i))
               edot_met(i) = -10._DKIND**log_cool(i)

!     Ignore CMB term if T >> T_CMB
               if ((icmbTfloor .eq. 1) .and. 
     &              ((log10tem(i) - log10_tCMB) .lt. 2._DKIND)) then
                  call interpolate_2D_g(log_n_h(i), log10_tCMB, 
     &                 clGridDim, clPar1, dclPar(1), clPar2, dclPar(2),
     &                 clDataSize, clCooling, log_cool_cmb(i))
                  edot_met(i) = edot_met(i) + 10._DKIND**log_cool_cmb(i)
               endif

               if (get_heat .eq. 1) then
               call interpolate_2D_g(log_n_h(i), log10tem(i), clGridDim,
     &              clPar1, dclPar(1), clPar2, dclPar(2),
     &              clDataSize, clHeating, log_heat(i))
                  edot_met(i) = edot_met(i) + 10._DKIND**log_heat(i)
               endif

!           Interpolate over density, redshift, and temperature.
            else if (clGridRank .eq. 3) then
               call interpolate_3Dz_g(log_n_h(i), zr, log10tem(i),
     &              clGridDim,
     &              clPar1, dclPar(1), 
     &              clPar2, zindex,
     &              clPar3, dclPar(3),
     &              clDataSize, clCooling, 
     &              end_int, log_cool(i))
               edot_met(i) = -10._DKIND**log_cool(i)

!     Ignore CMB term if T >> T_CMB
               if ((icmbTfloor .eq. 1) .and. 
     &              ((log10tem(i) - log10_tCMB) .lt. 2._DKIND)) then
                  call interpolate_3Dz_g(log_n_h(i), zr, log10_tCMB,
     &                 clGridDim,
     &                 clPar1, dclPar(1),
     &                 clPar2, zindex,
     &                 clPar3, dclPar(3),
     &                 clDataSize, clCooling, 
     &                 end_int, log_cool_cmb(i))
                  edot_met(i) = edot_met(i) + 10._DKIND**log_cool_cmb(i)
               endif

               if (get_heat .eq. 1) then
                  call interpolate_3Dz_g(log_n_h(i), zr, log10tem(i),
     &                 clGridDim,
     &                 clPar1, dclPar(1),
     &                 clPar2, zindex,
     &                 clPar3, dclPar(3),
     &                 clDataSize, clHeating, 
     &                 end_int, log_heat(i))
                  edot_met(i) = edot_met(i) + 10._DKIND**log_heat(i)
               endif

            else
#ifdef _OPENMP
!$omp critical
#endif
               write(*,*) "Maximum cooling data grid rank is 3!"
#ifdef _OPENMP
!$omp end critical
#endif
               return
            endif

!           Scale cooling by metallicity.

            if (iZscale .eq. 1) then
               edot_met(i) = edot_met(i) * metallicity(i)
            endif

            edot(i) = edot(i) + 
     &           (edot_met(i) * rhoH(i) * rhoH(i))

         end if
      enddo

      return
      end
