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

      logical itmask(in)

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

#ifdef USE_ORIGINAL_CODE

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
#else
      call calc_parameter_value_slopes(dclPar, clGridRank,
     &            clGridDim, clPar1, clPar2, clPar3)
#endif

      do i=is+1, ie+1
         if ( itmask(i) ) then

            log10tem(i) = logtem(i) * inv_log10

!           Calculate proper log(n_H)

            log_n_h(i) = log10(rhoH(i) * dom)

!           Calculate index for redshift dimension

            if (clGridRank .gt. 2) then

!           Get index for redshift dimension via bisection
#ifdef USE_ORIGINAL_CODE
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
#else
               call compute_redshift_dimension(zr, clGridDim(2), clPar2,
     &            zindex, end_int, get_heat)
#endif
            endif

!           Call interpolation functions to get heating/cooling

!           Interpolate over temperature.
            if (clGridRank .eq. 1) then
#ifdef USE_ORIGINAL_CODE
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
#else
               call cool_rank1_interpolation(get_heat, icmbTfloor,
     &              log10tem(i), log10_tCMB, clGridRank, clDataSize,
     &              clGridDim, clPar1, clCooling, clHeating,
     &              dclPar(1), edot_met(i))


#endif
!           Interpolate over density and temperature.
            else if (clGridRank .eq. 2) then
#ifdef USE_ORIGINAL_CODE
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
#else
               call cool_rank2_interpolation(get_heat, icmbTfloor,
     &            log_n_h(i), log10tem(i), log10_tCMB,
     &            clGridRank, clDataSize, clGridDim,
     &            clPar1, clPar2, clCooling, clHeating,
     &            dclPar(1), dclPar(2), edot_met(i))

#endif
!           Interpolate over density, redshift, and temperature.
            else if (clGridRank .eq. 3) then
#ifdef USE_ORIGINAL_CODE
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
#else
               call cool_rank3_interpolation(get_heat, icmbTfloor,
     &            log_n_h(i), zr, log10tem(i), log10_tCMB,
     &            clGridRank, clDataSize, clGridDim, clPar1,
     &            clPar2, clPar3, clCooling, clHeating, dclPar(1),
     &            zindex, dclPar(3), edot_met(i))
#endif

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

!=======================================================================
!///////////////  SUBROUTINE COMPUTE_REDSHIFT_DIMENSION  \\\\\\\\\\\\\\\

      subroutine compute_redshift_dimension(zr, clGridDim_2, clPar2,
     &   zindex, end_int, get_heat)

!
!  COMPUTE INDEX FOR REDSHIFT DIMENSION VIA BISECTION
!
!  written by:
!  date: May, 2024
!
!  PURPOSE:
!
!  INPUTS:
!
!  OUTPUTS:
!     zindex
!     end_int
!     get_heat
!
!  PARAMETERS:
!
!-----------------------------------------------------------------------

      implicit NONE
#include "grackle_fortran_types.def"

      real*8 zr
      integer*8 zindex, zmidpt, zhighpt, clGridDim_2
      real*8 clPar2(clGridDim_2)
      logical end_int
      integer get_heat

      if (zr .le. clPar2(1)) then
         zindex = 1
      else if (zr .ge. clPar2(clGridDim_2-1)) then
         zindex = clGridDim_2
         end_int = .true.
         get_heat = 0
      else if (zr .ge. clPar2(clGridDim_2-2)) then
         zindex = clGridDim_2 - 2
      else
         !TODO: create minimal bisection subroutine?
         zindex = 1
         zhighpt = clGridDim_2 - 2
         do while ((zhighpt - zindex) .gt. 1)
            zmidpt = int((zhighpt + zindex) / 2)
            if (zr .ge. clPar2(zmidpt)) then
               zindex = zmidpt
            else
               zhighpt = zmidpt
            endif
         enddo
      endif

      return
      end

!=======================================================================
!///////////////  SUBROUTINE COOL_RANK1_INTERPOLATION  \\\\\\\\\\\\\\\\\

      subroutine cool_rank1_interpolation(get_heat, icmbTfloor,
     &      log10tem_i, log10_tCMB, clGridRank, clDataSize,
     &      clGridDim, clPar1, clCooling, clHeating, dclPar_1,
     &      edot_met_i)

!
!  INTERPOLATE RANK 1 CASE (Interpolate over temperature.)
!
!  written by:
!  date: May, 2024
!
!  PURPOSE:
!
!  INPUTS:
!
!  OUTPUTS:
!    edot_met_i
!
!  PARAMETERS:
!
!-----------------------------------------------------------------------

      implicit NONE
#include "grackle_fortran_types.def"

      integer get_heat, icmbTfloor
      real*8 log10tem_i, edot_met_i
      integer*8 clGridRank, clDataSize, clGridDim(clGridRank)
      real*8 clPar1(clGridDim(1)), clCooling(clDataSize),
     &    clHeating(clDataSize), dclPar_1, log10_tCMB
      ! Local variables
      real*8 log_cool_i, log_cool_cmb_i, log_heat_i

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
!=======================================================================

      call interpolate_1D_g(log10tem_i, clGridDim, clPar1,
     &                      dclPar_1, clDataSize, clCooling, log_cool_i)
      edot_met_i = -10._DKIND**log_cool_i

!     Ignore CMB term if T >> T_CMB
      if ((icmbTfloor .eq. 1) .and.
     &     ((log10tem_i - log10_tCMB) .lt. 2._DKIND)) then
         call interpolate_1D_g(log10_tCMB, clGridDim, clPar1,
     &        dclPar_1, clDataSize, clCooling,
     &        log_cool_cmb_i)
         edot_met_i = edot_met_i + 10._DKIND**log_cool_cmb_i
      endif

      if (get_heat .eq. 1) then
         call interpolate_1D_g(log10tem_i, clGridDim, clPar1,
     &        clPar1(1), clDataSize, clHeating,
     &        log_heat_i)
         edot_met_i = edot_met_i + 10._DKIND**log_heat_i
      endif

      return
      end

!=======================================================================
!///////////////  SUBROUTINE COOL_RANK2_INTERPOLATION  \\\\\\\\\\\\\\\\\

      subroutine cool_rank2_interpolation(get_heat, icmbTfloor,
     &      log_n_h_i, log10tem_i, log10_tCMB,
     &      clGridRank, clDataSize, clGridDim,
     &      clPar1, clPar2, clCooling, clHeating,
     &      dclPar_1, dclPar_2, edot_met_i)

!
!  INTERPOLATE RANK 2 CASE (Interpolate over density and temperature.)
!
!  written by:
!  date: May, 2024
!
!  PURPOSE:
!
!  INPUTS:
!
!  OUTPUTS:
!    edot_met_i
!
!  PARAMETERS:
!
!-----------------------------------------------------------------------

      implicit NONE
#include "grackle_fortran_types.def"

      integer get_heat, icmbTfloor
      real*8 log10tem_i, edot_met_i, log_n_h_i
      integer*8 clGridRank, clDataSize, clGridDim(clGridRank)
      real*8 clPar1(clGridDim(1)), clPar2(clGridDim(2)),
     &   clCooling(clDataSize), clHeating(clDataSize), dclPar_1,
     &   dclPar_2, log10_tCMB
      ! Local variables
      real*8 log_cool_i, log_cool_cmb_i, log_heat_i

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
!=======================================================================

      call interpolate_2D_g(log_n_h_i, log10tem_i, clGridDim,
     &        clPar1, dclPar_1, clPar2, dclPar_2,
     &        clDataSize, clCooling, log_cool_i)
      edot_met_i = -10._DKIND**log_cool_i

!     Ignore CMB term if T >> T_CMB
      if ((icmbTfloor .eq. 1) .and.
     &     ((log10tem_i - log10_tCMB) .lt. 2._DKIND)) then
         call interpolate_2D_g(log_n_h_i, log10_tCMB,
     &      clGridDim, clPar1, dclPar_1, clPar2, dclPar_2,
     &      clDataSize, clCooling, log_cool_cmb_i)
         edot_met_i = edot_met_i + 10._DKIND**log_cool_cmb_i
      endif

      if (get_heat .eq. 1) then
         call interpolate_2D_g(log_n_h_i, log10tem_i, clGridDim,
     &      clPar1, dclPar_1, clPar2, dclPar_2,
     &      clDataSize, clHeating, log_heat_i)
         edot_met_i = edot_met_i + 10._DKIND**log_heat_i
      endif

      return
      end

!=======================================================================
!///////////////  SUBROUTINE COOL_RANK3_INTERPOLATION  \\\\\\\\\\\\\\\\\

      subroutine cool_rank3_interpolation(get_heat, icmbTfloor,
     &      log_n_h_i, zr, log10tem_i, log10_tCMB,
     &      clGridRank, clDataSize, clGridDim,
     &      clPar1, clPar2, clPar3,
     &      clCooling, clHeating,
     &      dclPar_1, zindex, dclPar_3, edot_met_i)

!
!  INTERPOLATE RANK 3 CASE (Interpolate over density, redshift, and temperature.)
!
!  written by:
!  date: May, 2024
!
!  PURPOSE:
!
!  INPUTS:
!
!  OUTPUTS:
!    edot_met_i
!
!  PARAMETERS:
!
!-----------------------------------------------------------------------

      implicit NONE
#include "grackle_fortran_types.def"

      integer get_heat, icmbTfloor
      real*8 log10tem_i, edot_met_i, log_n_h_i, zr
      integer*8 zindex, clGridRank, clDataSize, clGridDim(clGridRank)
      real*8 clPar1(clGridDim(1)), clPar2(clGridDim(2)),
     &   clPar3(clGridDim(3)), clCooling(clDataSize),
     &   clHeating(clDataSize), dclPar_1, dclPar_3, log10_tCMB
      logical end_int
      ! Local variables
      real*8 log_cool_i, log_cool_cmb_i, log_heat_i

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
!=======================================================================

      call interpolate_3Dz_g(log_n_h_i, zr, log10tem_i,
     &   clGridDim,
     &   clPar1, dclPar_1,
     &   clPar2, zindex,
     &   clPar3, dclPar_3,
     &   clDataSize, clCooling,
     &   end_int, log_cool_i)
      edot_met_i = -10._DKIND**log_cool_i

!     Ignore CMB term if T >> T_CMB
      if ((icmbTfloor .eq. 1) .and.
     &   ((log10tem_i - log10_tCMB) .lt. 2._DKIND)) then
         call interpolate_3Dz_g(log_n_h_i, zr, log10_tCMB,
     &      clGridDim,
     &      clPar1, dclPar_1,
     &      clPar2, zindex,
     &      clPar3, dclPar_3,
     &      clDataSize, clCooling,
     &      end_int, log_cool_cmb_i)
         edot_met_i = edot_met_i + 10._DKIND**log_cool_cmb_i
      endif

      if (get_heat .eq. 1) then
         call interpolate_3Dz_g(log_n_h_i, zr, log10tem_i,
     &      clGridDim,
     &      clPar1, dclPar_1,
     &      clPar2, zindex,
     &      clPar3, dclPar_3,
     &      clDataSize, clHeating,
     &      end_int, log_heat_i)
         edot_met_i = edot_met_i + 10._DKIND**log_heat_i
      endif

      return
      end
