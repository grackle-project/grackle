!=======================================================================
!//////////////////  SUBROUTINE CALC_TEMP1D_CLOUDY_G  \\\\\\\\\\\\\\\\\\

      subroutine calc_temp1d_cloudy_g(d, metal, e, rhoH,
     &                in, jn, kn, is, ie, j, k,
     &                tgas, mmw, dom, zr, 
     &                temstart, temend,
     &                gamma, utem, imetal,
     &                clGridRank, clGridDim,
     &                clPar1, clPar2, clPar3,
     &                clDataSize, clMMW,
     &                itmask)

!
!  CALCULATE TEMPERATURE AND MEAN MOLECULAR WEIGHT FROM CLOUDY TABLES
!
!  written by: Britton Smith
!  date: May, 2015
!
!  PURPOSE:
!    Calculate temperature and mean molecular weight for tabulated cooling.
!
!  INPUTS:
!    in,jn,kn - dimensions of 3D fields
!
!    d        - density field
!    metal    - metal density field
!    e        - internal energy field
!
!    rhoH     - total H mass density
!
!    is,ie    - start and end indices of active region (zero based)
!    tgas     - temperature values
!    mmw      - mean molecular weight values
!
!    dom      - unit conversion to proper number density in code units
!    zr       - current redshift
!    temstart, temend - start and end of temperature range for rate table
!    gamma    - adiabatic index
!    utem     - temperature units
!    imetal   - flag if metal field is active (0 = no, 1 = yes)
!
!    clGridRank - rank of cloudy cooling data grid
!    clGridDim  - array containing dimensions of cloudy data
!    clPar1, clPar2, clPar3 - arrays containing cloudy grid parameter values
!    clDataSize - total size of flattened 1D cooling data array
!    clMMW      - cloudy mmw data
!
!    itmask     - iteration mask
!
!  OUTPUTS:
!    fills temperature and mmw arrays
!
!  PARAMETERS:
!
!-----------------------------------------------------------------------

      implicit NONE
#include "grackle_fortran_types.def"

!  General Arguments

      integer in, jn, kn, is, ie, j, k,
     &     imetal

      real*8 dom, zr, gamma, temstart, temend, utem
      R_PREC d(in,jn,kn), metal(in,jn,kn), e(in,jn,kn)
      real*8 rhoH(in), tgas(in),
     &       mmw(in), logtem(in)

!  Cloudy parameters and data

      integer*8 clGridRank, clDataSize,
     &     clGridDim(clGridRank)
      real*8 clPar1(clGridDim(1)), clPar2(clGridDim(2)),
     &     clPar3(clGridDim(3)),
     &     clMMW(clDataSize)

!  Iteration mask

      logical itmask(in)

!  Parameters

      integer ti_max
      real*8 mu_metal
      parameter (mu_metal = 16._DKIND)    ! approx. mean molecular weight of metals
      parameter (ti_max = 20)

!  Locals

      integer i, ti
      integer*8 zindex, zmidpt, zhighpt
      real*8 dclPar(clGridRank), inv_log10,
     &     muold, munew
      logical end_int

!  Slice locals

      real*8 log_n_h(in), log10tem(in)

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
!=======================================================================

      end_int = .false.

      inv_log10 = 1._DKIND / log(10._DKIND)

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

!     Calculate index for redshift dimension

      if (clGridRank .gt. 2) then

!     Get index for redshift dimension via bisection

         if (zr .le. clPar2(1)) then
            zindex = 1
         else if (zr .ge. clPar2(clGridDim(2)-1)) then
            zindex = clGridDim(2)
            end_int = .true.
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

      do i=is+1, ie+1
         if ( itmask(i) ) then
!              Calculate proper log(n_H)
               log_n_h(i) = log10(rhoH(i) * dom)
         endif
      enddo

      do i=is+1, ie+1
         if ( itmask(i) ) then
            munew = 1._DKIND
            do ti = 1, ti_max
               muold = munew

               tgas(i) = max((gamma - 1._DKIND) * 
     &              e(i,j,k) * munew * utem,
     &              temstart)
               logtem(i) = log(tgas(i))

               log10tem(i) = logtem(i) * inv_log10

!              Call interpolation functions to get heating/cooling

!              Interpolate over temperature.
               if (clGridRank .eq. 1) then
                  call interpolate_1D_g(log10tem(i), clGridDim,
     &                 clPar1, dclPar(1),
     &                 clDataSize, clMMW, munew)

!              Interpolate over density and temperature.
               else if (clGridRank .eq. 2) then
                  call interpolate_2D_g(log_n_h(i), log10tem(i), 
     &                 clGridDim,
     &                 clPar1, dclPar(1),
     &                 clPar2, dclPar(2),
     &                 clDataSize, clMMW, munew)

!              Interpolate over density, redshift, and temperature.
               else if (clGridRank .eq. 3) then
                  call interpolate_3Dz_g(log_n_h(i), zr, log10tem(i),
     &                 clGridDim,
     &                 clPar1, dclPar(1), 
     &                 clPar2, zindex,
     &                 clPar3, dclPar(3),
     &                 clDataSize, clMMW, 
     &                 end_int, munew)

               else
#ifdef _OPENMP
!$omp critical
#endif
                  write(*,*) "Maximum mmw data grid rank is 3!"
#ifdef _OPENMP
!$omp end critical
#endif
                  return
               endif

               munew = 0.5 * (munew + muold)
               tgas(i) = tgas(i) * munew / muold

               if (abs((munew/muold) - 1._DKIND) .le. 
     &              1.e-2_DKIND) then
                  muold = munew

!     Add metal species to mean molecular weight

                  if (imetal .eq. 1) then
                     munew = d(i,j,k) / (d(i,j,k) / munew +
     &                    metal(i,j,k) / mu_metal)
                     tgas(i) = tgas(i) * munew / muold
                  endif

                  mmw(i) = munew
                  go to 9998
               endif

            enddo

            mmw(i) = munew
#ifdef _OPENMP
!$omp critical
#endif
            write(6,*) 'Mean molecular weight not converged!', 
     &           munew, muold, abs((munew/muold) - 1._DKIND)
#ifdef _OPENMP
!$omp end critical
#endif

 9998       continue

         end if
      enddo

      return
      end
