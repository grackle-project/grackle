#include "fortran.def"

!=======================================================================
!/////////////////////  SUBROUTINE MMW1D_CLOUDY_G  \\\\\\\\\\\\\\\\\\\\\

      subroutine mmw1d_cloudy_g(d, metal, e, rhoH,
     &                in, jn, kn, is, ie, j, k,
     &                logtem, mmw, dom, zr, imetal,
     &                clGridRank, clGridDim,
     &                clPar1, clPar2, clPar3,
     &                clDataSize, clMMW,
     &                itmask)

!
!  CALCULATE MEAN MOLECULAR WEIGHT AND TEMPERATURE FROM CLOUDY TABLES
!
!  written by: Britton Smith
!  date: September, 2009
!
!  PURPOSE:
!    Calculate mean molecular weight for tabulated cooling.
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
!    fills mmw and temperature array
!
!  PARAMETERS:
!
!-----------------------------------------------------------------------

      implicit NONE
#include "grackle_fortran_types.def"

!  General Arguments

      integer in, jn, kn, is, ie, j, k,
     &     imetal

      real*8 dom, zr
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

!  Locals

      integer i
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

      do i=is+1, ie+1
         if ( itmask(i) ) then

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
     &              dclPar(1), clDataSize, clMMW, mmw(i))

!           Interpolate over density and temperature.
            else if (clGridRank .eq. 2) then
               call interpolate_2D_g(log_n_h(i), log10tem(i), clGridDim,
     &              clPar1, dclPar(1), clPar2, dclPar(2),
     &              clDataSize, clMMW, mmw(i))

!           Interpolate over density, redshift, and temperature.
            else if (clGridRank .eq. 3) then
               call interpolate_3Dz_g(log_n_h(i), zr, log10tem(i),
     &              clGridDim,
     &              clPar1, dclPar(1), 
     &              clPar2, zindex,
     &              clPar3, dclPar(3),
     &              clDataSize, clMMW, 
     &              end_int, mmw(i))

            else
               write(*,*) "Maximum cooling data grid rank is 3!"
               return
            endif

         end if
      enddo

      return
      end
