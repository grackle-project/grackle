!=======================================================================
!//////////////////  SUBROUTINE CALC_TEMP1D_CLOUDY_G  \\\\\\\\\\\\\\\\\\

      subroutine calc_temp1d_cloudy_shim_g(d, metal, e, rhoH,
     &                in, jn, kn, is, ie, j, k,
     &                tgas, mmw, dom, zr, 
     &                temstart, temend,
     &                gamma, utem, imetal,
     &                clGridRank, clGridDim,
     &                clPar1, clPar2, clPar3,
     &                clDataSize, clMMW,
     &                itmask)

!
!  SHIM TO CALCULATE TEMPERATURE AND MEAN MOLECULAR WEIGHT FROM CLOUDY TABLES
!
!  written by: Matthew Abruzzo
!  date: April, 2023
!
!  PURPOSE:
!    Calculate temperature and mean molecular weight for tabulated cooling.
!    This converts itmask from a logical dtype to integer*4 since logical
!    dtypes can't be trivially passed to a C function.
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

!  Local variable
      integer*4 temp_itmask(in)
      integer i


      do i=is+1, ie+1
         if ( itmask(i) ) then
            temp_itmask(i) = 1
         else
            temp_itmask(i) = 0
         endif
      enddo

      call calc_temp1d_cloudy_g(d, metal, e, rhoH,
     &        in, jn, kn, is, ie, j, k,
     &        tgas, mmw, dom, zr, 
     &        temstart, temend,
     &        gamma, utem, imetal,
     &        clGridRank, clGridDim,
     &        clPar1, clPar2, clPar3,
     &        clDataSize, clMMW,
     &        temp_itmask)

      return
      end
