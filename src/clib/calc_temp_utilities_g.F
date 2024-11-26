!=======================================================================
!//////////////  SUBROUTINE CALC_PARAMETER_VALUE_SLOPES  \\\\\\\\\\\\\\\

      subroutine calc_parameter_value_slopes(dclPar, clGridRank, 
     &            clGridDim, clPar1, clPar2, clPar3, clpar4, clpar5)
   
!
!  CALCULATE PARAMETER VALUE SLOPES
!
!  written by: 
!  date: May
!
!  PURPOSE:
!    Calculate the slopes of the parameter values
!
!  INPUTS:
!
!    clGridRank - rank of cloudy cooling data grid
!    clGridDim  - array containing dimensions of cloudy data
!    clPar1, clPar2 (optional), clPar3 (optional) - arrays containing 
!           cloudy grid parameter values
!
!  OUTPUTS:
!    dclPar - parameter slopes
!
!  PARAMETERS:
!
!-----------------------------------------------------------------------
   
      implicit NONE
#include "grackle_fortran_types.def"
   
!  Cloudy parameters and data
   
      integer*8 clGridRank, clGridDim(clGridRank)
      real*8 clPar1(clGridDim(1))
      real*8, optional :: clPar2(clGridDim(2)), clPar3(clGridDim(3)),
     &     clPar4(clGridDim(4)), clPar5(clGridDim(5))
    
!  Slopes of parameter values   
      real*8 dclPar(clGridRank)
   
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
!=======================================================================

! TODO: check presence of clPar2/3/4/5 wrt clGridRank
      dclPar(1) = (clPar1(clGridDim(1)) - clPar1(1)) /
     &     real(clGridDim(1) - 1, DKIND)
      if (clGridRank .gt. 1) then 
!        if(present(clPar2)) then
          dclPar(2) = (clPar2(clGridDim(2)) - clPar2(1)) /
     &        real(clGridDim(2) - 1, DKIND)
!        else
!            write(*,*) 'Error in calc_parameter_value_slopes.'
!        endif
      endif
      if (clGridRank .gt. 2) then
!        if(present(clPar3)) then
          dclPar(3) = (clPar3(clGridDim(3)) - clPar3(1)) /
     &        real(clGridDim(3) - 1, DKIND)
!        else
!            write(*,*) 'Error in calc_parameter_value_slopes.'
!        endif
      endif
      if (clGridRank .gt. 3) then
!        if(present(clPar4)) then
            dclPar(4) = (clPar4(clGridDim(4)) - clPar4(1)) /
     &        real(clGridDim(4) - 1, DKIND)
!        else
!            write(*,*) 'Error in calc_parameter_value_slopes.'
!        endif
      endif
      if (clGridRank .gt. 4) then
!        if(present(clPar5)) then
            dclPar(5) = (clPar5(clGridDim(5)) - clPar5(1)) /
     &        real(clGridDim(5) - 1, DKIND)
!        else
!            write(*,*) 'Error in calc_parameter_value_slopes.'
!        endif
      endif

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
      integer, optional :: get_heat

      if (zr .le. clPar2(1)) then
          zindex = 1
      else if (zr .ge. clPar2(clGridDim_2-1)) then
          zindex = clGridDim_2
          end_int = .true.
          if (present(get_heat)) then
              get_heat = 0
          endif
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