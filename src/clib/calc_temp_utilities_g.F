!=======================================================================
!//////////////  SUBROUTINE CALC_PARAMETER_VALUE_SLOPES  \\\\\\\\\\\\\\\

      subroutine calc_parameter_value_slopes(dclPar, clGridRank, 
     &            clGridDim, clPar1, clPar2, clPar3)
   
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
      real*8, optional :: clPar2(clGridDim(2)), clPar3(clGridDim(3))
    
!  Slopes of parameter values   
      real*8 dclPar(clGridRank)
   
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
!=======================================================================
   
! TODO: check presence of clPar2/3 wrt clGridRank     
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

      return
      end
 
 