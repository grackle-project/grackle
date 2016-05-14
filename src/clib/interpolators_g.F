!=======================================================================
!////////////////////  SUBROUTINE INTERPOLATE_1D  \\\\\\\\\\\\\\\\\\\\\\

      subroutine interpolate_1D_g(
     &     input1, gridDim, gridPar1, dgridPar1, 
     &     dataSize, dataField, value)

      implicit NONE
#include "grackle_fortran_types.def"

!  General Arguments

      integer*8 dataSize
      integer*8 gridDim(1)
      real*8 input1, value
      real*8 gridPar1(gridDim(1)), dgridPar1
      real*8 dataField(dataSize)

!  Locals

      integer*8 index1
      real*8 slope

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
!=======================================================================

!     Calculate interpolation index

      index1 = min(gridDim(1)-1,
     &     max(1,int((input1-gridPar1(1))/dgridPar1, DIKIND)+1))

!     Interpolate over parameter 1

      slope = (dataField(index1+1) - dataField(index1)) /
     &     (gridPar1(index1+1) - gridPar1(index1))

      value = (input1 - gridPar1(index1)) * slope + dataField(index1)

      return
      end

!=======================================================================
!////////////////////  SUBROUTINE INTERPOLATE_2D  \\\\\\\\\\\\\\\\\\\\\\

      subroutine interpolate_2D_g(
     &     input1, input2, gridDim, 
     &     gridPar1, dgridPar1,
     &     gridPar2, dgridPar2,
     &     dataSize, dataField, value)

      implicit NONE
#include "grackle_fortran_types.def"

!  General Arguments

      integer*8 dataSize
      integer*8 gridDim(2)
      real*8 input1, input2, value
      real*8 gridPar1(gridDim(1)), dgridPar1,
     &       gridPar2(gridDim(2)), dgridPar2
      real*8 dataField(dataSize)

!  Locals

      integer*8 index1, index2, int_index
      integer q
      real*8 slope, value2(2)

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
!=======================================================================

!     Calculate interpolation indices

      index1 = min(gridDim(1)-1,
     &     max(1,int((input1-gridPar1(1))/dgridPar1, DIKIND)+1))
      index2 = min(gridDim(2)-1,
     &     max(1,int((input2-gridPar2(1))/dgridPar2, DIKIND)+1))

      do q=1, 2

!     interpolate over parameter 2

         int_index = (q+index1-2) * gridDim(2) + index2

         slope = (dataField(int_index+1) - dataField(int_index)) /
     &        (gridPar2(index2+1) - gridPar2(index2))

         value2(q) = (input2 - gridPar2(index2)) * slope + 
     &        dataField(int_index)

      enddo

!     interpolate over parameter 1

      slope = (value2(2) - value2(1)) / 
     &     (gridPar1(index1+1) - gridPar1(index1))

      value = (input1 - gridPar1(index1)) * slope + value2(1)

      return
      end

!=======================================================================
!////////////////////  SUBROUTINE INTERPOLATE_3D  \\\\\\\\\\\\\\\\\\\\\\

      subroutine interpolate_3D_g(
     &     input1, input2, input3, gridDim,
     &     gridPar1, dgridPar1,
     &     gridPar2, dgridPar2,
     &     gridPar3, dgridPar3,
     &     dataSize, dataField, value)

      implicit NONE
#include "grackle_fortran_types.def"

!  General Arguments

      integer*8 dataSize
      integer*8 gridDim(3)
      real*8 input1, input2, input3, value
      real*8 gridPar1(gridDim(1)), dgridPar1,
     &       gridPar2(gridDim(2)), dgridPar2,
     &       gridPar3(gridDim(3)), dgridPar3
      real*8 dataField(dataSize)

!  Locals

      integer*8 index1, index2, index3, int_index
      integer q, w
      real*8 slope, value3(2), value2(2)

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
!=======================================================================

!     Calculate interpolation indices

      index1 = min(gridDim(1)-1,
     &     max(1,int((input1-gridPar1(1))/dgridPar1, DIKIND)+1))
      index2 = min(gridDim(2)-1,
     &     max(1,int((input2-gridPar2(1))/dgridPar2, DIKIND)+1))
      index3 = min(gridDim(3)-1,
     &     max(1,int((input3-gridPar3(1))/dgridPar3, DIKIND)+1))

      do q=1, 2

         do w=1, 2

!     interpolate over parameter 3

            int_index = ((q+index1-2) * gridDim(2) + (w+index2-2)) * 
     &           gridDim(3) + index3

            slope = (dataField(int_index+1) - dataField(int_index)) /
     &           (gridPar3(index3+1) - gridPar3(index3))

            value3(w) = (input3 - gridPar3(index3)) * slope +
     &           dataField(int_index)

         enddo

!     interpolate over parameter 2

         slope = (value3(2) - value3(1)) / 
     &        (gridPar2(index2+1) - gridPar2(index2))

         value2(q) = (input2 - gridPar2(index2)) * slope + value3(1)

      enddo

!     interpolate over parameter 1

      slope = (value2(2) - value2(1)) /
     &     (gridPar1(index1+1) - gridPar1(index1))

      value = (input1 - gridPar1(index1)) * slope + value2(1)

      return
      end
!=======================================================================
!////////////////////  SUBROUTINE INTERPOLATE_3Dz  \\\\\\\\\\\\\\\\\\\\\
!
!     Similar to interpolate_3D except index2 is calculated
!     ahead of time because it is the redshift and will not 
!     change for the entire grid.
!
      subroutine interpolate_3Dz_g(
     &     input1, input2, input3, gridDim,
     &     gridPar1, dgridPar1,
     &     gridPar2, index2,
     &     gridPar3, dgridPar3,
     &     dataSize, dataField, 
     &     end_int, value)

      implicit NONE
#include "grackle_fortran_types.def"

!  General Arguments
      integer*8 dataSize, index2
      integer*8 gridDim(3)
      real*8 input1, input2, input3, value
      real*8 gridPar1(gridDim(1)), dgridPar1,
     &       gridPar2(gridDim(2)),
     &       gridPar3(gridDim(3)), dgridPar3
      real*8 dataField(dataSize)

!  Locals

      logical end_int
      integer*8 index1, index3, int_index
      integer q, w
      real*8 slope, value3(2), value2(2)

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
!=======================================================================

      if (end_int) then
         call interpolate_2Df3D_g(input1,
     &        input3, gridDim,
     &        gridPar1, dgridPar1,
     &        index2,
     &        gridPar3, dgridPar3,
     &        dataSize, dataField, 
     &        value)
         return
      endif

!     Calculate interpolation indices

      index1 = min(gridDim(1)-1,
     &     max(1,int((input1-gridPar1(1))/dgridPar1, DIKIND)+1))
      index3 = min(gridDim(3)-1,
     &     max(1,int((input3-gridPar3(1))/dgridPar3, DIKIND)+1))

      do q=1, 2

         do w=1, 2

!     interpolate over parameter 3

            int_index = ((q+index1-2) * gridDim(2) + (w+index2-2)) * 
     &           gridDim(3) + index3

            slope = (dataField(int_index+1) - dataField(int_index)) /
     &           (gridPar3(index3+1) - gridPar3(index3))

            value3(w) = (input3 - gridPar3(index3)) * slope +
     &           dataField(int_index)

         enddo

!     interpolate over parameter 2

         slope = (value3(2) - value3(1)) / 
     &        (gridPar2(index2+1) - gridPar2(index2))

         value2(q) = (input2 - gridPar2(index2)) * slope + value3(1)

      enddo

!     interpolate over parameter 1

      slope = (value2(2) - value2(1)) /
     &     (gridPar1(index1+1) - gridPar1(index1))

      value = (input1 - gridPar1(index1)) * slope + value2(1)

      return
      end

!=======================================================================
!///////////////////  SUBROUTINE INTERPOLATE_2Df3D  \\\\\\\\\\\\\\\\\\\\
!
!     Interpolation in 2 dimensions but with a 3D grid.
!     This is used for interpolating from just the last 
!     slice in the datacube before the redshift where 
!     the UV background turns on.
!
      subroutine interpolate_2Df3D_g(
     &     input1, input3, gridDim,
     &     gridPar1, dgridPar1,
     &     index2,
     &     gridPar3, dgridPar3,
     &     dataSize, dataField, 
     &     value)

      implicit NONE
#include "grackle_fortran_types.def"

!  General Arguments

      integer*8 dataSize, index2
      integer*8 gridDim(3)
      real*8 input1, input3, value
      real*8 gridPar1(gridDim(1)), dgridPar1,
     &       gridPar3(gridDim(3)), dgridPar3
      real*8 dataField(dataSize)

!  Locals

      integer*8 index1, index3, int_index
      integer q
      real*8 slope, value3(2)

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
!=======================================================================

!     Calculate interpolation indices

      index1 = min(gridDim(1)-1,
     &     max(1,int((input1-gridPar1(1))/dgridPar1, DIKIND)+1))
      index3 = min(gridDim(3)-1,
     &     max(1,int((input3-gridPar3(1))/dgridPar3, DIKIND)+1))

      do q=1, 2

!     interpolate over parameter 3

            int_index = ((q+index1-2) * gridDim(2) + (index2-1)) * 
     &           gridDim(3) + index3

            slope = (dataField(int_index+1) - dataField(int_index)) /
     &           (gridPar3(index3+1) - gridPar3(index3))

            value3(q) = (input3 - gridPar3(index3)) * slope +
     &           dataField(int_index)

      enddo

!     interpolate over parameter 1

      slope = (value3(2) - value3(1)) /
     &     (gridPar1(index1+1) - gridPar1(index1))

      value = (input1 - gridPar1(index1)) * slope + value3(1)

      return
      end

!=======================================================================
!////////////////////  SUBROUTINE INTERPOLATE_4D  \\\\\\\\\\\\\\\\\\\\\\

      subroutine interpolate_4D_g(
     &     input1, input2, input3, input4, 
     &     gridDim,
     &     gridPar1, dgridPar1,
     &     gridPar2, dgridPar2,
     &     gridPar3, dgridPar3,
     &     gridPar4, dgridPar4,
     &     dataSize, dataField, value)

      implicit NONE
#include "grackle_fortran_types.def"

!  General Arguments

      integer*8 dataSize
      integer*8 gridDim(4)
      real*8 input1, input2, input3, input4, value
      real*8 gridPar1(gridDim(1)), dgridPar1,
     &     gridPar2(gridDim(2)), dgridPar2,
     &     gridPar3(gridDim(3)), dgridPar3,
     &     gridPar4(gridDim(4)), dgridPar4
      real*8 dataField(dataSize)

!  Locals

      integer*8 index1, index2, index3, index4, int_index, q, w, e
      real*8 slope, value4(2), value3(2), value2(2)

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
!=======================================================================

!     Calculate interpolation indices

      index1 = min(gridDim(1)-1, max(1,
     &     int((input1-gridPar1(1))/dgridPar1,DIKIND)+1))
      index2 = min(gridDim(2)-1, max(1,
     &     int((input2-gridPar2(1))/dgridPar2,DIKIND)+1))
      index3 = min(gridDim(3)-1, max(1,
     &     int((input3-gridPar3(1))/dgridPar3,DIKIND)+1))
      index4 = min(gridDim(4)-1, max(1,
     &     int((input4-gridPar4(1))/dgridPar4,DIKIND)+1))

      do q=1, 2

         do w=1, 2

            do e=1, 2

!     interpolate over parameter 4

               int_index = (((q+index1-2) * gridDim(2) + (w+index2-2)) * 
     &              gridDim(3) + (e+index3-2)) * gridDim(4) + index4

               slope = (dataField(int_index+1) - dataField(int_index)) /
     &              (gridPar4(index4+1) - gridPar4(index4))

               value4(e) = (input4 - gridPar4(index4)) * slope + 
     &              dataField(int_index)

            enddo

!     interpolate over parameter 3

            slope = (value4(2) - value4(1)) /
     &           (gridPar3(index3+1) - gridPar3(index3))

            value3(w) = (input3 - gridPar3(index3)) * slope +
     &           value4(1)

         enddo

!     interpolate over parameter 2

         slope = (value3(2) - value3(1)) /
     &        (gridPar2(index2+1) - gridPar2(index2))

         value2(q) = (input2 - gridPar2(index2)) * slope + value3(1)

      enddo

!     interpolate over parameter 1

      slope = (value2(2) - value2(1)) /
     &     (gridPar1(index1+1) - gridPar1(index1))

      value = (input1 - gridPar1(index1)) * slope + value2(1)

      return
      end

!=======================================================================
!////////////////////  SUBROUTINE INTERPOLATE_5D  \\\\\\\\\\\\\\\\\\\\\\

      subroutine interpolate_5D_g(
     &     input1, input2, input3, input4, input5,
     &     gridDim,
     &     gridPar1, dgridPar1,
     &     gridPar2, dgridPar2,
     &     gridPar3, dgridPar3,
     &     gridPar4, dgridPar4,
     &     gridPar5, dgridPar5,
     &     dataSize, dataField, value)

      implicit NONE
#include "grackle_fortran_types.def"

!  General Arguments

      integer*8 dataSize
      integer*8 gridDim(5)
      real*8 input1, input2, input3, input4, input5, value
      real*8 gridPar1(gridDim(1)), dgridPar1,
     &     gridPar2(gridDim(2)), dgridPar2,
     &     gridPar3(gridDim(3)), dgridPar3,
     &     gridPar4(gridDim(4)), dgridPar4,
     &     gridPar5(gridDim(5)), dgridPar5
      real*8 dataField(dataSize)

!  Locals

      integer*8 index1, index2, index3, index4, index5, 
     &     int_index, q, w, e, r, midPt, highPt
      real*8 slope, value5(2), value4(2), value3(2), value2(2)

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
!=======================================================================

!     Calculate interpolation indices

      index1 = min(gridDim(1)-1, max(1,
     &     int((input1-gridPar1(1))/dgridPar1,DIKIND)+1))
      index2 = min(gridDim(2)-1, max(1,
     &     int((input2-gridPar2(1))/dgridPar2,DIKIND)+1))
      index3 = min(gridDim(3)-1, max(1,
     &     int((input3-gridPar3(1))/dgridPar3,DIKIND)+1))
#define INDEX_4_BISECTION
#ifdef INDEX_4_BISECTION
!     get index 4 with bisection, since not evenly spaced
      if (input4 <= gridPar4(1)) then
         index4 = 1
      else if (input4 >= gridPar4(gridDim(4)-1)) then
         index4 = gridDim(4) - 1
      else
         index4 = 1
         highPt = gridDim(4)
         do while ((highPt - index4) > 1)
            midPt = int((highPt + index4) / 2,DIKIND)
            if (input4 >= gridPar4(midPt)) then
               index4 = midPt
            else
               highPt = midPt
            endif
         enddo
      endif
#else
      index4 = min(gridDim(4)-1, max(1,
     &     int((input4-gridPar4(1))/dgridPar4,DIKIND)+1))
#endif /* INDEX_4_BISECTION */
      index5 = min(gridDim(5)-1, max(1,
     &     int((input5-gridPar5(1))/dgridPar5,DIKIND)+1))

      do q=1, 2

         do w=1, 2

            do e=1, 2

               do r=1, 2

!     interpolate over parameter 5

                  int_index = ((((q+index1-2) * gridDim(2) + 
     &                 (w+index2-2)) * gridDim(3) + (e+index3-2)) * 
     &                 gridDim(4) + (r+index4-2)) * gridDim(5) +
     &                 index5

                  slope = (dataField(int_index+1) - 
     &                 dataField(int_index)) /
     &                 (gridPar5(index5+1) - gridPar5(index5))

                  value5(r) = (input5 - gridPar5(index5)) * slope +
     &                 dataField(int_index)

               enddo

!     interpolate over parameter 4

               slope = (value5(2) - value5(1)) /
     &              (gridPar4(index4+1) - gridPar4(index4))

               value4(e) = (input4 - gridPar4(index4)) * slope +
     &              value5(1)

            enddo

!     interpolate over parameter 3

            slope = (value4(2) - value4(1)) /
     &           (gridPar3(index3+1) - gridPar3(index3))

            value3(w) = (input3 - gridPar3(index3)) * slope +
     &           value4(1)

         enddo

!     interpolate over parameter 2

         slope = (value3(2) - value3(1)) /
     &        (gridPar2(index2+1) - gridPar2(index2))

         value2(q) = (input2 - gridPar2(index2)) * slope +
     &        value3(1)

      enddo

!     interpolate over parameter 1

      slope = (value2(2) - value2(1)) /
     &     (gridPar1(index1+1) - gridPar1(index1))

      value = (input1 - gridPar1(index1)) * slope + value2(1)

      return
      end
