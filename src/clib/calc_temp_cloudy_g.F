#include "phys_const.def"

! This file originally defined the calc_temp_cloudy_g subroutine and a
! helper subroutine (scale_fields_table_g)
! - calc_temp_cloudy_g has now been transcribed to C++
! - scale_fields_table_g still need to be transcribed

c -----------------------------------------------------------
!   This routine scales the density fields from comoving to
!     proper densities (and back again).

      subroutine scale_fields_table_g(d, metal,
     &                        is, ie, js, je, ks, ke,
     &                        in, jn, kn, imetal, factor)
c -------------------------------------------------------------------

      implicit NONE
#include "grackle_fortran_types.def"

!     Arguments

      integer in, jn, kn, is, ie, js, je, ks, ke, imetal
      R_PREC  d(in,jn,kn), metal(in,jn,kn)
      real*8 factor

!     locals

      integer i, j, k

!     Multiply all fields by factor (1/a^3 or a^3)

      do k = ks+1, ke+1
         do j = js+1, je+1
            do i = is+1, ie+1
               d(i,j,k)        = d(i,j,k) * factor
            enddo
         enddo
      enddo

      if (imetal .eq. 1) then
         do k = ks+1, ke+1
            do j = js+1, je+1
               do i = is+1, ie+1
                  metal(i,j,k) = metal(i,j,k) * factor
               enddo
            enddo
         enddo
      endif

      return
      end
