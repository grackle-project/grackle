!=======================================================================
!////////////////////  SUBROUTINE SCALE_FIELDS_G  \\\\\\\\\\\\\\\\\\\\\\

      subroutine my_scale_fields_g(
     &     ispecies, imetal, idustfield, imchem,
     &     imabund, idspecies, immulti, igrgr, idsub, factor,
     &     is, ie, js, je, ks, ke, in, jn, kn,
     &     d, de, HI, HII, HeI, HeII, HeIII,
     &     HM, H2I, H2II, DI, DII, HDI, DM, HDII, HeHII,
     &     metal, dust,
     &     CI, CII, CO, CO2,
     &     OI, OH, H2O, O2,
     &     SiI, SiOI, SiO2I,
     &     CH, CH2, COII, OII,
     &     OHII, H2OII, H3OII, O2II,
     &     Mg, Al, S, Fe,
     &     SiM, FeM, Mg2SiO4, MgSiO3, Fe3O4,
     &     AC, SiO2D, MgO, FeS, Al2O3,
     &     reforg, volorg, H2Oice,
     &     metal_loc, metal_C13, metal_C20, metal_C25, metal_C30,
     &     metal_F13, metal_F15, metal_F50, metal_F80,
     &     metal_P170, metal_P200, metal_Y19)
c -------------------------------------------------------------------

      implicit NONE
#include "grackle_fortran_types.def"

      integer ispecies, imetal, idustfield, imchem,
     &     imabund, idspecies, immulti, igrgr, idsub,
     &     is, ie, js, je, ks, ke, in, jn, kn

      R_PREC factor

      R_PREC d(in,jn,kn), de(in,jn,kn), HI(in,jn,kn), HII(in,jn,kn),
     &     HeI(in,jn,kn), HeII(in,jn,kn), HeIII(in,jn,kn),
     &     HM(in,jn,kn), H2I(in,jn,kn), H2II(in,jn,kn),
     &     DI(in,jn,kn), DII(in,jn,kn), HDI(in,jn,kn),
     &     DM(in,jn,kn), HDII(in,jn,kn), HeHII(in,jn,kn),
     &     metal(in,jn,kn), dust(in,jn,kn),
     &     CI(in,jn,kn), CII(in,jn,kn), CO(in,jn,kn), CO2(in,jn,kn),
     &     OI(in,jn,kn), OH(in,jn,kn), H2O(in,jn,kn), O2(in,jn,kn),
     &     SiI(in,jn,kn), SiOI(in,jn,kn), SiO2I(in,jn,kn),
     &     CH(in,jn,kn), CH2(in,jn,kn), COII(in,jn,kn), OII(in,jn,kn),
     &     OHII(in,jn,kn), H2OII(in,jn,kn), H3OII(in,jn,kn),
     &     O2II(in,jn,kn),
     &     Mg(in,jn,kn), Al(in,jn,kn), S(in,jn,kn), Fe(in,jn,kn),
     &     SiM(in,jn,kn), FeM(in,jn,kn), Mg2SiO4(in,jn,kn),
     &     MgSiO3(in,jn,kn), Fe3O4(in,jn,kn),
     &     AC(in,jn,kn), SiO2D(in,jn,kn), MgO(in,jn,kn),
     &     FeS(in,jn,kn), Al2O3(in,jn,kn),
     &     reforg(in,jn,kn), volorg(in,jn,kn), H2Oice(in,jn,kn),
     &     metal_loc(in,jn,kn), metal_C13(in,jn,kn),
     &     metal_C20(in,jn,kn), metal_C25(in,jn,kn),
     &     metal_C30(in,jn,kn),
     &     metal_F13(in,jn,kn), metal_F15(in,jn,kn),
     &     metal_F50(in,jn,kn), metal_F80(in,jn,kn),
     &     metal_P170(in,jn,kn), metal_P200(in,jn,kn),
     &     metal_Y19(in,jn,kn)

      integer i, j, k, t, dj, dk

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
!=======================================================================
      dk = ke - ks + 1
      dj = je - js + 1

! parallelize the k and j loops with OpenMP
! flat j and k loops for better parallelism
#ifdef _OPENMP
!$omp parallel do schedule(runtime) private(
!$omp&   i, j, k )
#endif
      do t = 0, dk*dj-1
         k = t/dj      + ks+1
         j = mod(t,dj) + js+1

         do i = is+1, ie+1
            d(i,j,k) = d(i,j,k) * factor
         enddo

         if (ispecies .gt. 0) then
            do i = is+1, ie+1
               de(i,j,k)    = de(i,j,k)    * factor
               HI(i,j,k)    = HI(i,j,k)    * factor
               HII(i,j,k)   = HII(i,j,k)   * factor
               HeI(i,j,k)   = HeI(i,j,k)   * factor
               HeII(i,j,k)  = HeII(i,j,k)  * factor
               HeIII(i,j,k) = HeIII(i,j,k) * factor
            enddo
         endif
         if (ispecies .gt. 1) then
            do i = is+1, ie+1
               HM(i,j,k)   = HM(i,j,k)   * factor
               H2I(i,j,k)  = H2I(i,j,k)  * factor
               H2II(i,j,k) = H2II(i,j,k) * factor
            enddo
         endif
         if (ispecies .gt. 2) then
            do i = is+1, ie+1
               DI(i,j,k)  = DI(i,j,k)  * factor
               DII(i,j,k) = DII(i,j,k) * factor
               HDI(i,j,k) = HDI(i,j,k) * factor
            enddo
         endif
         if (ispecies .gt. 3) then
           do i = is+1, ie+1
              DM(i,j,k)    = DM(i,j,k)    * factor
              HDII(i,j,k)  = HDII(i,j,k)  * factor
              HeHII(i,j,k) = HeHII(i,j,k) * factor
           enddo
        endif

        if (imetal .eq. 1) then
           do i = is+1, ie+1
              metal(i,j,k) = metal(i,j,k) * factor
           enddo

           if (imchem .eq. 1) then
              do i = is+1, ie+1
                 CI(i,j,k)    = CI(i,j,k)    * factor
                 CII(i,j,k)   = CII(i,j,k)   * factor
                 CO(i,j,k)    = CO(i,j,k)    * factor
                 CO2(i,j,k)   = CO2(i,j,k)   * factor
                 OI(i,j,k)    = OI(i,j,k)    * factor
                 OH(i,j,k)    = OH(i,j,k)    * factor
                 H2O(i,j,k)   = H2O(i,j,k)   * factor
                 O2(i,j,k)    = O2(i,j,k)    * factor
                 SiI(i,j,k)   = SiI(i,j,k)   * factor
                 SiOI(i,j,k)  = SiOI(i,j,k)  * factor
                 SiO2I(i,j,k) = SiO2I(i,j,k) * factor
                 CH(i,j,k)    = CH(i,j,k)    * factor
                 CH2(i,j,k)   = CH2(i,j,k)   * factor
                 COII(i,j,k)  = COII(i,j,k)  * factor
                 OII(i,j,k)   = OII(i,j,k)   * factor
                 OHII(i,j,k)  = OHII(i,j,k)  * factor
                 H2OII(i,j,k) = H2OII(i,j,k) * factor
                 H3OII(i,j,k) = H3OII(i,j,k) * factor
                 O2II(i,j,k)  = O2II(i,j,k)  * factor
              enddo
           endif

           if (immulti .gt. 0) then
              do i = is+1, ie+1
                 metal_loc(i,j,k)  = metal_loc(i,j,k)  * factor
                 metal_C13(i,j,k)  = metal_C13(i,j,k)  * factor
                 metal_C20(i,j,k)  = metal_C20(i,j,k)  * factor
                 metal_C25(i,j,k)  = metal_C25(i,j,k)  * factor
                 metal_C30(i,j,k)  = metal_C30(i,j,k)  * factor
                 metal_F13(i,j,k)  = metal_F13(i,j,k)  * factor
                 metal_F15(i,j,k)  = metal_F15(i,j,k)  * factor
                 metal_F50(i,j,k)  = metal_F50(i,j,k)  * factor
                 metal_F80(i,j,k)  = metal_F80(i,j,k)  * factor
                 metal_P170(i,j,k) = metal_P170(i,j,k) * factor
                 metal_P200(i,j,k) = metal_P200(i,j,k) * factor
                 metal_Y19(i,j,k)  = metal_Y19(i,j,k)  * factor
              enddo
           endif

           if ( ( igrgr .eq. 1 ) .or. ( idsub .eq. 1 ) ) then
              if (idspecies .gt. 0) then
                 do i = is+1, ie+1
                    Mg(i,j,k) = Mg(i,j,k) * factor
                 enddo
              endif
              if (idspecies .gt. 1) then
                 do i = is+1, ie+1
                    Al(i,j,k) = Al(i,j,k) * factor
                    S(i,j,k)  = S(i,j,k)  * factor
                    Fe(i,j,k) = Fe(i,j,k) * factor
                 enddo
              endif
           endif

           if (idustfield .eq. 1) then
              do i = is+1, ie+1
                 dust(i,j,k) = dust(i,j,k) * factor
              enddo

              if ( ( igrgr .eq. 1 ) .or. ( idsub .eq. 1 ) ) then
                 if (idspecies .gt. 0) then
                    do i = is+1, ie+1
                       MgSiO3(i,j,k) = MgSiO3(i,j,k) * factor
                       AC(i,j,k)     = AC(i,j,k)     * factor
                    enddo
                 endif
                 if (idspecies .gt. 1) then
                    do i = is+1, ie+1
                       SiM(i,j,k)     = SiM(i,j,k)     * factor
                       FeM(i,j,k)     = FeM(i,j,k)     * factor
                       Mg2SiO4(i,j,k) = Mg2SiO4(i,j,k) * factor
                       Fe3O4(i,j,k)   = Fe3O4(i,j,k)   * factor
                       SiO2D(i,j,k)   = SiO2D(i,j,k)   * factor
                       MgO(i,j,k)     = MgO(i,j,k)     * factor
                       FeS(i,j,k)     = FeS(i,j,k)     * factor
                       Al2O3(i,j,k)   = Al2O3(i,j,k)   * factor
                    enddo
                 endif
                 if (idspecies .gt. 2) then
                    do i = is+1, ie+1
                       reforg(i,j,k) = reforg(i,j,k) * factor
                       volorg(i,j,k) = volorg(i,j,k) * factor
                       H2Oice(i,j,k) = H2Oice(i,j,k) * factor
                    enddo
                 endif
              endif

           endif  ! if (idustfield .eq. 1)

        endif ! if (imetal .eq. 1)

      enddo
#ifdef _OPENMP
!$omp end parallel do
#endif

      return
      end
