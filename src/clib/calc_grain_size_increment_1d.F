#include "phys_const.def"

!=======================================================================
!//////////////  CALC_GRAIN_SIZE_INCREMENT_SPECIES_1D \\\\\\\\\\\\\\\\\\

      subroutine calc_grain_size_increment_species_1d(
     &          igrgr, itmask, SN0_N
     &        , in, jn, kn, is, ie, j, k, dom, d
     &        , nSN, dsp, SN_metal, SN_fsp, SN_r0sp
     &        , ssp, sgsp, alsp
     &        , gr_N, gr_Size, gr_dT, gr_Td, SN_kp0sp
     &            )

!
!  Calculates properties that are derived from the grain size increment
!  for a single grain species.
!
!  NOTE: the subroutine's name should be more descriptive
!
!  INPUTS:
!    ssp - the density of an individual dust grain species. At the time
!          of writing it has units of g/cm^3
!
!
      implicit NONE

#include "grackle_fortran_types.def"

! input
      integer in, jn, kn, is, ie, j, k
      MASK_TYPE itmask(in)
      integer igrgr
      integer iSN, nSN, SN0_N
      real*8  dom
      R_PREC  d(in,jn,kn)
      R_PREC  dsp(in,jn,kn), SN_metal(in, SN0_N)
! table
      real*8  SN_fsp(SN0_N), SN_r0sp(3,SN0_N)
      real*8  ssp
! opacity table
      integer gr_N(2), gr_Size
      real*8  gr_dT, gr_Td(gr_N(2))
      real*8  SN_kp0sp(gr_Size,SN0_N)
! output
      real*8 sgsp(in), alsp(gr_N(2),in)
! local
      integer i
      real*8 coef0, coef1, coef2, coef3
      real*8 dsp_inject_sum
      real*8 dsp0, SN_sgsp, SN_kpsp
      real*8 SN_dsp0(SN0_N), SN_nsp0(SN0_N)
      real*8 drsp(in)
      real*8, parameter :: pi = pi_val
      real*8, parameter :: mh = mass_h
      integer iTd, iTd0
 
      do i = is+1, ie+1
         if ( itmask(i) .ne. MASK_FALSE ) then

! Step 1: compute the total mass density of the current grain species
!         that was injected (by summing the amounts injected by each
!         injection pathway)
            do iSN = 1, nSN
               if(SN_fsp(iSN) .gt. 0.d0) then
                  SN_dsp0(iSN) = SN_fsp(iSN) * SN_metal(i, iSN)
               endif
            enddo

! Step 2: Compute the size increment for the grain species
!
! Let's go into detail:
! - todo: maybe we move most of this description into the docstring
!   and possibly move some of it into the narrative documentation
! - let's define some notation:
!    - for some time t, let φ(r,t) denote the grain species's
!      differential size distribution at time t, ρ(t) denote the grain
!      species's mass density, and n(t) denote it number density.
!    - we define <rᵖ>(t) as the `p`th order moment of `φ(r,t)`. In
!      other words, <rᵖ>ⱼ = ∫ rᵖ φⱼ dr, where the integral is taken
!      from 0 to infinity
!    - for additional context:
!      - φ(r,t) is normalized such that the p=0 moment is 1
!      - n(t)*φ(r,t) specifies the number density of the grain with
!        radii between r and r+dr
!      - given the mass density of a single grain, ζ, the variables
!        are related via ρ(t) = (4π/3) ζ n(t) <r³>(t).
!    - we define Φⱼ(r) as the grain species's initial differential
!      size distribution when distributed by injection pathway j.
!    - for notational convenience, we define <rᵖ>ⱼ as the `p`th order
!      moment of `Φⱼ(r)`
!    - NOTE: we don't explicitly track φ(r,t) or Φⱼ(r), instead we
!      track moments
! - The underlying model assumes that all injection pathways inject
!   the grain species at a single time t_inj.
!   - In other words:  φ(r,t_inj) = (∑ⱼ nⱼ Φⱼ(r)) / (∑ⱼ nⱼ)
!   - Consequently:
!         <rᵖ>(t_inj) = (∑ⱼ nⱼ <rᵖ>ⱼ) / (∑ⱼ nⱼ)
!         OR
!         <rᵖ>(t_inj) = (∑ⱼ <rᵖ>ⱼ ρⱼ/<r³>ⱼ) / (∑ⱼ ρⱼ/<r³>ⱼ)
! - Furthermore we assume that:
!   - grains can only be created via injection (i.e. the number
!     density of the grain species is frozen)
!   - φ(r,t) can't be deformed. However it can be translated (which
!     is exactly what happens when grains undergo growth).
! - Thus, when we model growth, φ(r,t) = φ(r - δr(t), t_inj) where
!   we call δr(t) the "size increment"
!
! Anybody reading this might notice that several limitations to this
! model. This is discussed in detail in GH Issue 444 (the description
! should be made part of the narrative docs

            if(igrgr .eq. 0) then
! the model is constructed such that δr(t) = 0 if we aren't
! modelling grain growth

            drsp(i) = 0.d0

            else

! here we calculate the current grain size increment, δr(t), using
! conservation of grain number (i.e. grains are only created via
! injection pathways)
! - recall that to model growth we effectively assume that each
!   occurrence of a grain species was injected at a single time, t_inj
! - conservation of grain number tells us that
!   <r³>(t) = <r³>(t_inj) * ρ(t) / ρ(t_inj),
!   where ρ(t_inj) = ∑ⱼ ρⱼ
! - Equation A3 of Chiaki & Wise 2019 indicates that the grain-size
!   increment is the root of the following cubic equation
!    [<r³>(t_inj) - <r³>(t)] + [3 * <r²>(t_inj)] * δr(t)
!       + [3 * <r¹>(t_inj)] * δr²(t) + δr³(t) = 0
! - Let's rewrite the expression (<r³>(t_inj) - <r³>(t)). First,
!   let's plug in the equation from conservation of grain number:
!     <r³>(t_inj) - <r³>(t) = <r³>(t_inj) * [ 1 - (ρ(t)/∑ⱼρⱼ)]
!   Now, let's plug in <r³>(t_inj) = (∑ⱼρⱼ) / (∑ⱼ ρⱼ/<r³>ⱼ), which
!   comes from the generic formula for <rᵖ>(t_inj) defined earlier:
!     <r³>(t_inj) - <r³>(t) = [(∑ⱼρⱼ) - ρ(t)] / (∑ⱼ ρⱼ/<r³>ⱼ)
! - If we plug this into the preceding cubic equation and multiply the
!   resulting equation by (∑ⱼ ρⱼ/<r³>ⱼ), we get:
!       coef0 + coef1 * δr(t) + coef2 * δr²(t) + coef3 * δr³(t) = 0
!   where:
!       - coef0 = (∑ⱼ ρⱼ) - ρ(t)
!       - coef1 = ∑ⱼ 3 <r²>ⱼ ρⱼ/<r³>ⱼ
!       - coef2 = ∑ⱼ 3 <r¹>ⱼ ρⱼ/<r³>ⱼ
!       - coef3 = ∑ⱼ ρⱼ/<r³>ⱼ

! Let's compute these coefficients
            dsp_inject_sum = 0.d0
! We'll use dsp_inject_sum to compute coef0
            coef1 = 0.d0
            coef2 = 0.d0
            coef3 = 0.d0

!  Loop over each injection pathway
            do iSN = 1, nSN
               if(SN_fsp(iSN) .gt. 0.d0) then
!  Calculate 4πζnⱼ/3 = ρⱼ/<r³>ⱼ
!  -> recall: that ζ is the mass density of a single grain of the
!     current grain species (i.e. it's a constant)
!  -> TODO: its very confusing that we store the result within SN_nsp0
!     since that variable is later reused to directly hold nⱼ
                  SN_nsp0(iSN) = SN_dsp0(iSN) / SN_r0sp(3,iSN)

                  dsp_inject_sum = dsp_inject_sum + SN_dsp0(iSN)
                  coef1 = coef1 + 3.d0 * SN_nsp0(iSN) * SN_r0sp(2,iSN)
                  coef2 = coef2 + 3.d0 * SN_nsp0(iSN) * SN_r0sp(1,iSN)
                  coef3 = coef3 +        SN_nsp0(iSN)
               endif
            enddo

            coef0 = dsp_inject_sum - dsp(i,j,k)

            if(coef3 .eq. 0.d0 .or. coef3 .ne. coef3) then
               if(i .ne. 1 .or. j .ne. 1 .or. k .ne. 1) then
                  write(*,*) "Warning: coef3 is zero or NaN for i=", i, "j=", j, "k=", k, "coef3=", coef3
               endif
            endif
!  Let's actually solve for the δr(t), the root of the cubic equation
!  - before we do that, we divide both sides of the equation by the
!    value of coef3. This ensures that the rescaled value of coef3 is
!    exactly 1 (a requirement of the solve_cubic_equation routine)
!  - todo: consider not modifying the variables in-place or performing
!    this operation within the solve_cubic_equation
            coef0 = coef0 / coef3
            coef1 = coef1 / coef3
            coef2 = coef2 / coef3

            call solve_cubic_equation(coef2, coef1, coef0, drsp(i))
            if(i .eq. 1 .and. j .eq. 1 .and. k .eq. 1) then
               drsp(i) = 0.d0
            endif

            drsp(i) = max(drsp(i), 0.d0)

            endif


! Step 3: calculate number density (code_density / g)

            do iSN = 1, nSN
               if(SN_fsp(iSN) .gt. 0.d0) then
                  SN_nsp0(iSN) = SN_dsp0(iSN)
     &             / (4.d0*pi/3.d0 * ssp * SN_r0sp(3,iSN))
               else
                  SN_nsp0(iSN) = 0.d0
               endif
            enddo

! Step 4: calculate geometrical cross-section per unit gas mass
! -> units of cm^2/g
            sgsp(i) = 0.d0
            do iSN = 1, nSN
               if( SN_fsp(iSN) .gt. 0.d0) then
                  SN_sgsp = pi *
     &               (        SN_r0sp(2,iSN)
     &               + 2.d0 * SN_r0sp(1,iSN) * drsp(i)
     &               +                         drsp(i)**2
     &               )
               else
                  SN_sgsp = 0.d0
               endif
               sgsp(i) = sgsp(i) + SN_nsp0(iSN) * SN_sgsp
            enddo
            sgsp(i) = sgsp(i) / d(i,j,k)
    
! Step 5: calculate optical opacity related quantities
! -> we are effectively constructing a 1d table of values, at various
!    possible grain temperature for each injection pathway (in other
!    words, its a 2D table)
! -> I'm pretty sure that the values we are the opacity (per unit gas
!    mass)
! -> I'm pretty confident that the units are independent of code units
!    (I think the units are cm^2/g)
            do iTd = 1, gr_N(2)
               iTd0 = (iTd-1)*gr_N(1)
               alsp(iTd,i) = 0.d0
               do iSN = 1, nSN
                  if( SN_fsp(iSN) .gt. 0.d0) then
                     SN_kpsp = 4.d0*pi/3.d0 * ssp *
     &                  (        SN_kp0sp(iTd0+4,iSN)
     &                  + 3.d0 * SN_kp0sp(iTd0+3,iSN) * drsp(i)
     &                  + 3.d0 * SN_kp0sp(iTd0+2,iSN) * drsp(i)**2
     &                  +        SN_kp0sp(iTd0+1,iSN) * drsp(i)**3
     &                  )
                  else
                     SN_kpsp = 0.d0
                  endif
                  alsp(iTd,i) = alsp(iTd,i) + SN_nsp0(iSN) * SN_kpsp
               enddo
               alsp(iTd,i) = alsp(iTd,i) / d(i,j,k)
            enddo


         endif ! itmask
      enddo

      return
      end

!=======================================================================
!////////////////  SUBROUTINE SOLVE_CUBIC_EQUATION \\\\\\\\\\\\\\\\\\\\\

      subroutine solve_cubic_equation(a, b, c, root)

      implicit NONE
#include "grackle_fortran_types.def"

!  General Arguments

      real*8 a, b, c
      real*8 root
      real*8 root1, root2, root3

!  Locals

      real*8 q, r, m
      real*8 th
      real*8 s,t
      real*8, parameter :: pi = pi_val

      q = (a*a - 3.d0*b)/9.d0
      r = (2.d0*a*a*a - 9.d0*a*b + 27.d0*c)/54.d0
      m = r*r - q*q*q
    
      if (m .lt. 0.d0) then !! three real roots 
        th = acos( r / sqrt(q*q*q) )
        root1 = -(2.d0*sqrt(q)*cos( th         /3.d0))-a/3.d0
        root2 = -(2.d0*sqrt(q)*cos((th+2.d0*pi)/3.d0))-a/3.d0
        root3 = -(2.d0*sqrt(q)*cos((th-2.d0*pi)/3.d0))-a/3.d0
        write(*,*) "three real roots", root1, root2, root3
        stop
      else !! one real root
        if(r .gt. 0.d0) then
          s = -( r+sqrt(m))**(1.d0/3.d0)
        else
          s =  (-r+sqrt(m))**(1.d0/3.d0)
        endif
        t = q/s
        root = s+t-a/3.d0
      endif

      return
      end
