!#################################################################
!##### O(as) and O(as^2) soft contribution to e+e- -> qqb type of
!##### processes. Calculates soft1(-1:1) and soft2(-1:3), i.e.
!##### contributions to delta(tau1), L0(tau1), L1(tau1), and so on.
!##### Factor of alphas/(2*pi) and (alphas/(2*pi))**2 have been 
!##### extracted out. Recovers softqqbis.f when L12->0.
!##### L12 is defined as Log[n1*n2/2] with n1 and n2 unit vectors
!##### of hard final-state jets.
!#################################################################
!---  rotuine by R. Mondini, edited by CW.
      
      subroutine soft_qqdec(order,soft1,soft2,L12)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'scale.f'
      include 'nf.f'
      include 'scet_const.f'
      integer :: order
      real(dp) :: soft1(-1:1),soft2(-1:3)
      real(dp) L12,be0NF 
      logical nfonly
      common/nfonly/nfonly
!$omp threadprivate(/nfonly/)

c---  evaluation of L12

      soft1(:)=zip
      soft2(:)=zip
      
      if (order < 1) return

c--- coefficient of (alphas/2/pi)

!---- contribution to delta(tau)
      soft1(-1) = CF*(zeta2-L12**2)
!---- contribution to L0(tau)
      soft1(0) = 4*CF*L12
!---- contribution to L1(tau)
      soft1(1) =-8*CF

      if (order < 2) return

c--- coefficient of (alphas/2/pi)^2

!---- contribution to delta(tau)
      soft2(-1) =
     & +CF**2*(-27/ten*zeta2**2 +L12**4/two 
     &         -9*zeta2*L12**2 -32*zeta3*L12)
     & +CF*(-5/27._dp*be0 - 160/27._dp*CA - 37/12._dp*be0*zeta2
     &      +2*CA*zeta2 + 22/five*CA*zeta2**2 + 29/six*be0*zeta3
     &      -14/nine*be0*L12 -16/nine*CA*L12 -5/six*be0*L12**2
     &      -2/three*CA*L12**2 -be0/six*L12**3 +be0*zeta2*L12
     &      +CA*zeta2*L12**2 +7*CA*zeta3*L12)
!---- contribution to L0(tau)
      soft2(0) = 
     & + CF**2*(64*zeta3 -4*L12**3 +36*zeta2*L12)
     & + CF*(28/nine*be0 - 2*be0*zeta2 + 32/nine*CA - 14*CA*zeta3
     &       +10/three*be0*L12 +8/three*CA*L12 +be0*L12**2 -4*CA*zeta2*L12)
!---- contribution to L1(tau)
      soft2(1) = 
     & + CF**2*(-72*zeta2 +24*L12**2)
     & + CF*(- 20/three*be0 - 16/three*CA + 8*CA*zeta2 -4*be0*L12)
!---- contribution to L2(tau)
      soft2(2) =
     & + CF**2*(-48*L12)
     & + 4*CF*be0
!---- contribution to L3(tau)
      soft2(3) = 32*CF**2

!===== uncomment this bit for nf only
!      if(nfonly) then
!          be0NF = -4/three*TR*NF
!          soft1(:)=zip
!          soft2(-1)=
!     &        +CF*(-5/27._dp*be0NF  - 37/12._dp*be0NF*zeta2
!     &        + 29/six*be0NF*zeta3
!     &        -14/nine*be0NF*L12  -5/six*be0NF*L12**2
!     &        -be0NF/six*L12**3 +be0NF*zeta2*L12)
!         soft2(0) = 
!     &        + CF*(28/nine*be0NF - 2*be0NF*zeta2
!     &        +10/three*be0NF*L12 +be0NF*L12**2)
!         soft2(1) = 
!     & + CF*(- 20/three*be0NF -4*be0NF*L12)
!         soft2(2) =4*CF*be0NF
!         soft2(3) =zip
!      endif

!======uncomment this bit for no nf
!         be0NF = -4/three*TR*NF
!         soft2(-1)=soft2(-1)-
!     &        ( +CF*(-5/27._dp*be0NF  - 37/12._dp*be0NF*zeta2
!     &        + 29/six*be0NF*zeta3
!     &        -14/nine*be0NF*L12  -5/six*be0NF*L12**2
!     &        -be0NF/six*L12**3 +be0NF*zeta2*L12))
!         soft2(0) = soft2(0)-(
!     &        + CF*(28/nine*be0NF - 2*be0NF*zeta2
!     &        +10/three*be0NF*L12 +be0NF*L12**2))
!         soft2(1) = soft2(1)-(
!     & + CF*(- 20/three*be0NF -4*be0NF*L12))
!        soft2(2) =soft2(2)-4*CF*be0NF
      
      
      return
      end


