      subroutine soft_nab_qgq(order,y12,y23,y31,
     & preIijm,j1,j2,j3,j4,j5,j6,soft2)
      implicit none
!     non-abelian piece  q(1) g(2) q(3)
!     returns coefficients of [log(tau)^j/tau_+],and for j=-1, delta(tau)
!     in units of as/2/pi
!     y12=n1.n2/2 so that 0<y12<1
!     y23=n2.n3/2 so that 0<y23<1
!     y31=n3.n1/2 so that 0<y31<1
!     order indicates the order in [as/2/pi]
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'scet_const.f'
      real(dp)::y12,y23,y31,soft2(-1:3),preIijm(5),
     & Iijm(6),Lss13,Lss23
      integer::order,j1,j2,j3,j4,j5,j6
!      integer,parameter::i(6)=(/1,2,2,3,3,1/)
!      integer,parameter::j(6)=(/2,1,3,2,1,3/)
!      integer,parameter::m(6)=(/3,3,1,1,2,2/)
!Id,I12x3=Iijm(1);
!Id,I21x3=Iijm(2);
!Id,I23x1=Iijm(3);
!Id,I32x1=Iijm(4);
!Id,I31x2=Iijm(5);
!Id,I13x2=Iijm(6);

      soft2(:)=zip
      if (order < 2) return

!      y(:,:)=zip

!      y(1,2)=y12
!      y(2,1)=y12
!      y(2,3)=y23
!      y(3,2)=y23
!      y(1,3)=y31
!      y(3,1)=y31

!      do n=1,6
!      al=y(j(n),m(n))/y(i(n),j(n))
!      be=y(i(n),m(n))/y(i(n),j(n))
!      I0=I0JSTW(al,be)
!      I1=I1JSTW(al,be)
!      Iijm(n)=I0*log(al)+I1
!      enddo

      Iijm(1)=preIijm(j1)
      Iijm(2)=preIijm(j2)
      Iijm(3)=preIijm(j3)
      Iijm(4)=preIijm(j4)
      Iijm(5)=preIijm(j5)
      Iijm(6)=preIijm(j6)

      soft2(3)=zip
      soft2(2)=2*be0*(CA+2*CF)
      soft2(1)=+CA*(-half*Ga1
     & -two*log(y12)*be0+two*log(y31)*be0-two*log(y23)*be0)
     & +CF*(-Ga1-four*log(y31)*be0)
      soft2(0)=CA*(
     &     - quarter*gams1
     &     - half*zeta2*be0
     &     + quarter*log(y12)*Ga1
     &     + half*log(y12)**2*be0
     &     - quarter*log(y31)*Ga1
     &     - half*log(y31)**2*be0
     &     + quarter*log(y23)*Ga1
     &     + half*log(y23)**2*be0
     &     + Iijm(1)*be0
     &     + Iijm(2)*be0
     &     + Iijm(3)*be0
     &     + Iijm(4)*be0
     &     - Iijm(5)*be0
     &     - Iijm(6)*be0
     &     )
      soft2(0) = soft2(0) + CF * (
     &     - half*gams1
     &     - zeta2*be0
     &     + half*log(y31)*Ga1
     &     + log(y31)**2*be0
     &     + two*Iijm(5)*be0
     &     + two*Iijm(6)*be0
     &     )

      Lss13=log(y31)
      Lss23=log(y23)
      soft2(-1)=-200.096_dp-223.007_dp*Lss13-210.796_dp*Lss23
     & -112.644_dp*Lss13**2-118.705_dp*Lss23**2
     & -30.8755_dp*Lss13**3-32.7046_dp*Lss23**3
     & -5.57982_dp*Lss13**4-5.43977_dp*Lss23**4
     & -0.666669_dp*Lss13**5-0.63292_dp*Lss23**5
     & -0.0533767_dp*Lss13**6-0.0484011_dp*Lss23**6
     & -0.0028547_dp*Lss13**7-0.00240717_dp*Lss23**7
     & -0.0000952163_dp*Lss13**8-0.0000728699_dp*Lss23**8
     & -1.45696e-6_dp*Lss13**9-9.83475e-7_dp*Lss23**9

!      soft2(-1)=soft2(-1)/four ! DEBUG

      return
      end
