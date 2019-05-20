!#################################################################
!##### O(as) and O(as^2) soft contribution to qg->q type of
!##### processes. Calculates soft1(-1:1) and soft2(-1:3), i.e.
!##### contributions to delta(tau1), L0(tau1), L1(tau1), and so on.
!##### Factor of alphas/(2*pi) and (alphas/(2*pi))**2 have been 
!##### extracted out.
!#################################################################
!---  RM, May 19
      
      subroutine soft_qgq(order,y12,y31,y23,soft1,soft2)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'scale.f'
      include 'nf.f'
      include 'scet_const.f'
      integer,intent(in)::order
      real(dp),intent(in)::y12,y31,y23
      real(dp),intent(out)::soft1(-1:1),soft2(-1:3)
      real(dp)::Iijm(6),Lss13,Lss23
      real(dp)::soft2ab(-1:3),soft2nab(-1:3)
      real(dp)::c00,c01,c02,c03,c10,c20,c30

      soft1(:)=zip
      soft2(:)=zip
      soft2ab(:)=zip
      soft2nab(:)=zip

      if (order < 1) return

      call computeIijmnobug(y12,y23,y31,Iijm)

!---- contribution to delta(tau)
      soft1(-1)=CA*(half*zeta2
     & -half*log(y12)**2+half*log(y31)**2-half*log(y23)**2
     & -Iijm(1)-Iijm(2)-Iijm(3)-Iijm(4)+Iijm(5)+Iijm(6))
     & +CF*(zeta2-log(y31)**2-2*Iijm(5)-2*Iijm(6))
!---- contribution to L0(tau)
      soft1(0)=2*CA*(log(y12)+log(y23)-log(y31))+4*log(y31)*CF
!---- contribution to L1(tau)
      soft1(1)=-4*CA-8*CF

      if (order < 2) return

!---- abelian contribution
      soft2ab(3)=8*(CA+2*CF)**2
      soft2ab(2)=-12*(CA+2*CF)*(CA*(log(y12)+log(y23)-log(y31))
     & +2*log(y31)*CF)

      soft2ab(1) =  + CA**2 * (
     &     - 18._dp*zeta2
     &     + 6._dp*log(y12)**2
     &     - 8._dp*log(y12)*log(y31)
     &     + 8._dp*log(y12)*log(y23)
     &     + 2._dp*log(y31)**2
     &     - 8._dp*log(y31)*log(y23)
     &     + 6._dp*log(y23)**2
     &     + 4._dp*Iijm(1)
     &     + 4._dp*Iijm(2)
     &     + 4._dp*Iijm(3)
     &     + 4._dp*Iijm(4)
     &     - 4._dp*Iijm(5)
     &     - 4._dp*Iijm(6)
     &     )
      soft2ab(1) = soft2ab(1) + CF*CA * (
     &     - 72._dp*zeta2
     &     + 4._dp*log(y12)**2
     &     + 16._dp*log(y12)*log(y31)
     &     - 16._dp*log(y31)**2
     &     + 16._dp*log(y31)*log(y23)
     &     + 4._dp*log(y23)**2
     &     + 8._dp*Iijm(1)
     &     + 8._dp*Iijm(2)
     &     + 8._dp*Iijm(3)
     &     + 8._dp*Iijm(4)
     &     )
      soft2ab(1) = soft2ab(1) + CF**2 * (
     &     - 72._dp*zeta2
     &     + 24._dp*log(y31)**2
     &     + 16._dp*Iijm(5)
     &     + 16._dp*Iijm(6)
     &     )

      soft2ab(0) =  + CA**2 * (
     &     + 16._dp*zeta3
     &     + 9._dp*log(y12)*zeta2
     &     - log(y12)**3
     &     + log(y12)**2*log(y31)
     &     - log(y12)**2*log(y23)
     &     + log(y12)*log(y31)**2
     &     - log(y12)*log(y23)**2
     &     - 2._dp*log(y12)*Iijm(1)
     &     - 2._dp*log(y12)*Iijm(2)
     &     - 2._dp*log(y12)*Iijm(3)
     &     - 2._dp*log(y12)*Iijm(4)
     &     + 2._dp*log(y12)*Iijm(5)
     &     + 2._dp*log(y12)*Iijm(6)
     &     - 9._dp*log(y31)*zeta2
     &     )
      soft2ab(0) = soft2ab(0) + CA**2 * (
     &     - log(y31)**3
     &     + log(y31)**2*log(y23)
     &     + log(y31)*log(y23)**2
     &     + 2._dp*log(y31)*Iijm(1)
     &     + 2._dp*log(y31)*Iijm(2)
     &     + 2._dp*log(y31)*Iijm(3)
     &     + 2._dp*log(y31)*Iijm(4)
     &     - 2._dp*log(y31)*Iijm(5)
     &     - 2._dp*log(y31)*Iijm(6)
     &     + 9._dp*log(y23)*zeta2
     &     - log(y23)**3
     &     - 2._dp*log(y23)*Iijm(1)
     &     - 2._dp*log(y23)*Iijm(2)
     &     - 2._dp*log(y23)*Iijm(3)
     &     )
      soft2ab(0) = soft2ab(0) + CA**2 * (
     &     - 2._dp*log(y23)*Iijm(4)
     &     + 2._dp*log(y23)*Iijm(5)
     &     + 2._dp*log(y23)*Iijm(6)
     &     )
      soft2ab(0) = soft2ab(0) + CF*CA * (
     &     + 64._dp*zeta3
     &     + 18._dp*log(y12)*zeta2
     &     - 2._dp*log(y12)**2*log(y31)
     &     - 2._dp*log(y12)*log(y31)**2
     &     - 4._dp*log(y12)*Iijm(5)
     &     - 4._dp*log(y12)*Iijm(6)
     &     + 4._dp*log(y31)**3
     &     - 2._dp*log(y31)**2*log(y23)
     &     - 2._dp*log(y31)*log(y23)**2
     &     - 4._dp*log(y31)*Iijm(1)
     &     - 4._dp*log(y31)*Iijm(2)
     &     - 4._dp*log(y31)*Iijm(3)
     &     - 4._dp*log(y31)*Iijm(4)
     &     + 8._dp*log(y31)*Iijm(5)
     &     )
      soft2ab(0) = soft2ab(0) + CF*CA * (
     &     + 8._dp*log(y31)*Iijm(6)
     &     + 18._dp*log(y23)*zeta2
     &     - 4._dp*log(y23)*Iijm(5)
     &     - 4._dp*log(y23)*Iijm(6)
     &     )
      soft2ab(0) = soft2ab(0) + CF**2 * (
     &     + 64._dp*zeta3
     &     + 36._dp*log(y31)*zeta2
     &     - 4._dp*log(y31)**3
     &     - 8._dp*log(y31)*Iijm(5)
     &     - 8._dp*log(y31)*Iijm(6)
     &     )

      soft2ab(-1) =  + CA**2 * (
     &     - 27._dp/40._dp*zeta2**2
     &     - 8._dp*log(y12)*zeta3
     &     - 9._dp/4._dp*log(y12)**2*zeta2
     &     + 1._dp/8._dp*log(y12)**4
     &     - 1._dp/4._dp*log(y12)**2*log(y31)**2
     &     + 1._dp/4._dp*log(y12)**2*log(y23)**2
     &     + 1._dp/2._dp*log(y12)**2*Iijm(1)
     &     + 1._dp/2._dp*log(y12)**2*Iijm(2)
     &     + 1._dp/2._dp*log(y12)**2*Iijm(3)
     &     + 1._dp/2._dp*log(y12)**2*Iijm(4)
     &     - 1._dp/2._dp*log(y12)**2*Iijm(5)
     &     - 1._dp/2._dp*log(y12)**2*Iijm(6)
     &     + 4._dp*log(y12)*log(y31)*zeta2
     &     - 4._dp*log(y12)*log(y23)*zeta2
     &     )
      soft2ab(-1) = soft2ab(-1) + CA**2 * (
     &     + 8._dp*log(y31)*zeta3
     &     - 7._dp/4._dp*log(y31)**2*zeta2
     &     + 1._dp/8._dp*log(y31)**4
     &     - 1._dp/4._dp*log(y31)**2*log(y23)**2
     &     - 1._dp/2._dp*log(y31)**2*Iijm(1)
     &     - 1._dp/2._dp*log(y31)**2*Iijm(2)
     &     - 1._dp/2._dp*log(y31)**2*Iijm(3)
     &     - 1._dp/2._dp*log(y31)**2*Iijm(4)
     &     + 1._dp/2._dp*log(y31)**2*Iijm(5)
     &     + 1._dp/2._dp*log(y31)**2*Iijm(6)
     &     + 4._dp*log(y31)*log(y23)*zeta2
     &     - 8._dp*log(y23)*zeta3
     &     - 9._dp/4._dp*log(y23)**2*zeta2
     &     + 1._dp/8._dp*log(y23)**4
     &     )
      soft2ab(-1) = soft2ab(-1) + CA**2 * (
     &     + 1._dp/2._dp*log(y23)**2*Iijm(1)
     &     + 1._dp/2._dp*log(y23)**2*Iijm(2)
     &     + 1._dp/2._dp*log(y23)**2*Iijm(3)
     &     + 1._dp/2._dp*log(y23)**2*Iijm(4)
     &     - 1._dp/2._dp*log(y23)**2*Iijm(5)
     &     - 1._dp/2._dp*log(y23)**2*Iijm(6)
     &     - 1._dp/2._dp*Iijm(1)*zeta2
     &     + 1._dp/2._dp*Iijm(1)**2
     &     + Iijm(1)*Iijm(2)
     &     + Iijm(1)*Iijm(3)
     &     + Iijm(1)*Iijm(4)
     &     - Iijm(1)*Iijm(5)
     &     - Iijm(1)*Iijm(6)
     &     - 1._dp/2._dp*Iijm(2)*zeta2
     &     )
      soft2ab(-1) = soft2ab(-1) + CA**2 * (
     &     + 1._dp/2._dp*Iijm(2)**2
     &     + Iijm(2)*Iijm(3)
     &     + Iijm(2)*Iijm(4)
     &     - Iijm(2)*Iijm(5)
     &     - Iijm(2)*Iijm(6)
     &     - 1._dp/2._dp*Iijm(3)*zeta2
     &     + 1._dp/2._dp*Iijm(3)**2
     &     + Iijm(3)*Iijm(4)
     &     - Iijm(3)*Iijm(5)
     &     - Iijm(3)*Iijm(6)
     &     - 1._dp/2._dp*Iijm(4)*zeta2
     &     + 1._dp/2._dp*Iijm(4)**2
     &     - Iijm(4)*Iijm(5)
     &     - Iijm(4)*Iijm(6)
     &     )
      soft2ab(-1) = soft2ab(-1) + CA**2 * (
     &     + 1._dp/2._dp*Iijm(5)*zeta2
     &     + 1._dp/2._dp*Iijm(5)**2
     &     + Iijm(5)*Iijm(6)
     &     + 1._dp/2._dp*Iijm(6)*zeta2
     &     + 1._dp/2._dp*Iijm(6)**2
     &     )
      soft2ab(-1) = soft2ab(-1) + CF*CA * (
     &     - 27._dp/10._dp*zeta2**2
     &     - 16._dp*log(y12)*zeta3
     &     - 1._dp/2._dp*log(y12)**2*zeta2
     &     + 1._dp/2._dp*log(y12)**2*log(y31)**2
     &     + log(y12)**2*Iijm(5)
     &     + log(y12)**2*Iijm(6)
     &     - 8._dp*log(y12)*log(y31)*zeta2
     &     + 8._dp*log(y31)**2*zeta2
     &     - 1._dp/2._dp*log(y31)**4
     &     + 1._dp/2._dp*log(y31)**2*log(y23)**2
     &     + log(y31)**2*Iijm(1)
     &     + log(y31)**2*Iijm(2)
     &     + log(y31)**2*Iijm(3)
     &     + log(y31)**2*Iijm(4)
     &     )
      soft2ab(-1) = soft2ab(-1) + CF*CA * (
     &     - 2._dp*log(y31)**2*Iijm(5)
     &     - 2._dp*log(y31)**2*Iijm(6)
     &     - 8._dp*log(y31)*log(y23)*zeta2
     &     - 16._dp*log(y23)*zeta3
     &     - 1._dp/2._dp*log(y23)**2*zeta2
     &     + log(y23)**2*Iijm(5)
     &     + log(y23)**2*Iijm(6)
     &     - Iijm(1)*zeta2
     &     + 2._dp*Iijm(1)*Iijm(5)
     &     + 2._dp*Iijm(1)*Iijm(6)
     &     - Iijm(2)*zeta2
     &     + 2._dp*Iijm(2)*Iijm(5)
     &     + 2._dp*Iijm(2)*Iijm(6)
     &     - Iijm(3)*zeta2
     &     )
      soft2ab(-1) = soft2ab(-1) + CF*CA * (
     &     + 2._dp*Iijm(3)*Iijm(5)
     &     + 2._dp*Iijm(3)*Iijm(6)
     &     - Iijm(4)*zeta2
     &     + 2._dp*Iijm(4)*Iijm(5)
     &     + 2._dp*Iijm(4)*Iijm(6)
     &     - 2._dp*Iijm(5)**2
     &     - 4._dp*Iijm(5)*Iijm(6)
     &     - 2._dp*Iijm(6)**2
     &     )
      soft2ab(-1) = soft2ab(-1) + CF**2 * (
     &     - 27._dp/10._dp*zeta2**2
     &     - 32._dp*log(y31)*zeta3
     &     - 9._dp*log(y31)**2*zeta2
     &     + 1._dp/2._dp*log(y31)**4
     &     + 2._dp*log(y31)**2*Iijm(5)
     &     + 2._dp*log(y31)**2*Iijm(6)
     &     - 2._dp*Iijm(5)*zeta2
     &     + 2._dp*Iijm(5)**2
     &     + 4._dp*Iijm(5)*Iijm(6)
     &     - 2._dp*Iijm(6)*zeta2
     &     + 2._dp*Iijm(6)**2
     &     )

!---- non-abelian contribution
      soft2nab(3)=zip
      soft2nab(2)=2*be0*(CA+2*CF)
      soft2nab(1)=+CA*(-half*Ga1
     & -two*log(y12)*be0+two*log(y31)*be0-two*log(y23)*be0)
     & +CF*(-Ga1-four*log(y31)*be0)
      soft2nab(0)=CA*(
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
      soft2nab(0) = soft2nab(0) + CF * (
     &     - half*gams1
     &     - zeta2*be0
     &     + half*log(y31)*Ga1
     &     + log(y31)**2*be0
     &     + two*Iijm(5)*be0
     &     + two*Iijm(6)*be0
     &     )

c---- BLP fit 1504.02540
!      Lss13=log(y31)
!      Lss23=log(y23)
!      soft2nab(-1)=-200.096_dp-223.007_dp*Lss13-210.796_dp*Lss23
!     & -112.644_dp*Lss13**2-118.705_dp*Lss23**2
!     & -30.8755_dp*Lss13**3-32.7046_dp*Lss23**3
!     & -5.57982_dp*Lss13**4-5.43977_dp*Lss23**4
!     & -0.666669_dp*Lss13**5-0.63292_dp*Lss23**5
!     & -0.0533767_dp*Lss13**6-0.0484011_dp*Lss23**6
!     & -0.0028547_dp*Lss13**7-0.00240717_dp*Lss23**7
!     & -0.0000952163_dp*Lss13**8-0.0000728699_dp*Lss23**8
!     & -1.45696e-6_dp*Lss13**9-9.83475e-7_dp*Lss23**9

c---- CEMW fit 1711.09984
      c00=39.101_dp
      c10=13.726_dp
      c20=-2.737_dp
      c30=0.016_dp
      c01=25.591_dp
      c02=-8.749_dp
      c03=-2.126_dp

      soft2nab(-1)=c00
     & +c10*log(y31)+c20*log(y31)**2+c30*log(y31)**3
     & +c01*log(one-y31)+c02*log(one-y31)**2+c03*log(one-y31)**3

!---- contribution to delta(tau)
      soft2(-1) = soft2ab(-1)+soft2nab(-1)
!---- contribution to L0(tau)
      soft2(0) = soft2ab(0)+soft2nab(0)
!---- contribution to L1(tau)
      soft2(1) = soft2ab(1)+soft2nab(1)
!---- contribution to L2(tau)
      soft2(2) = soft2ab(2)+soft2nab(2)
!---- contribution to L3(tau)
      soft2(3) = soft2ab(3)+soft2nab(3)

!      soft2(-1) = soft2nab(-1) ! debug
!      soft2(0) = soft2nab(0) ! debug
!      soft2(1) = soft2nab(1) ! debug
!      soft2(2) = soft2nab(2) ! debug
!      soft2(3) = soft2nab(3) ! debug

      return
      end

