      subroutine soft_ab_qgq(order,y12,y23,y31,
     & preIijm,j1,j2,j3,j4,j5,j6,soft1,soft2)
      implicit none
!     Abelian piece q(1) g(2) q(3)
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
      real(dp)::y12,y23,y31,soft1(-1:1),soft2(-1:3),preIijm(5),
     & Iijm(6)
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

      soft1(:)=zip
      soft2(:)=zip
      if (order < 1) return

!      y(:,:)=zip
!      y(1,2)=y12
!      y(2,1)=y12
!      y(2,3)=y23
!      y(3,2)=y23
!      y(1,3)=y31
!      y(3,1)=y31

c      if ((y12 < zip) .or. (y12 > one)) then
c        write(6,*) 'error: y12 outside range, y12 = ',y12
c      endif
c      if ((y23 < zip) .or. (y23 > one)) then
c        write(6,*) 'error: y23 outside range, y23 = ',y23
c      endif
c      if ((y31 < zip) .or. (y31 > one)) then
c        write(6,*) 'error: y31 outside range, y31 = ',y31
c      endif

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

      soft1(-1)=CA*(half*zeta2
     & -half*log(y12)**2+half*log(y31)**2-half*log(y23)**2
     & -Iijm(1)-Iijm(2)-Iijm(3)-Iijm(4)+Iijm(5)+Iijm(6))
     & +CF*(zeta2-log(y31)**2-2*Iijm(5)-2*Iijm(6))
      soft1(0)=2*CA*(log(y12)+log(y23)-log(y31))+4*log(y31)*CF
      soft1(1)=-4*CA-8*CF
      
      if (order < 2) return
      
!-------------------
      soft2(3)=8*(CA+2*CF)**2
      soft2(2)=-12*(CA+2*CF)*(CA*(log(y12)+log(y23)-log(y31))
     & +2*log(y31)*CF)

      soft2(1) =  + CA**2 * (
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
      soft2(1) = soft2(1) + CF*CA * (
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
      soft2(1) = soft2(1) + CF**2 * (
     &     - 72._dp*zeta2
     &     + 24._dp*log(y31)**2
     &     + 16._dp*Iijm(5)
     &     + 16._dp*Iijm(6)
     &     )

      soft2(0) =  + CA**2 * (
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
      soft2(0) = soft2(0) + CA**2 * (
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
      soft2(0) = soft2(0) + CA**2 * (
     &     - 2._dp*log(y23)*Iijm(4)
     &     + 2._dp*log(y23)*Iijm(5)
     &     + 2._dp*log(y23)*Iijm(6)
     &     )
      soft2(0) = soft2(0) + CF*CA * (
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
      soft2(0) = soft2(0) + CF*CA * (
     &     + 8._dp*log(y31)*Iijm(6)
     &     + 18._dp*log(y23)*zeta2
     &     - 4._dp*log(y23)*Iijm(5)
     &     - 4._dp*log(y23)*Iijm(6)
     &     )
      soft2(0) = soft2(0) + CF**2 * (
     &     + 64._dp*zeta3
     &     + 36._dp*log(y31)*zeta2
     &     - 4._dp*log(y31)**3
     &     - 8._dp*log(y31)*Iijm(5)
     &     - 8._dp*log(y31)*Iijm(6)
     &     )

      soft2(-1) =  + CA**2 * (
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
      soft2(-1) = soft2(-1) + CA**2 * (
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
      soft2(-1) = soft2(-1) + CA**2 * (
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
      soft2(-1) = soft2(-1) + CA**2 * (
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
      soft2(-1) = soft2(-1) + CA**2 * (
     &     + 1._dp/2._dp*Iijm(5)*zeta2
     &     + 1._dp/2._dp*Iijm(5)**2
     &     + Iijm(5)*Iijm(6)
     &     + 1._dp/2._dp*Iijm(6)*zeta2
     &     + 1._dp/2._dp*Iijm(6)**2
     &     )
      soft2(-1) = soft2(-1) + CF*CA * (
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
      soft2(-1) = soft2(-1) + CF*CA * (
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
      soft2(-1) = soft2(-1) + CF*CA * (
     &     + 2._dp*Iijm(3)*Iijm(5)
     &     + 2._dp*Iijm(3)*Iijm(6)
     &     - Iijm(4)*zeta2
     &     + 2._dp*Iijm(4)*Iijm(5)
     &     + 2._dp*Iijm(4)*Iijm(6)
     &     - 2._dp*Iijm(5)**2
     &     - 4._dp*Iijm(5)*Iijm(6)
     &     - 2._dp*Iijm(6)**2
     &     )
      soft2(-1) = soft2(-1) + CF**2 * (
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
      return
      end
