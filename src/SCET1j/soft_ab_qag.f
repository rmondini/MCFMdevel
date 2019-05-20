      subroutine soft_ab_qag(order,y12,y23,y31,
     & preIijm,j1,j2,j3,j4,j5,j6,soft1,soft2)
      implicit none
!---   Abelian piece q(1) qbar(2) g(3)
!--- returns coefficients of [log(tau)^j/tau_+],and for j=-1, delta(tau)
!--- in units of as/2/pi
!     y12=n1.n2/2 so that 0<y12<1
!     y23=n2.n3/2 so that 0<y23<1
!     y31=n3.n1/2 so that 0<y31<1
!     order indicates the order in [as/2/pi]
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'scet_const.f'
      real(dp)::y12,y23,y31,soft1(-1:1),soft2(-1:3),preIijm(5),
     & Iijm(6),Lss,Lsd,Isum12,Isum3456
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

!      if ((y12 < zip) .or. (y12 > one)) then
!        write(6,*) 'error: y12 outside range, y12 = ',y12
!      endif
!      if ((y23 < zip) .or. (y23 > one)) then
!        write(6,*) 'error: y23 outside range, y23 = ',y23
!      endif
!      if ((y31 < zip) .or. (y31 > one)) then
!        write(6,*) 'error: y31 outside range, y31 = ',y31
!      endif
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

      Lss=log(y31)+log(y23)
      Lsd=log(y31)-log(y23)
      Isum12=Iijm(1)+Iijm(2)
      Isum3456=Iijm(3)+Iijm(4)+Iijm(5)+Iijm(6)

      soft1(1)=-four*CA-eight*CF
      soft1(0)=two*CA*(Lss-log(y12))+four*log(y12)*CF

      soft1(-1)=CA*(-0.25_dp*(Lsd**2+Lss**2)+Isum12-Isum3456
     & +half*zeta2+half*log(y12)**2)
     & +CF*(-two*Isum12+zeta2-log(y12)**2)

      if (order < 2) return

      soft2(3)=8._dp*(CA+2*CF)**2

      soft2(2)=-12._dp*CA**2*Lss-24._dp*CF*CA*Lss

      soft2(1) = + CA**2*Lsd**2
     &  + 5._dp*CA**2*Lss**2
     &  - 4._dp*CA**2*Isum12
     &  + 4._dp*CA**2*Isum3456
     &  - 18._dp*CA**2*zeta2
     &  + 2._dp*CF*CA*Lsd**2
     &  + 2._dp*CF*CA*Lss**2
     &  + 8._dp*CF*CA*Isum3456
     &  - 72._dp*CF*CA*zeta2
     &  + 16._dp*CF**2*Isum12
     &  - 72._dp*CF**2*zeta2

      soft2(0) = - 1._dp/2._dp*CA**2*Lss*Lsd**2
     &  - 1._dp/2._dp*CA**2*Lss**3
     &  + 2._dp*CA**2*Isum12*Lss
     &  - 2._dp*CA**2*Isum3456*Lss
     &  + 16._dp*CA**2*zeta3
     &  + 9._dp*CA**2*zeta2*Lss
     &  - 4._dp*CF*CA*Isum12*Lss
     &  + 64._dp*CF*CA*zeta3
     &  + 18._dp*CF*CA*zeta2*Lss
     &  + 64._dp*CF**2*zeta3

      soft2(-1) = + 1._dp/32._dp*CA**2*Lsd**4
     &  + 1._dp/16._dp*CA**2*Lss**2*Lsd**2
     &  + 1._dp/32._dp*CA**2*Lss**4
     &  - 1._dp/4._dp*CA**2*Isum12*Lsd**2
     &  - 1._dp/4._dp*CA**2*Isum12*Lss**2
     &  + 1._dp/2._dp*CA**2*Isum12**2
     &  + 1._dp/4._dp*CA**2*Isum3456*Lsd**2
     &  + 1._dp/4._dp*CA**2*Isum3456*Lss**2
     &  - CA**2*Isum3456*Isum12
     &  + 1._dp/2._dp*CA**2*Isum3456**2
     &  - 8._dp*CA**2*zeta3*Lss
     &  - 1._dp/8._dp*CA**2*zeta2*Lsd**2
     &  - 17._dp/8._dp*CA**2*zeta2*Lss**2
     &  + 1._dp/2._dp*CA**2*zeta2*Isum12
     &  - 1._dp/2._dp*CA**2*zeta2*Isum3456
      soft2(-1) = soft2(-1) - 27._dp/40._dp*CA**2*zeta2**2
     &  + 1._dp/2._dp*CF*CA*Isum12*Lsd**2
     &  + 1._dp/2._dp*CF*CA*Isum12*Lss**2
     &  - 2._dp*CF*CA*Isum12**2
     &  + 2._dp*CF*CA*Isum3456*Isum12
     &  - 16._dp*CF*CA*zeta3*Lss
     &  - 1._dp/4._dp*CF*CA*zeta2*Lsd**2
     &  - 1._dp/4._dp*CF*CA*zeta2*Lss**2
     &  - CF*CA*zeta2*Isum3456
     &  - 27._dp/10._dp*CF*CA*zeta2**2
     &  + 2._dp*CF**2*Isum12**2
     &  - 2._dp*CF**2*zeta2*Isum12
     &  - 27._dp/10._dp*CF**2*zeta2**2

      return
      end
