      subroutine soft_nab_qag(order,y12,y23,y31,
     & preIijm,j1,j2,j3,j4,j5,j6,soft2)
      implicit none
!--- Non-abelian piece q(1) qbar(2) g(3)
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
      real(dp)::y12,y23,y31,soft2(-1:3),preIijm(5),
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

      soft2(:)=zip
      if (order < 2) return

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


      soft2(3)= 0
      soft2(2)=2._dp*CA*be0 + 4._dp*CF*be0
      soft2(1)= - half*CA*Ga1 - 2._dp*CA*be0*Lss - CF*Ga1
      soft2(0)= - quarter*CA*gams1 + quarter*CA*Ga1*Lss + CA*be0*
     & Isum3456 - CA*be0*Isum12 + quarter*CA*be0*Lsd**2 + quarter*
     & CA*be0*Lss**2 - half*CA*zeta2*be0 - half*CF*gams1 
     & + 2._dp*CF*be0*Isum12 - CF*zeta2*be0

      soft2(-1)=60.1426_dp+41.0237_dp*Lss-2.57381_dp*Lss**2
     & -0.692668_dp*Lss**3+0.226885_dp*Lss**4
     & +0.0209989_dp*Lss**5+0.00090579_dp*Lss**6
     & +0.0000132008_dp*Lss**7
     
!      soft2(-1)=soft2(-1)/four ! DEBUG

      return
      end
