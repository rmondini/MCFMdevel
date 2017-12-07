      subroutine soft_nab_ggg(order,y12,y23,y31,
     & preIijm,j1,j2,j3,j4,j5,j6,soft2)
      implicit none
!     non-abelian piece
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
     & Iijm(6),Lss
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
      soft2(2) =6._dp*CA*be0
      soft2(1) = soft2(1) 
     &     - 3._dp/2._dp*CA*Ga1
     &     - 2._dp*log(y12)*CA*be0
     &     - 2._dp*log(y31)*CA*be0
     &     - 2._dp*log(y23)*CA*be0

      soft2(0) =
     &     - 3._dp/4._dp*CA*gams1
     &     - 3._dp/2._dp*CA*zeta2*be0
     &     + quarter*log(y12)*CA*Ga1
     &     + half*log(y12)**2*CA*be0
     &     + quarter*log(y31)*CA*Ga1
     &     + half*log(y31)**2*CA*be0
     &     + quarter*log(y23)*CA*Ga1
     &     + half*log(y23)**2*CA*be0
     &     + Iijm(1)*CA*be0
     &     + Iijm(2)*CA*be0
     &     + Iijm(3)*CA*be0
     &     + Iijm(4)*CA*be0
     &     + Iijm(5)*CA*be0
     &     + Iijm(6)*CA*be0

      Lss=log(y31)+log(y23)
      soft2(-1)=74.1772_dp+30.2529_dp*Lss-18.7425_dp*Lss**2
     & -5.61366_dp*Lss**3-0.584101_dp*Lss**4
     & -0.0428711_dp*Lss**5-0.000301809_dp*Lss**6
     & +0.00010388_dp*Lss**7+(2.55309e-6_dp)*Lss**8

      return
      end
