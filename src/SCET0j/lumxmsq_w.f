      subroutine lumxmsq_w(p,xx,z1,z2,QB,order,xmsq)
      implicit none
      include 'types.f'
c----Matrix element for W production
C----averaged over initial colours and spins
C For nwz=+1
c     u(-p1)+dbar(-p2)-->W^+(n(p3)+e^+(p4))
C For nwz=-1
c     d(-p1)+ubar(-p2)-->W^-(e^-(p3)+nbar(p4))
c---
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'scale.f'
      include 'facscale.f'
      include 'ckm.f'
      include 'scet_const.f'
      include 'taucut.f'
      integer:: j,k,ih1,ih2,m,n,order
      real(dp),intent(out):: xmsq
      real(dp):: p(mxpart,4),s,fac,qqb,qbq,tauc,getdynamictau,
     & xx(2),soft1(-1:1),soft2(-1:3),hard(2),
     & beama0(-5:5),beamb0(-5:5),
     & beama1(-5:5,-1:1),beamb1(-5:5,-1:1),
     & beama2(-5:5,-1:3),beamb2(-5:5,-1:3),
     & z1,z2,QB(2),lum0,lum1(-1:1),lum2(-1:3),bit,lobit,
     & msq(-nf:nf,-nf:nf),assemble,msqpow(-nf:nf,-nf:nf)
      common/density/ih1,ih2
      s(j,k)=two*(p(j,4)*p(k,4)-p(j,1)*p(k,1)
     &           -p(j,2)*p(k,2)-p(j,3)*p(k,3))

      fac=gw**4*xn
      fac=aveqq*fac/((s(3,4)-wmass**2)**2+(wmass*wwidth)**2)
c---case dbar-u or ubar-d
      qqb=fac*s(1,4)**2
      qbq=fac*s(2,4)**2

      call softqqbis(order,soft1,soft2)
      call hardqq(s(1,2),musq,hard)

      if (order >= 0) then
      call fdist_select(ih1,1,xx(1),facscale,beama0)
      call fdist_select(ih2,2,xx(2),facscale,beamb0)
      endif
      if (order >= 1) then
      call xbeam1bis(ih1,1,z1,xx(1),QB(1),beama1)
      call xbeam1bis(ih2,2,z2,xx(2),QB(2),beamb1)
      endif
      if (order >= 2) then
      call xbeam2bis(ih1,1,z1,xx(1),QB(1),beama2)
      call xbeam2bis(ih2,2,z2,xx(2),QB(2),beamb2)
      endif

! compute power corrections if necessary
      call tau0_powcorr_qa(order,xx(1),xx(2),sqrt(s(3,4)),beama0,beamb0,msqpow)

      if (dynamictau) then
        tauc=getdynamictau(p)
      else
        tauc=taucut
      endif

      xmsq=zip
      do j=-nf,nf
      do k=-nf,nf
      if (j*k >= 0) cycle ! skip gluons, qq, aa
      
      bit=assemble(order,tauc,
     & beama0(j),beamb0(k),beama1(j,:),beamb1(k,:),
     & beama2(j,:),beamb2(k,:),soft1,soft2,hard)

      if ((j > 0) .and. (k < 0)) then
        lobit=Vsq(j,k)*qqb
      elseif ((j < 0) .and. (k > 0)) then
        lobit=Vsq(j,k)*qbq
      else
        lobit=zip
      endif
      
! include power corrections
!      bit=bit*zip+(msqpow(j,k)+msqpow(j,0)+msqpow(0,k))
      
      xmsq=xmsq+bit*lobit

      enddo
      enddo

      return
      end
