      subroutine qqb_dirgam_ew(p,msq)
! Implementation of EW corrections to direct photon production according to
! J. Kuhn, A. Kulesza, S. Pozzorini, M. Schulze, hep-ph/0508253
      implicit none
      include 'types.f'
C-----Matrix element for f(-p1)+f(-p2)->gamma(p3)+g(p4)
      include 'constants.f'
      include 'epinv.f'
      include 'masses.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'ewcharge.f'
      include 'sprods_com.f'
      integer:: j,k
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),fac,H0,H1A,H1N,mv2,
     & shat,that,uhat,g0A,g1A,g2A,g0N,g1N,g2N,aew,aewon2pi,sw,cw,
     & qa(5),qg(5),gq(5),IZIZ(5),IWIW(5),T3(5),dCN
! Statement functions
! Lowest order kinematic structure
      H0(that,uhat)=that/uhat+uhat/that
! Equation (16)
      H1A(shat,that,uhat,mv2)=
     & +g0A(shat,that,uhat,mv2)*H0(that,uhat)
     & +g1A(shat,that,uhat)*(that/uhat-uhat/that)
     & +g2A(shat,that,uhat)
      H1N(shat,that,uhat,mv2)=
     & +g0N(shat,that,uhat,mv2)*H0(that,uhat)
     & +g1N(shat,that,uhat)*(that/uhat-uhat/that)
     & +g2N(shat,that,uhat)
     
      aew=esq/fourpi
      aewon2pi=aew/twopi!*zip ! zero corrections

! Additional factor of 4 here since Eq.(15) already includes 1/4 for
! average over initial-state polarizations
      fac=four*(eight*pisq*aew*as*V)

      sw=sqrt(xw)
      cw=sqrt(one-xw)
      T3(:)=tau(1:5)/two
! Coupling structures according to Eq. (10), summed over chiralities
      IWIW(:)=half/xw
      IZIZ(:)=(cw/sw*T3(:)-sw/cw*(Q(1:5)-T3(:)))**2 ! L
     &       +(           -sw/cw* Q(1:5))**2        ! R
      

! Counter-term, Eq.(14) (only keep pole in DeltaUVbar as check, should cancel)
      dCN=-half/sw/cw*(aewon2pi/two)*four*cw/sw*(epinv-log(wmass**2/zmass**2))

! Kinematic quantities
      call dotem(3,p,s)
      shat=s(1,2)
      that=s(2,3) ! =s(1,4)
      uhat=s(1,3) ! =s(2,3)
      
!      fac=4._dp*V*gsq*esq
!      qa=fac*aveqq*(s(1,3)/s(2,3)+s(2,3)/s(1,3))
!      aq=qa
!      qg=-fac*aveqg*(s(1,3)/s(1,2)+s(1,2)/s(1,3))
!      ag=qg
!      gq=-fac*aveqg*(s(1,2)/s(2,3)+s(2,3)/s(1,2))
!      ga=gq
      qa(:)=aveqq*fac*(
     &      Q(1:5)**2*(two*H0(that,uhat)+aewon2pi*(
     &        IZIZ(:)*H1A(shat,that,uhat,zmass**2)
     &       +IWIW(:)*H1A(shat,that,uhat,wmass**2)))
     &       +T3(:)*Q(1:5)*(two*H0(that,uhat)*dCN
     &                   +aewon2pi/xw*H1N(shat,that,uhat,wmass**2)))
     
      qg(:)=-aveqg*fac*(
     &      Q(1:5)**2*(two*H0(shat,uhat)+aewon2pi*(
     &        IZIZ(:)*H1A(that,shat,uhat,zmass**2)
     &       +IWIW(:)*H1A(that,shat,uhat,wmass**2)))
     &       +T3(:)*Q(1:5)*(two*H0(shat,uhat)*dCN
     &                   +aewon2pi/xw*H1N(that,shat,uhat,wmass**2)))
     
      gq(:)=-aveqg*fac*(
     &      Q(1:5)**2*(two*H0(that,shat)+aewon2pi*(
     &        IZIZ(:)*H1A(uhat,that,shat,zmass**2)
     &       +IWIW(:)*H1A(uhat,that,shat,wmass**2)))
     &       +T3(:)*Q(1:5)*(two*H0(that,shat)*dCN
     &                   +aewon2pi/xw*H1N(uhat,that,shat,wmass**2)))
     
      do j=-nf,nf
      do k=-nf,nf
c--set msq=0 to initalize
      msq(j,k)=0._dp
C--qa      
      if     ((j > 0) .and. (k == -j)) then
        msq(j,k)=qa(j)
C--aq      
      elseif ((j < 0) .and. (k == -j)) then
        msq(j,k)=qa(k)
C--qg
      elseif ((j > 0) .and. (k == 0)) then
        msq(j,k)=qg(j)
C--ag      
      elseif ((j < 0) .and. (k == 0)) then
        msq(j,k)=qg(-j)
C--gq
      elseif ((j == 0) .and. (k > 0)) then
        msq(j,k)=gq(k)
C--ga      
      elseif ((j == 0) .and. (k < 0)) then
        msq(j,k)=gq(-k)
      endif

      enddo
      enddo


      return
      end


      
