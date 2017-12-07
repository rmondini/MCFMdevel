      subroutine lumxmsq_h1jet(p,xx,z1,z2,QB,order,xmsq)
      implicit none
      include 'types.f'
c---- Matrix element for H production
c---- in the heavy quark (mt=Infinity) limit.
C---- averaged over initial colours and spins
c     g(-p1) + g(-p2) --> H(-->b(p3)+b~(p4)) + g(p5)
c---
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'scale.f'
      include 'facscale.f'
      include 'hdecaymode.f'
      include 'scet_const.f'
      include 'taucut.f'
      integer:: j,k,ih1,ih2,m,n,order
      real(dp),intent(out):: xmsq
      real(dp):: p(mxpart,4),fac,Asq,xx(2),hard(2),
     & soft1gg(-1:1),soft2gg(-1:3),soft2gg_nab(-1:3),
     & soft1qg(-1:1),soft2qg(-1:3),soft2qg_nab(-1:3),
     & soft1gq(-1:1),soft2gq(-1:3),soft2gq_nab(-1:3),
     & soft1qa(-1:1),soft2qa(-1:3),soft2qa_nab(-1:3),
     & beama0(-5:5),beamb0(-5:5),tauc,getdynamictau,
     & beama1(-5:5,-1:1),beamb1(-5:5,-1:1),
     & beama2(-5:5,-1:3),beamb2(-5:5,-1:3),
     & jet1q(-1:1),jet2q(-1:3),jet1g(-1:1),jet2g(-1:3),
     & z1,z2,QB(2),lum0,lum1(-1:1),lum2(-1:3),bit,
     & msqgg,msqqa,msqqg,msqgq,assemblejet,y12,y15,y25,Iijm(5),
     & hdecay,msqgamgam,ss,tt,uu,mhsq,s(mxpart,mxpart),s34
      real(dp):: gg,qqb,qg,gqb
      real(dp):: Hgg(2),Hqqb(2),Hqg(2),Hgqb(2)
      common/density/ih1,ih2
      integer, parameter:: iglue=5

      call dotem(iglue,p,s)

c--- note: important to use this statement (rather than dot product) in
c--- order to retain the possibility of massive decay products
      s34=(p(3,4)+p(4,4))**2
     &   -(p(3,1)+p(4,1))**2-(p(3,2)+p(4,2))**2-(p(3,3)+p(4,3))**2

C   Deal with Higgs decay
      if (hdecaymode == 'tlta') then
          call htautaudecay(p,3,4,hdecay)
      elseif (hdecaymode == 'bqba') then
          call hbbdecay(p,3,4,hdecay)
      elseif (hdecaymode == 'gaga') then
          hdecay=msqgamgam(hmass)
      else
      write(6,*) 'Unimplemented process in gg_hg'
      stop
      endif
      hdecay=hdecay/((s34-hmass**2)**2+(hmass*hwidth)**2)

      call dotem(iglue,p,s)
c--- compute matrix elements
      call Hampgggsq(order,5,1,2,3,4,gg,Hgg(1),Hgg(2),qg,Hqg(1),Hqg(2))
      call Hampgggsq(order,2,5,1,3,4,gg,Hgg(1),Hgg(2),gqb,Hgqb(1),Hgqb(2))
      call Hampgggsq(order,2,1,5,3,4,gg,Hgg(1),Hgg(2),qqb,Hqqb(1),Hqqb(2))
      
      Asq=(as/(3._dp*pi))**2/vevsq
      fac=Asq*gsq*hdecay

      msqgg=avegg*fac*V*xn*gg
      msqqa=+aveqq*fac*V/two*qqb
      msqgq=-aveqg*fac*V/two*gqb
      msqqg=-aveqg*fac*V/two*qg

c--- beam functions
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

c--- jet functions
      call jetq(order,two*p(iglue,4),jet1q,jet2q)
      call jetg(order,two*p(iglue,4),jet1g,jet2g)
      
c---- soft functions
      y12=s(1,2)/p(1,4)/p(2,4)/four
      y15=s(1,iglue)/p(1,4)/p(iglue,4)/four
      y25=s(2,iglue)/p(2,4)/p(iglue,4)/four
!      call soft_ab_ggg_old(order,y12,y15,y25,soft1gg,soft2gg)
!      call soft_ab_qgq_old(order,y12,y25,y15,soft1qg,soft2qg)
!      call soft_ab_qgq_old(order,y12,y15,y25,soft1gq,soft2gq)
!      call soft_ab_qag_old(order,y12,y15,y25,soft1qa,soft2qa)
c      call soft_ab_qgq(order,y25,y12,y15,soft1qg,soft2qg)
c      call soft_ab_qgq(order,y15,y12,y25,soft1gq,soft2gq)
c      call soft_ab_qgq(order,y25,y15,y12,soft1qa,soft2qa)

!      if (order > 1) then
!        call soft_nab_ggg_old(order,y12,y25,y15,soft2gg_nab)
!        call soft_nab_qgq_old(order,y12,y25,y15,soft2qg_nab)
!        call soft_nab_qgq_old(order,y12,y15,y25,soft2gq_nab)
!        call soft_nab_qag_old(order,y12,y15,y25,soft2qa_nab)
!        soft2gg(:)=soft2gg(:)+soft2gg_nab(:)
!        soft2qg(:)=soft2qg(:)+soft2qg_nab(:)
!        soft2gq(:)=soft2gq(:)+soft2gq_nab(:)
!        soft2qa(:)=soft2qa(:)+soft2qa_nab(:)
!      endif
      
!      oldsoft2gg(:)=soft2gg(:)
!      oldsoft2qg(:)=soft2qg(:)
!      oldsoft2gq(:)=soft2gq(:)
!      oldsoft2qa(:)=soft2qa(:)

      call computeIijm(y12,y25,y15,Iijm)
      call soft_ab_ggg(order,y12,y25,y15,Iijm,1,2,3,4,5,3,soft1gg,soft2gg)
      call soft_ab_qgq(order,y12,y25,y15,Iijm,1,2,3,4,5,3,soft1qg,soft2qg)
      call soft_ab_qgq(order,y12,y15,y25,Iijm,2,1,3,5,4,3,soft1gq,soft2gq)
      call soft_ab_qag(order,y12,y25,y15,Iijm,1,2,3,4,5,3,soft1qa,soft2qa)
      if (order > 1) then
        call soft_nab_ggg(order,y12,y25,y15,Iijm,1,2,3,4,5,3,soft2gg_nab)
        call soft_nab_qgq(order,y12,y25,y15,Iijm,1,2,3,4,5,3,soft2qg_nab)
        call soft_nab_qgq(order,y12,y15,y25,Iijm,2,1,3,5,4,3,soft2gq_nab)
        call soft_nab_qag(order,y12,y25,y15,Iijm,1,2,3,4,5,3,soft2qa_nab)
        soft2gg(:)=soft2gg(:)+soft2gg_nab(:)
        soft2qg(:)=soft2qg(:)+soft2qg_nab(:)
        soft2gq(:)=soft2gq(:)+soft2gq_nab(:)
        soft2qa(:)=soft2qa(:)+soft2qa_nab(:)
      endif

!      write(6,*)'soft2gg',soft2gg/oldsoft2gg
!      write(6,*)'soft2qg',soft2qg/oldsoft2qg
!      write(6,*)'soft2gq',soft2gq/oldsoft2gq
!      write(6,*)'soft2qa',soft2qa/oldsoft2qa
!      pause

!      soft2qa(-1)=zip
!      jet2g(-1)=zip
!      Hqqb(2)=zip
!      beama2(:,-1)=zip
!      beamb2(:,-1)=zip
      
      if (dynamictau) then
        tauc=getdynamictau(p)
      else
        tauc=taucut
      endif

      xmsq=zip
      do j=-nf,nf
      do k=-nf,nf
      
      if ((j .ne. 0) .and. (k .ne. 0)
     & .and. (j .ne. -k)) cycle
      
      if (((j == 0) .and. (k == 0))
     & .or. (j*k .ne. 0)) then !  gg -> g or qa -> g => jet is a gluon
        if ((j == 0) .and. (k == 0)) then
          bit=assemblejet(order,tauc,
     &     beama0(j),beamb0(k),beama1(j,:),beamb1(k,:),
     &     beama2(j,:),beamb2(k,:),soft1gg,soft2gg,jet1g,jet2g,Hgg)
        else
          bit=assemblejet(order,tauc,
     &     beama0(j),beamb0(k),beama1(j,:),beamb1(k,:),
     &     beama2(j,:),beamb2(k,:),soft1qa,soft2qa,jet1g,jet2g,Hqqb)
        endif
      else                     ! qg -> q  => jet is a quark
        if (j == 0) then
          bit=assemblejet(order,tauc,
     &     beama0(j),beamb0(k),beama1(j,:),beamb1(k,:),
     &     beama2(j,:),beamb2(k,:),soft1gq,soft2gq,jet1q,jet2q,Hgqb)
        elseif (k == 0) then
          bit=assemblejet(order,tauc,
     &     beama0(j),beamb0(k),beama1(j,:),beamb1(k,:),
     &     beama2(j,:),beamb2(k,:),soft1qg,soft2qg,jet1q,jet2q,Hqg)
        endif
      endif
      
      if     ((j == 0) .and. (k == 0)) then
        bit=bit*msqgg
      elseif (j == 0) then
        bit=bit*msqgq
      elseif (k == 0) then
        bit=bit*msqqg
      elseif (j == -k) then
        bit=bit*msqqa
      else
        bit=zip
      endif

      xmsq=xmsq+bit
c      write(6,*) 'j,k,bit',j,k,bit

      enddo
      enddo
c      pause

      return
      end
