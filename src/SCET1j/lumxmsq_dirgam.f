      subroutine lumxmsq_dirgam(p,xx,z1,z2,QB,order,xmsq)
      implicit none
      include 'types.f'
C-----Authors: John Campbell, Keith Ellis
C-----June 2000 and December 2001
c----Matrix element for Z production
C----averaged over initial colours and spins
c     q(-p1)+qbar(-p2)-->(e^-(p3)+e^+(p4))+g(p5)
c---
       include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'zcouple.f'
      include 'ewcharge.f'
      include 'qcdcouple.f'
      include 'scale.f'
      include 'facscale.f'
      include 'ckm.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      include 'scet_const.f'
      include 'taucut.f'
      integer:: j,k,ih1,ih2,m,n,order
      real(dp),intent(out):: xmsq
      real(dp):: p(mxpart,4),fac,facgg,xx(2),hard(2),
     & soft1qg(-1:1),soft2qg(-1:3),soft2qg_nab(-1:3),
     & soft1gq(-1:1),soft2gq(-1:3),soft2gq_nab(-1:3),
     & soft1qa(-1:1),soft2qa(-1:3),soft2qa_nab(-1:3),
     & beama0(-5:5),beamb0(-5:5),tauc,getdynamictau,
     & beama1(-5:5,-1:1),beamb1(-5:5,-1:1),
     & beama2(-5:5,-1:3),beamb2(-5:5,-1:3),
     & jet1q(-1:1),jet2q(-1:3),jet1g(-1:1),jet2g(-1:3),
     & z1,z2,QB(2),lum0,lum1(-1:1),lum2(-1:3),bit,
     & msqhard1(-nf:nf,-nf:nf),assemblejet,y12,y14,y24,Iijm(5)
      real(dp):: msq0(-nf:nf,-nf:nf),msq1(-nf:nf,-nf:nf),msq2(-nf:nf,-nf:nf)
      common/density/ih1,ih2
      
      msq0(:,:)=zip
      msq1(:,:)=zip
      msq2(:,:)=zip
      call qqb_dirgam_vbis(p,msq0,0)
      if (order > 0) then
        call qqb_dirgam_vbis(p,msq1,1)
      endif
      if (order > 1) then
        call qqb_dirgam_vbis(p,msq2,2)
      endif

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
      call jetq(order,two*p(4,4),jet1q,jet2q)
      call jetg(order,two*p(4,4),jet1g,jet2g)

c---- soft functions
      y12=s(1,2)/p(1,4)/p(2,4)/four
      y14=s(1,4)/p(1,4)/p(4,4)/four
      y24=s(2,4)/p(2,4)/p(4,4)/four
      
      call computeIijm(y12,y24,y14,Iijm)
      call soft_ab_qgq(order,y12,y24,y14,Iijm,1,2,3,4,5,3,soft1qg,soft2qg)
      call soft_ab_qgq(order,y12,y14,y24,Iijm,2,1,3,5,4,3,soft1gq,soft2gq)
      call soft_ab_qag(order,y12,y14,y24,Iijm,1,2,3,4,5,3,soft1qa,soft2qa)
      if (order > 1) then
        call soft_nab_qgq(order,y12,y24,y14,Iijm,1,2,3,4,5,3,soft2qg_nab)
        call soft_nab_qgq(order,y12,y14,y24,Iijm,2,1,3,5,4,3,soft2gq_nab)
        call soft_nab_qag(order,y12,y14,y24,Iijm,1,2,3,4,5,3,soft2qa_nab)
        soft2qg(:)=soft2qg(:)+soft2qg_nab(:)
        soft2gq(:)=soft2gq(:)+soft2gq_nab(:)
        soft2qa(:)=soft2qa(:)+soft2qa_nab(:)
      endif

      if (dynamictau) then
        tauc=getdynamictau(p)
      else
        tauc=taucut
      endif

      xmsq=zip
      do j=-nf,nf
      do k=-nf,nf
      if( j .ne. 0 .and. k .ne. 0 .and. j .ne. -k) cycle

      if (order > 0) then
        hard(1)=msq1(j,k)/msq0(j,k)/ason2pi
      else
        hard(1)=zip
      endif
      if (order > 1) then
        hard(2)=msq2(j,k)/msq0(j,k)/ason2pi**2
      else
        hard(2)=zip
      endif
      
      if ((j == 0) .and. (k == 0)) then
c--- gg->gag only has a hard(2) contribution
        bit=beama0(j)*beamb0(k)*msq2(j,k)
        goto 20
      elseif (j*k == 0) then ! initial state gluon -> jet is a quark
        if (j == 0) then
          bit=assemblejet(order,tauc,
     &     beama0(j),beamb0(k),beama1(j,:),beamb1(k,:),
     &     beama2(j,:),beamb2(k,:),soft1gq,soft2gq,jet1q,jet2q,hard)
        elseif (k == 0) then
          bit=assemblejet(order,tauc,
     &     beama0(j),beamb0(k),beama1(j,:),beamb1(k,:),
     &     beama2(j,:),beamb2(k,:),soft1qg,soft2qg,jet1q,jet2q,hard)
        endif
      else               ! initial state quarks -> jet is a gluon
          bit=assemblejet(order,tauc,
     &     beama0(j),beamb0(k),beama1(j,:),beamb1(k,:),
     &     beama2(j,:),beamb2(k,:),soft1qa,soft2qa,jet1g,jet2g,hard)
      endif
      
      bit=bit*msq0(j,k)

 20   continue

      xmsq=xmsq+bit

      enddo
      enddo
      
      return
      end
