      subroutine lumxmsq_w1jet(p,xx,z1,z2,QB,order,xmsq)
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
      include 'sprods_com.f'
      include 'zprods_decl.f'
      include 'scet_const.f'
      include 'taucut.f'
      integer:: j,k,ih1,ih2,m,n,order
      real(dp),intent(out):: xmsq
      real(dp):: p(mxpart,4),fac,xx(2),hard(2),
     & soft1qg(-1:1),soft2qg(-1:3),soft2qg_nab(-1:3),
     & soft1gq(-1:1),soft2gq(-1:3),soft2gq_nab(-1:3),
     & soft1qa(-1:1),soft2qa(-1:3),soft2qa_nab(-1:3),
     & beama0(-5:5),beamb0(-5:5),tauc,getdynamictau,
     & beama1(-5:5,-1:1),beamb1(-5:5,-1:1),
     & beama2(-5:5,-1:3),beamb2(-5:5,-1:3),
     & jet1q(-1:1),jet2q(-1:3),jet1g(-1:1),jet2g(-1:3),
     & z1,z2,QB(2),lum0,lum1(-1:1),lum2(-1:3),bit,propsq,
     & msqhard1(-nf:nf,-nf:nf),assemblejet,y12,y15,y25,Iijm(5)
      real(dp):: qqbWg,qbqWg,qgWq,qbgWqb,gqbWqb,gqWq
      real(dp):: Hqqb(2),Hqg(2),Hgqb(2),Hqbq(2),Hqbg(2),Hgq(2)
      common/density/ih1,ih2

c--- compute hard function
      call spinoru(5,p,za,zb)
c--- calculate the propagator
      propsq=(s(3,4)-wmass**2)**2+(wmass*wwidth)**2
      fac=gwsq**2*gsq*V/propsq

c--- compute matrix elements
      call Wampqqbgsq(order,2,1,5,3,4,za,zb,qqbWg,Hqqb(1),Hqqb(2))
      call Wampqqbgsq(order,5,1,2,3,4,za,zb,qgWq,Hqg(1),Hqg(2))
      call Wampqqbgsq(order,2,5,1,3,4,za,zb,gqbWqb,Hgqb(1),Hgqb(2))

      call Wampqqbgsq(order,1,2,5,3,4,za,zb,qbqWg,Hqbq(1),Hqbq(2))
      call Wampqqbgsq(order,1,5,2,3,4,za,zb,qbgWqb,Hqbg(1),Hqbg(2))
      call Wampqqbgsq(order,5,2,1,3,4,za,zb,gqWq,Hgq(1),Hgq(2))

c--- apply overall factor
      qqbWg= aveqq*fac*qqbWg
      qgWq  =aveqg*fac*qgWq
      gqbWqb=aveqg*fac*gqbWqb
      qbqWg =aveqq*fac*qbqWg
      qbgWqb=aveqg*fac*qbgWqb
      gqWq  =aveqg*fac*gqWq

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
      call jetq(order,two*p(5,4),jet1q,jet2q)
      call jetg(order,two*p(5,4),jet1g,jet2g)

c---- soft functions
      y12=s(1,2)/p(1,4)/p(2,4)/four
      y15=s(1,5)/p(1,4)/p(5,4)/four
      y25=s(2,5)/p(2,4)/p(5,4)/four
c      call soft_ab_qgq(order,y25,y12,y15,soft1qg,soft2qg)
c      call soft_ab_qgq(order,y15,y12,y25,soft1gq,soft2gq)
c      call soft_ab_qgq(order,y25,y15,y12,soft1qa,soft2qa)

!      call soft_ab_qgq(order,y12,y25,y15,soft1qg,soft2qg)
!      call soft_ab_qgq(order,y12,y15,y25,soft1gq,soft2gq)
!      call soft_ab_qag(order,y12,y15,y25,soft1qa,soft2qa)
!      if (order > 1) then
!        call soft_nab_qgq(order,y12,y25,y15,soft2qg_nab)
!        call soft_nab_qgq(order,y12,y15,y25,soft2gq_nab)
!        call soft_nab_qag(order,y12,y15,y25,soft2qa_nab)
!        soft2qg(:)=soft2qg(:)+soft2qg_nab(:)
!        soft2gq(:)=soft2gq(:)+soft2gq_nab(:)
!        soft2qa(:)=soft2qa(:)+soft2qa_nab(:)
!      endif
      
      call computeIijm(y12,y25,y15,Iijm)
      call soft_ab_qgq(order,y12,y25,y15,Iijm,1,2,3,4,5,3,soft1qg,soft2qg)
      call soft_ab_qgq(order,y12,y15,y25,Iijm,2,1,3,5,4,3,soft1gq,soft2gq)
      call soft_ab_qag(order,y12,y15,y25,Iijm,1,2,3,4,5,3,soft1qa,soft2qa)
      if (order > 1) then
        call soft_nab_qgq(order,y12,y25,y15,Iijm,1,2,3,4,5,3,soft2qg_nab)
        call soft_nab_qgq(order,y12,y15,y25,Iijm,2,1,3,5,4,3,soft2gq_nab)
        call soft_nab_qag(order,y12,y15,y25,Iijm,1,2,3,4,5,3,soft2qa_nab)
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
      if ((j == 0) .and. (k == 0)) cycle
      if (j*k > 0) cycle
      
      if (j*k == 0) then ! initial state gluon -> jet is a quark
        if (j == 0) then
          if (k > 0) then
            hard(:)=Hgq(:)
          else
            hard(:)=Hgqb(:)
          endif
          bit=assemblejet(order,tauc,
     &     beama0(j),beamb0(k),beama1(j,:),beamb1(k,:),
     &     beama2(j,:),beamb2(k,:),soft1gq,soft2gq,jet1q,jet2q,hard)
        elseif (k == 0) then
          if (j > 0) then
            hard(:)=Hqg(:)
          else
            hard(:)=Hqbg(:)
          endif
          bit=assemblejet(order,tauc,
     &     beama0(j),beamb0(k),beama1(j,:),beamb1(k,:),
     &     beama2(j,:),beamb2(k,:),soft1qg,soft2qg,jet1q,jet2q,hard)
        endif
      else               ! initial state quarks -> jet is a gluon
          if (j > 0) then
            hard(:)=Hqqb(:)
          else
            hard(:)=Hqbq(:)
          endif
          bit=assemblejet(order,tauc,
     &     beama0(j),beamb0(k),beama1(j,:),beamb1(k,:),
     &     beama2(j,:),beamb2(k,:),soft1qa,soft2qa,jet1g,jet2g,hard)
      endif
      
      if     ((j > 0) .and. (k < 0)) then
        bit=bit*Vsq(j,k)*qqbWg
      elseif ((j < 0) .and. (k > 0)) then
        bit=bit*Vsq(j,k)*qbqWg
      elseif ((j > 0) .and. (k == 0)) then
        bit=bit*
     &   (Vsq(j,-1)+Vsq(j,-2)+Vsq(j,-3)+Vsq(j,-4)+Vsq(j,-5))*qgWq
      elseif ((j < 0) .and. (k == 0)) then
        bit=bit*
     &    (Vsq(j,+1)+Vsq(j,+2)+Vsq(j,+3)+Vsq(j,+4)+Vsq(j,+5))*qbgWqb
      elseif ((j == 0) .and. (k > 0)) then
        bit=bit*
     &    (Vsq(-1,k)+Vsq(-2,k)+Vsq(-3,k)+Vsq(-4,k)+Vsq(-5,k))*gqWq
      elseif ((j == 0) .and. (k < 0)) then
        bit=bit*
     &    (Vsq(+1,k)+Vsq(+2,k)+Vsq(+3,k)+Vsq(+4,k)+Vsq(+5,k))*gqbWqb
      else
        bit=zip
      endif

      xmsq=xmsq+bit

      enddo
      enddo
      
      return
      end
