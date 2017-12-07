      subroutine qqb_w1jet_vbis(p,msq,order)
      implicit none
      include 'types.f'
c----Matrix element for W production at O[alphas^(1+order)]
C----averaged over initial colours and spins
C for nwz=+1
c     u(-p1)+dbar(-p2)--> W^+(n(p3)+e^+(p4))   + g(p5)
C For nwz=-1
c     d(-p1)+ubar(-p2)--> W^-(e^-(p3)+nbar(p4))+ g(p5) 
c---
c---- Computed using amplitudes extracted from PeTer
      include 'constants.f'
      include 'mxpart.f'
      include 'nf.f'
      include 'ckm.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'masses.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      integer:: j,k,order
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),fac,propsq
      real(dp):: qqbWg,qbqWg,qgWq,qbgWqb,gqbWqb,gqWq
      real(dp):: Hqqb(2),Hqg(2),Hgqb(2),Hqbq(2),Hqbg(2),Hgq(2)

      msq(:,:)=zip

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

c--- pick out higher-order coefficients
      if (order == 1) then
        qqbWg =ason2pi*Hqqb(1)*qqbWg
        qgWq  =ason2pi*Hqg(1)*qgWq
        gqbWqb=ason2pi*Hgqb(1)*gqbWqb
        qbqWg =ason2pi*Hqbq(1)*qbqWg
        qbgWqb=ason2pi*Hqbg(1)*qbgWqb
        gqWq  =ason2pi*Hgq(1)*gqWq
      endif
      if (order == 2) then
        qqbWg =ason2pi**2*Hqqb(2)*qqbWg
        qgWq  =ason2pi**2*Hqg(2)*qgWq
        gqbWqb=ason2pi**2*Hgqb(2)*gqbWqb
        qbqWg =ason2pi**2*Hqbq(2)*qbqWg
        qbgWqb=ason2pi**2*Hqbg(2)*qbgWqb
        gqWq  =ason2pi**2*Hgq(2)*gqWq
      endif

c--- numerical cross-check
c      write(6,*) 'Hqqb(1)/twopi,Hqqb(2)/twopi**2',Hqqb(1)/twopi,Hqqb(2)/twopi**2
c      write(6,*) 'Hqg(1)/twopi, Hqg(2)/twopi**2 ',Hqg(1)/twopi, Hqg(2)/twopi**2

      do j=-nf,nf
      do k=-nf,nf
      if     ((j > 0) .and. (k < 0)) then
          msq(j,k)=Vsq(j,k)*qqbWg
      elseif ((j < 0) .and. (k > 0)) then
          msq(j,k)=Vsq(j,k)*qbqWg
      elseif ((j > 0) .and. (k == 0)) then
          msq(j,k)=
     &   (Vsq(j,-1)+Vsq(j,-2)+Vsq(j,-3)+Vsq(j,-4)+Vsq(j,-5))*qgWq
      elseif ((j < 0) .and. (k == 0)) then
          msq(j,k)=
     &    (Vsq(j,+1)+Vsq(j,+2)+Vsq(j,+3)+Vsq(j,+4)+Vsq(j,+5))*qbgWqb
      elseif ((j == 0) .and. (k > 0)) then
          msq(j,k)=
     &    (Vsq(-1,k)+Vsq(-2,k)+Vsq(-3,k)+Vsq(-4,k)+Vsq(-5,k))*gqWq
      elseif ((j == 0) .and. (k < 0)) then
          msq(j,k)=
     &    (Vsq(+1,k)+Vsq(+2,k)+Vsq(+3,k)+Vsq(+4,k)+Vsq(+5,k))*gqbWqb
      endif

      enddo
      enddo
      
      return
      end
 
