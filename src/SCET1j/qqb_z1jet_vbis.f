      subroutine qqb_z1jet_vbis(p,msq,order)
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
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'zcouple.f'
      include 'ewcharge.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      include 'nflav.f'
      integer:: j,k,order
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),fac,facgg
      complex(dp):: prop
c      real(dp):: 
c     & AqqbZg2(2,2),AqbqZg2(2,2),AqgZq2(2,2),AqbgZqb2(2,2)
c     & AgqbZqb2(2,2),AgqZq2(2,2)
      real(dp):: qqbZg(2,2),qbqZg(2,2),qgZq(2,2),
     & qbgZqb(2,2),gqbZqb(2,2),gqZq(2,2),ggZg(2)
      real(dp)::
     & Hqqb1(2,2),Hqg1(2,2),Hgqb1(2,2),Hqbq1(2,2),Hqbg1(2,2),Hgq1(2,2),
     & Hqqb2(2,2),Hqg2(2,2),Hgqb2(2,2),Hqbq2(2,2),Hqbg2(2,2),Hgq2(2,2),
     & Hqqb2NFV(2,2),Hqg2NFV(2,2),Hgqb2NFV(2,2),
     & Hqbq2NFV(2,2),Hqbg2NFV(2,2),Hgq2NFV(2,2),Qsumg,QsumZ

      msq(:,:)=zip
      
      call spinoru(5,p,za,zb)
      
      prop=s(3,4)/cmplx((s(3,4)-zmass**2),zmass*zwidth,dp)
      fac=four*V*esq**2*gsq
      facgg=(40._dp/3._dp)*esq**2*gsq

c--- remove photon contribution for debugging with/without decay
c      q1=zip

c--- compute matrix elements
c--- for Zggg, note that 4,3 are flipped to reflect different notation in amp
      call Zampgggsq(order,1,2,5,4,3,za,zb,ggZg)

      call Zampqqbgsq(order,2,1,5,3,4,zip,za,zb,qqbZg,Hqqb1,Hqqb2)
      call Zampqqbgsq(order,5,1,2,3,4,zip,za,zb,qgZq,Hqg1,Hqg2)
      call Zampqqbgsq(order,2,5,1,3,4,zip,za,zb,gqbZqb,Hgqb1,Hgqb2)

      call Zampqqbgsq(order,1,2,5,3,4,zip,za,zb,qbqZg,Hqbq1,Hqbq2)
      call Zampqqbgsq(order,1,5,2,3,4,zip,za,zb,qbgZqb,Hqbg1,Hqbg2)
      call Zampqqbgsq(order,5,2,1,3,4,zip,za,zb,gqZq,Hgq1,Hgq2)

c--- extract term proportional to NFV piece of 2-loop calculation
      if (order == 2) then
      call Zampqqbgsq(order,2,1,5,3,4,one,za,zb,qqbZg,Hqqb1,Hqqb2NFV)
      call Zampqqbgsq(order,5,1,2,3,4,one,za,zb,qgZq,Hqg1,Hqg2NFV)
      call Zampqqbgsq(order,2,5,1,3,4,one,za,zb,gqbZqb,Hgqb1,Hgqb2NFV)

      call Zampqqbgsq(order,1,2,5,3,4,one,za,zb,qbqZg,Hqbq1,Hqbq2NFV)
      call Zampqqbgsq(order,1,5,2,3,4,one,za,zb,qbgZqb,Hqbg1,Hqbg2NFV)
      call Zampqqbgsq(order,5,2,1,3,4,one,za,zb,gqZq,Hgq1,Hgq2NFV)

      Hqqb2NFV(:,:)=Hqqb2NFV(:,:) - Hqqb2(:,:)
      Hqg2NFV(:,:) =Hqg2NFV(:,:)  - Hqg2(:,:) 
      Hgqb2NFV(:,:)=Hgqb2NFV(:,:) - Hgqb2(:,:)
      Hqbq2NFV(:,:)=Hqbq2NFV(:,:) - Hqbq2(:,:)
      Hqbg2NFV(:,:)=Hqbg2NFV(:,:) - Hqbg2(:,:)
      Hgq2NFV(:,:) =Hgq2NFV(:,:)  - Hgq2(:,:) 
       
      else
       
      Hqqb2NFV(:,:)=zip
      Hqg2NFV(:,:)=zip
      Hgqb2NFV(:,:)=zip
      Hqbq2NFV(:,:)=zip
      Hqbg2NFV(:,:)=zip
      Hgq2NFV(:,:)=zip 
       
      endif
      
c--- apply overall factor
      qqbZg= aveqq*fac*qqbZg
      qgZq  =aveqg*fac*qgZq
      gqbZqb=aveqg*fac*gqbZqb
      qbqZg =aveqq*fac*qbqZg
      qbgZqb=aveqg*fac*qbgZqb
      gqZq  =aveqg*fac*gqZq

      ggZg(:)=avegg*facgg*(ason2pi)**2*ggZg(:)
      
c--- pick out higher-order coefficients
      if (order == 1) then
        qqbZg(:,:) =ason2pi*Hqqb1(:,:)*qqbZg(:,:)
        qgZq(:,:)  =ason2pi*Hqg1(:,:)*qgZq(:,:)
        gqbZqb(:,:)=ason2pi*Hgqb1(:,:)*gqbZqb(:,:)
        qbqZg(:,:) =ason2pi*Hqbq1(:,:)*qbqZg(:,:)
        qbgZqb(:,:)=ason2pi*Hqbg1(:,:)*qbgZqb(:,:)
        gqZq(:,:)  =ason2pi*Hgq1(:,:)*gqZq(:,:)
      endif
      if (order == 2) then
        Hqqb2NFV(:,:)=ason2pi**2*Hqqb2NFV(:,:)*qqbZg(:,:)
        Hqg2NFV(:,:) =ason2pi**2*Hqg2NFV(:,:)*qgZq(:,:)
        Hgqb2NFV(:,:)=ason2pi**2*Hgqb2NFV(:,:)*gqbZqb(:,:)
        Hqbq2NFV(:,:)=ason2pi**2*Hqbq2NFV(:,:)*qbqZg(:,:)
        Hqbg2NFV(:,:)=ason2pi**2*Hqbg2NFV(:,:)*qbgZqb(:,:)
        Hgq2NFV(:,:) =ason2pi**2*Hgq2NFV(:,:)*gqZq(:,:)
       
        qqbZg(:,:) =ason2pi**2*Hqqb2(:,:)*qqbZg(:,:)
        qgZq(:,:)  =ason2pi**2*Hqg2(:,:)*qgZq(:,:)
        gqbZqb(:,:)=ason2pi**2*Hgqb2(:,:)*gqbZqb(:,:)
        qbqZg(:,:) =ason2pi**2*Hqbq2(:,:)*qbqZg(:,:)
        qbgZqb(:,:)=ason2pi**2*Hqbg2(:,:)*qbgZqb(:,:)
        gqZq(:,:)  =ason2pi**2*Hgq2(:,:)*gqZq(:,:)

      endif

! couplings for gg->Zg according to Eqs.(2.19) and (2.20) of 1302.2630;
! note that this is also the coupling for the NFV pieces as well,
!  c.f. Eqs. (3.8) and (3.9) of 1112.1531
c--- five light quark flavors are included in loop
      Qsumg=2._dp*Qu+3._dp*Qd
      QsumZ=half*(2._dp*(L(2)+R(2))+3._dp*(L(1)+R(1)))
      
      do j=-nflav,nflav
      do k=-nflav,nflav

      if( (j .ne. 0) .and. (k .ne. 0) .and. (j .ne. -k)) goto 20

      if     ((j == 0) .and. (k == 0)) then
          msq(j,k)=abs(Qsumg*q1+QsumZ*l1*prop)**2*ggZg(1)
     &            +abs(Qsumg*q1+QsumZ*r1*prop)**2*ggZg(2)
          msq(j,k)=msq(j,k)/s(3,4)**2
      elseif ((j > 0) .and. (k < 0)) then
          msq(j,k)=abs(Q(j)*q1+L(j)*l1*prop)**2*qqbZg(1,1)
     &            +abs(Q(j)*q1+L(j)*r1*prop)**2*qqbZg(1,2)
     &            +abs(Q(j)*q1+R(j)*l1*prop)**2*qqbZg(2,1)
     &            +abs(Q(j)*q1+R(j)*r1*prop)**2*qqbZg(2,2)
          
          msq(j,k)=msq(j,k)
     &     +real(conjg(Q(j)*q1+L(j)*l1*prop)*(Qsumg*q1+QsumZ*l1*prop),dp)*Hqqb2NFV(1,1)
     &     +real(conjg(Q(j)*q1+L(j)*r1*prop)*(Qsumg*q1+QsumZ*r1*prop),dp)*Hqqb2NFV(1,2)
     &     +real(conjg(Q(j)*q1+R(j)*l1*prop)*(Qsumg*q1+QsumZ*l1*prop),dp)*Hqqb2NFV(2,1)
     &     +real(conjg(Q(j)*q1+R(j)*r1*prop)*(Qsumg*q1+QsumZ*r1*prop),dp)*Hqqb2NFV(2,2)
     
          msq(j,k)=msq(j,k)/s(3,4)**2
      elseif ((j < 0) .and. (k > 0)) then
          msq(j,k)=abs(Q(k)*q1+L(k)*l1*prop)**2*qbqZg(1,1)
     &            +abs(Q(k)*q1+L(k)*r1*prop)**2*qbqZg(1,2)
     &            +abs(Q(k)*q1+R(k)*l1*prop)**2*qbqZg(2,1)
     &            +abs(Q(k)*q1+R(k)*r1*prop)**2*qbqZg(2,2)
          
          msq(j,k)=msq(j,k)
     &     +real(conjg(Q(k)*q1+L(k)*l1*prop)*(Qsumg*q1+QsumZ*l1*prop),dp)*Hqbq2NFV(1,1)
     &     +real(conjg(Q(k)*q1+L(k)*r1*prop)*(Qsumg*q1+QsumZ*r1*prop),dp)*Hqbq2NFV(1,2)
     &     +real(conjg(Q(k)*q1+R(k)*l1*prop)*(Qsumg*q1+QsumZ*l1*prop),dp)*Hqbq2NFV(2,1)
     &     +real(conjg(Q(k)*q1+R(k)*r1*prop)*(Qsumg*q1+QsumZ*r1*prop),dp)*Hqbq2NFV(2,2)
     
          msq(j,k)=msq(j,k)/s(3,4)**2
      elseif ((j > 0) .and. (k == 0)) then
          msq(j,k)=abs(Q(j)*q1+L(j)*l1*prop)**2*qgZq(1,1)
     &            +abs(Q(j)*q1+L(j)*r1*prop)**2*qgZq(1,2)
     &            +abs(Q(j)*q1+R(j)*l1*prop)**2*qgZq(2,1)
     &            +abs(Q(j)*q1+R(j)*r1*prop)**2*qgZq(2,2)
          
          msq(j,k)=msq(j,k)
     &     +real(conjg(Q(j)*q1+L(j)*l1*prop)*(Qsumg*q1+QsumZ*l1*prop),dp)*Hqg2NFV(1,1)
     &     +real(conjg(Q(j)*q1+L(j)*r1*prop)*(Qsumg*q1+QsumZ*r1*prop),dp)*Hqg2NFV(1,2)
     &     +real(conjg(Q(j)*q1+R(j)*l1*prop)*(Qsumg*q1+QsumZ*l1*prop),dp)*Hqg2NFV(2,1)
     &     +real(conjg(Q(j)*q1+R(j)*r1*prop)*(Qsumg*q1+QsumZ*r1*prop),dp)*Hqg2NFV(2,2)
     
          msq(j,k)=msq(j,k)/s(3,4)**2
      elseif ((j < 0) .and. (k == 0)) then
          msq(j,k)=abs(Q(-j)*q1+L(-j)*l1*prop)**2*qbgZqb(1,1)
     &            +abs(Q(-j)*q1+L(-j)*r1*prop)**2*qbgZqb(1,2)
     &            +abs(Q(-j)*q1+R(-j)*l1*prop)**2*qbgZqb(2,1)
     &            +abs(Q(-j)*q1+R(-j)*r1*prop)**2*qbgZqb(2,2)
          
          msq(j,k)=msq(j,k)
     &     +real(conjg(Q(-j)*q1+L(-j)*l1*prop)*(Qsumg*q1+QsumZ*l1*prop),dp)*Hqbg2NFV(1,1)
     &     +real(conjg(Q(-j)*q1+L(-j)*r1*prop)*(Qsumg*q1+QsumZ*r1*prop),dp)*Hqbg2NFV(1,2)
     &     +real(conjg(Q(-j)*q1+R(-j)*l1*prop)*(Qsumg*q1+QsumZ*l1*prop),dp)*Hqbg2NFV(2,1)
     &     +real(conjg(Q(-j)*q1+R(-j)*r1*prop)*(Qsumg*q1+QsumZ*r1*prop),dp)*Hqbg2NFV(2,2)
     
          msq(j,k)=msq(j,k)/s(3,4)**2
      elseif ((j == 0) .and. (k > 0)) then
          msq(j,k)=abs(Q(k)*q1+L(k)*l1*prop)**2*gqZq(1,1)
     &            +abs(Q(k)*q1+L(k)*r1*prop)**2*gqZq(1,2)
     &            +abs(Q(k)*q1+R(k)*l1*prop)**2*gqZq(2,1)
     &            +abs(Q(k)*q1+R(k)*r1*prop)**2*gqZq(2,2)
          
          msq(j,k)=msq(j,k)
     &     +real(conjg(Q(k)*q1+L(k)*l1*prop)*(Qsumg*q1+QsumZ*l1*prop),dp)*Hgq2NFV(1,1)
     &     +real(conjg(Q(k)*q1+L(k)*r1*prop)*(Qsumg*q1+QsumZ*r1*prop),dp)*Hgq2NFV(1,2)
     &     +real(conjg(Q(k)*q1+R(k)*l1*prop)*(Qsumg*q1+QsumZ*l1*prop),dp)*Hgq2NFV(2,1)
     &     +real(conjg(Q(k)*q1+R(k)*r1*prop)*(Qsumg*q1+QsumZ*r1*prop),dp)*Hgq2NFV(2,2)
     
          msq(j,k)=msq(j,k)/s(3,4)**2
      elseif ((j == 0) .and. (k < 0)) then
          msq(j,k)=abs(Q(-k)*q1+L(-k)*l1*prop)**2*gqbZqb(1,1)
     &            +abs(Q(-k)*q1+L(-k)*r1*prop)**2*gqbZqb(1,2)
     &            +abs(Q(-k)*q1+R(-k)*l1*prop)**2*gqbZqb(2,1)
     &            +abs(Q(-k)*q1+R(-k)*r1*prop)**2*gqbZqb(2,2)
          
          msq(j,k)=msq(j,k)
     &     +real(conjg(Q(-k)*q1+L(-k)*l1*prop)*(Qsumg*q1+QsumZ*l1*prop),dp)*Hgqb2NFV(1,1)
     &     +real(conjg(Q(-k)*q1+L(-k)*r1*prop)*(Qsumg*q1+QsumZ*r1*prop),dp)*Hgqb2NFV(1,2)
     &     +real(conjg(Q(-k)*q1+R(-k)*l1*prop)*(Qsumg*q1+QsumZ*l1*prop),dp)*Hgqb2NFV(2,1)
     &     +real(conjg(Q(-k)*q1+R(-k)*r1*prop)*(Qsumg*q1+QsumZ*r1*prop),dp)*Hgqb2NFV(2,2)
     
          msq(j,k)=msq(j,k)/s(3,4)**2
      endif

 20   continue
      enddo
      enddo
      return
      end
 

