      subroutine qqb_gam2j_gvec(p,n,in,msqv)
! Wrapper routine to qqb_gam2j_gvecx_new below
! that only returns msqv (and msqv_cs via common block)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'zcouple.f'
      include 'ewcharge.f'
      include 'pp.f'
      include 'msqv_cs.f'
C ip is the label of the emitter
C in is the label of the contracted line
      integer:: in
      real(dp):: msqv(-nf:nf,-nf:nf),p(mxpart,4),n(4)
      real(dp):: ppmsqvx(ppmax),fillmsqv_cs(0:2,-nf:nf,-nf:nf)
!     & ,fillmsqv(-nf:nf,-nf:nf)

      msqv(:,:)=zip
      msqv_cs(:,:,:)=zip

      call qqb_gam2j_gvecx_new(p,n,in,msqv,fillmsqv_cs,ppmsqvx)

!      msqv(:,:)=fillmsqv(:,:)
      msqv_cs(:,:,:)=fillmsqv_cs(:,:,:)
      
      return
      end
      
      


      subroutine qqb_gam2j_gvecx_new(p,n,in,msqv,msqv_cs,ppmsqvx)
      implicit none
      include 'types.f'
C-----June 2016: checked that gives the right matrix element
C     when summed over polarizations.

c----Matrix element for gamma+2jet production
C----averaged over initial colours and spins
c    line 6 contracted with the vector n(mu)
c     q(-p1)+qbar(-p2)--> ga(p3) + g(p4) + g(p5)

c---- Note that this routine is the counterpart of qqb_gam2j_gvec
c---- and should be used when we are calculating the 4Q contribution.
c---- This routine is necessary because it returns the contribution
c---- for the correct ordering of particles in the final state - for
c---- example, gg -> qqb as well as gg -> qbq

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'zcouple.f'
      include 'ewcharge.f'
      include 'zprods_decl.f'
      include 'pp.f'

C ip is the label of the emitter
C in is the label of the contracted line
      integer:: j,k,nquark,in,icol
      real(dp):: fac,n(4)
      complex(dp):: zab(mxpart,mxpart),zba(mxpart,mxpart)
!      complex(dp):: amp1(2,2,2),amp2(2,2,2),amplo(2,2,2,2),amp1s(2,2,2),amp2s(2,2,2),amplos(2,2,2,2)
      real(dp):: msqv(-nf:nf,-nf:nf),p(mxpart,4),
     & qqb_gagg(0:2),qbq_gagg(0:2),qg_gaqg(0:2),gq_gaqg(0:2),
     & qbg_gaqbg(0:2),gqb_gaqbg(0:2),gg_gaqbq(0:2),gg_gaqqb(0:2),tmp
      real(dp):: msqv_cs(0:2,-nf:nf,-nf:nf)
      real(dp):: ppmsqvx(ppmax),tmpmsqvx(0:ppmax)
      
      msqv(:,:)=zip
      msqv_cs(:,:,:)=zip
      tmpmsqvx=zip

      qqb_gagg(:) =zip
      qbq_gagg(:) =zip
      gq_gaqg(:)  =zip
      qg_gaqg(:)  =zip
      gqb_gaqbg(:)=zip
      qbg_gaqbg(:)=zip
      gg_gaqbq(:) =zip
      gg_gaqqb(:) =zip

      call spinoru(5,p,za,zb)
C   zab=<i-|k|j-> zba=<i+|k|j+> where k is an arbitrary 4-vector 
      call spinork(5,p,zab,zba,n)

c--- check validity of vector n
      call checkndotp(p,n,in)

! DEBUG
!      n(:)=zip
!      n(1)=one
!      if (in > 2) then
!        n(1)=p(in,2)/sqrt(p(in,1)**2+p(in,2)**2)
!        n(2)=-p(in,1)/sqrt(p(in,1)**2+p(in,2)**2)
!        n(3)=0._dp
!        n(4)=0._dp       
!      endif
!      call spinork(5,p,zab,zba,n)
!c--- check validity of vector n
!      call checkndotp(p,n,in)
!      call namp_qqbggga(5,1,2,4,3,za,zb,zab,zba,amp1)
!      call namp_qqbggga(1,5,2,4,3,zb,za,zba,zab,amp1s)
!      n(:)=zip
!      n(2)=one
!      if (in > 2) then
!        n(1)=p(in,1)*p(in,3)/p(in,4)/sqrt(p(in,1)**2+p(in,2)**2)
!        n(2)=p(in,2)*p(in,3)/p(in,4)/sqrt(p(in,1)**2+p(in,2)**2)
!        n(3)=-sqrt(p(in,1)**2+p(in,2)**2)/p(in,4)
!        n(4)=0._dp       
!      endif
!      call spinork(5,p,zab,zba,n)
!c--- check validity of vector n
!      call checkndotp(p,n,in)
!      call namp_qqbggga(5,1,2,4,3,za,zb,zab,zba,amp2)
!      call namp_qqbggga(1,5,2,4,3,zb,za,zba,zab,amp2s)
!      
!      call qqbgg_ga_amp(5,1,2,4,3,za,zb,amplo)
!      call qqbgg_ga_amp(5,1,4,2,3,za,zb,amplos)
!      
!      do h1=1,2
!      do h2=1,2
!      do h3=1,2
!      write(6,*) h1,h2,h3,abs(amp1(h1,h2,h3)+amp1s(h1,3-h2,3-h3))**2+abs(amp2(h1,h2,h3)+amp2s(h1,3-h2,3-h3))**2,
!     & (abs(amplo(h1,h2,1,h3)+amplos(h1,1,h2,h3))**2+abs(amplo(h1,h2,2,h3)+amplos(h1,2,h2,h3))**2)
!      enddo
!      enddo
!      enddo
!      pause
! DEBUG
      
      if (in == 1) then
Cargument 1-4 represent (1) incoming quark line
C                       (2) incoming anti-quark line
C                       (3) outgoing photon
C                       (4) outgoing gluon 
C                       (5) outgoing gluon contracted with n
        call gam2jsqn(4,5,3,2,1,za,zb,zab,zba,gg_gaqbq)
        call gam2jsqn(2,4,3,5,1,za,zb,zab,zba,gq_gaqg)
        call gam2jsqn(4,2,3,5,1,za,zb,zab,zba,gqb_gaqbg)
c--- extra term for gg -> qqb
        call gam2jsqn(5,4,3,2,1,za,zb,zab,zba,gg_gaqqb)
      elseif (in == 2)  then
        call gam2jsqn(4,5,3,1,2,za,zb,zab,zba,gg_gaqbq)
        call gam2jsqn(4,1,3,5,2,za,zb,zab,zba,qbg_gaqbg)
        call gam2jsqn(1,4,3,5,2,za,zb,zab,zba,qg_gaqg)
c--- extra term for gg -> qqb
        call gam2jsqn(5,4,3,1,2,za,zb,zab,zba,gg_gaqqb)
c--- interchange color structures for gg contributions since 1 and 2 exchanged
        tmp=gg_gaqbq(1)
        gg_gaqbq(1)=gg_gaqbq(2)
        gg_gaqbq(2)=tmp
        tmp=gg_gaqqb(1)
        gg_gaqqb(1)=gg_gaqqb(2)
        gg_gaqqb(2)=tmp
      elseif (in == 4) then
        call gam2jsqn(1,2,3,5,4,za,zb,zab,zba,qqb_gagg)
        call gam2jsqn(2,1,3,5,4,za,zb,zab,zba,qbq_gagg)
        call gam2jsqn(1,5,3,2,4,za,zb,zab,zba,qg_gaqg)
        call gam2jsqn(2,5,3,1,4,za,zb,zab,zba,gq_gaqg)
        call gam2jsqn(5,1,3,2,4,za,zb,zab,zba,qbg_gaqbg)
        call gam2jsqn(5,2,3,1,4,za,zb,zab,zba,gqb_gaqbg)
      elseif (in == 5) then
        call gam2jsqn(1,2,3,4,5,za,zb,zab,zba,qqb_gagg)
        call gam2jsqn(2,1,3,4,5,za,zb,zab,zba,qbq_gagg)
        call gam2jsqn(1,4,3,2,5,za,zb,zab,zba,qg_gaqg)
        call gam2jsqn(2,4,3,1,5,za,zb,zab,zba,gq_gaqg)
        call gam2jsqn(4,1,3,2,5,za,zb,zab,zba,qbg_gaqbg)
        call gam2jsqn(4,2,3,1,5,za,zb,zab,zba,gqb_gaqbg)
      else
        write(6,*) 'Invalid value of in in qqb_gam2j_gvecx_new.f:  in=',in
        stop
      endif
      
      fac=two*esq*gsq**2*V*xn
c--- apply overall factors
      qqb_gagg(:) =half*aveqq*fac*qqb_gagg(:)
      qbq_gagg(:) =half*aveqq*fac*qbq_gagg(:)
      gq_gaqg(:)  =aveqg*fac*gq_gaqg(:)
      qg_gaqg(:)  =aveqg*fac*qg_gaqg(:)
      gqb_gaqbg(:)=aveqg*fac*gqb_gaqbg(:)
      qbg_gaqbg(:)=aveqg*fac*qbg_gaqbg(:)
      gg_gaqbq(:) =avegg*fac*gg_gaqbq(:)
      gg_gaqqb(:) =avegg*fac*gg_gaqqb(:)

      do j=-nf,nf
      do k=-nf,nf
      
      if((j .ne. 0) .and. (k .ne. 0) .and. (j .ne. -k)) cycle

      if     ((j == 0) .and. (k == 0)) then
        do icol=0,2
          do nquark=1,2
          if (icol == 0) then
            tmpmsqvx(pp(j,k,-nquark,nquark))=zip
            tmpmsqvx(pp(j,k,nquark,-nquark))=zip
          endif
          tmpmsqvx(pp(j,k,-nquark,nquark))=tmpmsqvx(pp(j,k,-nquark,nquark))
     &     +Q(nquark)**2*gg_gaqbq(icol)
          tmpmsqvx(pp(j,k,nquark,-nquark))=tmpmsqvx(pp(j,k,nquark,-nquark))
     &     +Q(nquark)**2*gg_gaqqb(icol)
          enddo
          msqv_cs(icol,j,k)=(3._dp*Q(1)**2+2._dp*Q(2)**2)*gg_gaqbq(icol)
        enddo
      elseif ((j > 0) .and. (k < 0)) then
        msqv_cs(:,j,k)=Q(j)**2*qqb_gagg(:)
      elseif ((j < 0) .and. (k > 0)) then
        msqv_cs(:,j,k)=Q(k)**2*qbq_gagg(:)
      elseif ((j > 0) .and. (k == 0)) then
c--- normal case
        msqv_cs(:,j,k)=Q(j)**2*qg_gaqg(:)
             if (j <= 4) then
             tmpmsqvx(pp(j,k,j,k))=msqv_cs(0,j,k)+msqv_cs(1,j,k)+msqv_cs(2,j,k)
             tmpmsqvx(pp(j,k,k,j))=tmpmsqvx(pp(j,k,j,k))
          endif
      elseif ((j < 0) .and. (k == 0)) then
        msqv_cs(:,j,k)=Q(-j)**2*qbg_gaqbg(:)
             if (j >= -4) then
             tmpmsqvx(pp(j,k,j,k))=msqv_cs(0,j,k)+msqv_cs(1,j,k)+msqv_cs(2,j,k)
             tmpmsqvx(pp(j,k,k,j))=tmpmsqvx(pp(j,k,j,k))
             endif
      elseif ((j == 0) .and. (k > 0)) then
        msqv_cs(:,j,k)=Q(k)**2*gq_gaqg(:)
             if (k <= 4) then
             tmpmsqvx(pp(j,k,k,j))=msqv_cs(0,j,k)+msqv_cs(1,j,k)+msqv_cs(2,j,k)
             tmpmsqvx(pp(j,k,j,k))=tmpmsqvx(pp(j,k,k,j))
             endif
      elseif ((j == 0) .and. (k < 0)) then
        msqv_cs(:,j,k)=Q(-k)**2*gqb_gaqbg(:)
             if (k >= -4) then
             tmpmsqvx(pp(j,k,k,j))=msqv_cs(0,j,k)+msqv_cs(1,j,k)+msqv_cs(2,j,k)
             tmpmsqvx(pp(j,k,j,k))=tmpmsqvx(pp(j,k,k,j))
             endif
      endif
      msqv(j,k)=msqv_cs(0,j,k)+msqv_cs(1,j,k)+msqv_cs(2,j,k)
      enddo
      enddo

c--- fill ppmsqx array to be returned
      ppmsqvx(1:ppmax)=tmpmsqvx(1:ppmax)

      return
      end



      subroutine gam2jsqn(i1,i2,i3,i4,i5,za,zb,zab,zba,msqn_cs)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      integer i1,i2,i3,i4,i5,h1,h2,h3
      real(dp):: msqn_cs(0:2)
      complex(dp):: amp34(2,2,2),amp43(2,2,2),ampQED,
     & zab(mxpart,mxpart),zba(mxpart,mxpart),swap

c--- note that ordering is permuted wrt namp_qqbgga
      call namp_qqbggga(i2,i1,i4,i5,i3,za,zb,zab,zba,amp34)
c--- obtain other ordering by c.c. and relabelling
      call namp_qqbggga(i1,i2,i4,i5,i3,zb,za,zba,zab,amp43)
      do h1=1,2
      swap=amp43(h1,1,1)
      amp43(h1,1,1)=amp43(h1,2,2)
      amp43(h1,2,2)=swap
      swap=amp43(h1,1,2)
      amp43(h1,1,2)=amp43(h1,2,1)
      amp43(h1,2,1)=swap
      enddo

!----- sum over polarizations
      msqn_cs(:)=zip
      do h1=1,2
      do h2=1,2
      do h3=1,2
        ampQED=amp34(h1,h2,h3)+amp43(h1,h2,h3)
        msqn_cs(1)=msqn_cs(1)+real(amp34(h1,h2,h3)*conjg(amp34(h1,h2,h3)),dp)
        msqn_cs(2)=msqn_cs(2)+real(amp43(h1,h2,h3)*conjg(amp43(h1,h2,h3)),dp)
        msqn_cs(0)=msqn_cs(0)+real(ampQED*conjg(ampQED),dp)
      enddo
      enddo
      enddo

      msqn_cs(0)=-msqn_cs(0)/xnsq
      
!      write(6,*) 'gvec: amp 34 bit',msqn_cs(1)
!      write(6,*) 'gvec: amp 43 bit',msqn_cs(2)
!      write(6,*) 'gvec: QED bit',msqn_cs(0)

      return
      end
      

