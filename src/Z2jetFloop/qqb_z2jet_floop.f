      subroutine qqb_z2jet_floop(p,msq)
      implicit none
      include 'types.f'
*                                                                      *
*    Fermion loops only                                                *
*                                                                      *
*                                                                      *
*                                                                      *
      
************************************************************************
*     Author: R.K. Ellis                                               *
*     September 2001.                                                  *
*                                                                      *
*     Calculate the virtual matrix element squared and                 *
*     subtraction terms for the process                                *
*                                                                      *
*     q(-p1) + qbar(-p2) --> Z + j(p5) + j(p6)                         *
*                            |                                         *
*                            --> e^-(p3) + e^+(p4)                     *
*                                                                      *
*     where the partons are either q(p5) and qbar(p6) [Qflag = .true.] *
*                               or g(p5) and g(p6)    [Gflag = .true.] *
*                                                                      *
*     The value of COLOURCHOICE determines which colour structures     *
*     are included in the terms for the QQGG piece                     *
*                                                                      *
************************************************************************

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'ewcharge.f'
      include 'zcouple.f'
      include 'masses.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      include 'epinv.f'
      include 'scheme.f'
      include 'flags.f'
      include 'msq_cs.f'
      include 'lc.f'
      include 'ppmax.f'
      include 'mpicommon.f'
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),pswap(mxpart,4),fac
      complex(dp):: mmsq_qqb_ax(2,2),mmsq_qbq_ax(2,2),
     & mmsq_gq_ax(2,2),mmsq_qg_ax(2,2),mmsq_qbg_ax(2,2),
     & mmsq_gqb_ax(2,2),mmsq_gg_ax(2,2),
     & mmsq_qqb_vec0(2,2),mmsq_qbq_vec0(2,2),
     & mmsq_gq_vec0(2,2),mmsq_qg_vec0(2,2),mmsq_qbg_vec0(2,2),
     & mmsq_gqb_vec0(2,2),mmsq_gg_vec0(2,2),
     & mmsq_qqb_vect(2,2),mmsq_qbq_vect(2,2),
     & mmsq_gq_vect(2,2),mmsq_qg_vect(2,2),mmsq_qbg_vect(2,2),
     & mmsq_gqb_vect(2,2),mmsq_gg_vect(2,2)
      complex(dp):: prop,vcouple(2),tcouple(2)
      integer:: nu,j,k,polq,polz,nup,ndo
      integer, parameter :: jj(-nf:nf)=(/-1,-2,-1,-2,-1,0,1,2,1,2,1/)
      real(dp):: v2(2),vQ(nf,2)

      parameter (nup=2,ndo=nf-nup)
! discardvecax: if true, discard diagrams where Z couples directly
! to a closed quark loop, through vector and axial-vector couplings
      logical, parameter:: discardvecax=.false., madloopcheck=.false.
      include 'cplx.h'

      msq(:,:)=zip

      prop=s(3,4)/cplx2(s(3,4)-zmass**2,zmass*zwidth)
      
      v2(1)=l1
      v2(2)=r1

      do j=1,nf
        vQ(j,1)=L(j)
        vQ(j,2)=R(j)
      enddo

      if (madloopcheck) then
c--- comparison with MadLoop, arXiv: 1103.0621, section A.2.6
        call ps_check(p,5)
c--- set up couplings
        L(1)=-one
        R(1)=+one
        L(2)=+one
        R(2)=-one
        L(3)=L(1)
        R(3)=R(1)
        L(4)=L(2)
        R(4)=R(2)
        L(5)=L(1)
        R(5)=R(1)
        v2(1)=one
        v2(2)=one
        do j=1,nf
          vQ(j,1)=L(j)
          vQ(j,2)=R(j)
        enddo
c--- set up crossing for e+ e- -> d d~ g g
        do nu=1,4
        pswap(1,nu)=p(3,nu)
        pswap(2,nu)=p(5,nu)
        pswap(3,nu)=p(6,nu)
        pswap(4,nu)=p(4,nu)
        pswap(5,nu)=p(2,nu)
        pswap(6,nu)=p(1,nu)
        enddo
        call spinoru(6,pswap,za,zb)
        prop=s(5,6)/cplx2(s(5,6)-zmass**2,zmass*zwidth)
        call xzqqgg_floop(pswap,mmsq_qqb_vec0,mmsq_qqb_vect,mmsq_qqb_ax)
        mmsq_qqb_vec0=mmsq_qqb_vec0/esq**2*(quarter/avegg) ! correct averaging, remove couplings
        mmsq_qqb_vect=mmsq_qqb_vect/esq**2*(quarter/avegg) ! correct averaging, remove couplings
        mmsq_qqb_ax  =mmsq_qqb_ax  /esq**2*(quarter/avegg) ! correct averaging, remove couplings
        fac=zip
        do polq=1,2
        do polz=1,2
! remember identical factor of half for two gluons in final state, no intermediate Z
         fac=fac+half*real(mmsq_qqb_vec0(polq,polz)
     &          *conjg(1d0+vQ(j,polq)*v2(polz)*prop*0)*1d0)
        enddo
        enddo
        write(6,*) 'vector mq=0',fac
        fac=zip
        do polq=1,2
        do polz=1,2
! remember identical factor of half for two gluons in final state
! no intermediate photon or vector Z, tree-level coupling is to down quarks
         fac=fac+half*real(mmsq_qqb_ax(polq,polz)
     &          *conjg(Q(j)*q1*0+vQ(1,polq)*v2(polz)*prop)
     &            *(two*v2(polz)*prop))
        enddo
        enddo
        write(6,*) 'axial',fac
        pause
      endif

c--- Now calculate the relevant lowest-order matrix elements
c--- for each possible initial state from the QQGG contribution
      
c---  calculate the qqb terms
c-BDKW  0--> q(p1)+g(p2)+g(p3)+qb(p4)+lbar(p5)+l(p6)
CALL    0--> q(p2)+g(p6)+g(p5)+qb(p1)+lbar(p4)+l(p3)
      do nu=1,4
      pswap(1,nu)=p(2,nu)
      pswap(2,nu)=p(6,nu)
      pswap(3,nu)=p(5,nu)
      pswap(4,nu)=p(1,nu)
      pswap(5,nu)=p(4,nu)
      pswap(6,nu)=p(3,nu)
      enddo
      call spinoru(6,pswap,za,zb)
      call xzqqgg_floop(pswap,mmsq_qqb_vec0,mmsq_qqb_vect,mmsq_qqb_ax)

c--- obtain qbq from qqb by symmetry
      do polq=1,2
      do polz=1,2
        mmsq_qbq_vec0(polq,polz)=mmsq_qqb_vec0(3-polq,polz)
        mmsq_qbq_vect(polq,polz)=mmsq_qqb_vect(3-polq,polz)
        mmsq_qbq_ax(polq,polz)=-mmsq_qqb_ax(3-polq,polz)
      enddo
      enddo     
c---  calculate the qbq terms
c-BDKW  0--> q(p1)+g(p2)+g(p3)+qb(p4)+lbar(p5)+l(p6)
CALL    0--> q(p1)+g(p6)+g(p5)+qb(p2)+lbar(p4)+l(p3)
c      do nu=1,4
c      pswap(1,nu)=p(1,nu)
c      pswap(2,nu)=p(6,nu)
c      pswap(3,nu)=p(5,nu)
c      pswap(4,nu)=p(2,nu)
c      pswap(5,nu)=p(4,nu)
c      pswap(6,nu)=p(3,nu)
c      enddo
c      call spinoru(6,pswap,za,zb)
c      call xzqqgg_floop(mmsq_qbq,mmsq_qbq_vec,mmsq_qbq_ax)
c      do polq=1,2
c      do polz=1,2
c        write(6,*) mmsq_qbq(polq,polz)-mmsq_qqb(3-polq,polz)
c        write(6,*) mmsq_qbq_vec0(polq,polz)-mmsq_qqb_vec0(3-polq,polz)
c        write(6,*) mmsq_qbq_ax(polq,polz)+mmsq_qqb_ax(3-polq,polz)
c      enddo
c      enddo
c      pause

c---  calculate the gq terms
c-BDKW  0--> q(p1)+g(p2)+g(p3)+qb(p4)+lbar(p5)+l(p6)
CALL    0--> q(p5)+g(p6)+g(p1)+qb(p2)+lbar(p4)+l(p3)
      do nu=1,4
      pswap(1,nu)=p(5,nu)
      pswap(2,nu)=p(6,nu)
      pswap(3,nu)=p(1,nu)
      pswap(4,nu)=p(2,nu)
      pswap(5,nu)=p(4,nu)
      pswap(6,nu)=p(3,nu)
      enddo
      call spinoru(6,pswap,za,zb)
      call xzqqgg_floop(pswap,mmsq_gq_vec0,mmsq_gq_vect,mmsq_gq_ax)

c--- obtain gqb from gq by symmetry
      do polq=1,2
      do polz=1,2
        mmsq_gqb_vec0(polq,polz)=mmsq_gq_vec0(3-polq,polz)
        mmsq_gqb_vect(polq,polz)=mmsq_gq_vect(3-polq,polz)
        mmsq_gqb_ax(polq,polz)=-mmsq_gq_ax(3-polq,polz)
      enddo
      enddo
c---  calculate the gqb terms
c-BDKW  0--> q(p1)+g(p2)+g(p3)+qb(p4)+lbar(p5)+l(p6)
CALL    0--> q(p2)+g(p6)+g(p1)+qb(p5)+lbar(p4)+l(p3)
c      do nu=1,4
c      pswap(1,nu)=p(2,nu)
c      pswap(2,nu)=p(6,nu)
c      pswap(3,nu)=p(1,nu)
c      pswap(4,nu)=p(5,nu)
c      pswap(5,nu)=p(4,nu)
c      pswap(6,nu)=p(3,nu)
c      enddo
c      call spinoru(6,pswap,za,zb)
c      call xzqqgg_floop(mmsq_gqb,mmsq_gqb_vec,mmsq_gqb_ax)

c---  calculate the qg terms
c-BDKW  0--> q(p1)+g(p2)+g(p3)+qb(p4)+lbar(p5)+l(p6)
CALL    0--> q(p5)+g(p6)+g(p2)+qb(p1)+lbar(p4)+l(p3)
      do nu=1,4
      pswap(1,nu)=p(5,nu)
      pswap(2,nu)=p(6,nu)
      pswap(3,nu)=p(2,nu)
      pswap(4,nu)=p(1,nu)
      pswap(5,nu)=p(4,nu)
      pswap(6,nu)=p(3,nu)
      enddo
      call spinoru(6,pswap,za,zb)
      call xzqqgg_floop(pswap,mmsq_qg_vec0,mmsq_qg_vect,mmsq_qg_ax)

c--- obtain qbg from qg by symmetry
      do polq=1,2
      do polz=1,2
        mmsq_qbg_vec0(polq,polz)=mmsq_qg_vec0(3-polq,polz)
        mmsq_qbg_vect(polq,polz)=mmsq_qg_vect(3-polq,polz)
        mmsq_qbg_ax(polq,polz)=-mmsq_qg_ax(3-polq,polz)
      enddo
      enddo
c---  calculate the qbg terms
c-BDKW  0--> q(p1)+g(p2)+g(p3)+qb(p4)+lbar(p5)+l(p6)
CALL    0--> q(p1)+g(p6)+g(p2)+qb(p5)+lbar(p4)+l(p3)
c      do nu=1,4
c      pswap(1,nu)=p(1,nu)
c      pswap(2,nu)=p(6,nu)
c      pswap(3,nu)=p(2,nu)
c      pswap(4,nu)=p(5,nu)
c      pswap(5,nu)=p(4,nu)
c      pswap(6,nu)=p(3,nu)
c      enddo
c      call spinoru(6,pswap,za,zb)
c      call xzqqgg_floop(mmsq_qbg,mmsq_qbg_vec,mmsq_qbg_ax)

c--- calculate the gg terms
c-BDKW  0--> q(p1)+g(p2)+g(p3)+qb(p4)+lbar(p5)+l(p6)
CALL    0--> q(p6)+g(p1)+g(p2)+qb(p5)+lbar(p4)+l(p3)
      do nu=1,4
      pswap(1,nu)=p(6,nu)
      pswap(2,nu)=p(1,nu)
      pswap(3,nu)=p(2,nu)
      pswap(4,nu)=p(5,nu)
      pswap(5,nu)=p(4,nu)
      pswap(6,nu)=p(3,nu)
      enddo
      call spinoru(6,pswap,za,zb)
      call xzqqgg_floop(pswap,mmsq_gg_vec0,mmsq_gg_vect,mmsq_gg_ax)
      
      if (discardvecax) then
        mmsq_qqb_ax(:,:)=zip
        mmsq_qbq_ax(:,:)=zip
        mmsq_gq_ax(:,:)=zip
        mmsq_qg_ax(:,:)=zip
        mmsq_qbg_ax(:,:)=zip
        mmsq_gqb_ax(:,:)=zip
        mmsq_gg_ax(:,:)=zip
        mmsq_qqb_vec0(:,:)=zip
        mmsq_qbq_vec0(:,:)=zip
        mmsq_gq_vec0(:,:)=zip
        mmsq_qg_vec0(:,:)=zip
        mmsq_qbg_vec0(:,:)=zip
        mmsq_gqb_vec0(:,:)=zip
        mmsq_gg_vec0(:,:)=zip
        mmsq_qqb_vect(:,:)=zip
        mmsq_qbq_vect(:,:)=zip
        mmsq_gq_vect(:,:)=zip
        mmsq_qg_vect(:,:)=zip
        mmsq_qbg_vect(:,:)=zip
        mmsq_gqb_vect(:,:)=zip
        mmsq_gg_vect(:,:)=zip
      endif

c--- compute correct vector-like coupling for diagrams with Z coupled to a loop
      vcouple(1)=czip
      vcouple(2)=czip
      do polz=1,2
      do j=1,nf
      vcouple(polz)=vcouple(polz)
     & +Q(j)*q1+0.5_dp*(vQ(j,1)+vQ(j,2))*v2(polz)*prop
      enddo
      tcouple(polz)=Q(2)*q1+0.5_dp*(vQ(2,1)+vQ(2,2))*v2(polz)*prop
      enddo
      
      do j=-nf,nf
      do k=-nf,nf
      do polq=1,2
      do polz=1,2

c--- quark-antiquark
      if ((j > 0) .and. (k < 0)) then
        if (j == -k) 
     &  msq(j,k)=msq(j,k)+half*(aveqq/avegg)*(
     &                 +real(mmsq_qqb_vec0(polq,polz)
     &        *conjg(Q(j)*q1+vQ(j,polq)*v2(polz)*prop)*vcouple(polz))
     &                 +real(mmsq_qqb_vect(polq,polz)
     &        *conjg(Q(j)*q1+vQ(j,polq)*v2(polz)*prop)*tcouple(polz))
     &                 +real(mmsq_qqb_ax(polq,polz)
     &        *conjg(Q(j)*q1+vQ(j,polq)*v2(polz)*prop)
     &          *(v2(polz)*prop)/sin2w))

c--- antiquark-quark
      elseif ((j < 0) .and. (k > 0)) then
        if (j == -k)
     &  msq(j,k)=msq(j,k)+half*(aveqq/avegg)*(
     &                 +real(mmsq_qbq_vec0(polq,polz)
     &        *conjg(Q(k)*q1+vQ(k,polq)*v2(polz)*prop)*vcouple(polz))
     &                 +real(mmsq_qbq_vect(polq,polz)
     &        *conjg(Q(k)*q1+vQ(k,polq)*v2(polz)*prop)*tcouple(polz))
     &                 +real(mmsq_qbq_ax(polq,polz)
     &        *conjg(Q(k)*q1+vQ(k,polq)*v2(polz)*prop)
     &          *(v2(polz)*prop)/sin2w))


c--- quark-gluon
      elseif ((j > 0) .and. (k == 0)) then
            
      msq(j,k)=msq(j,k)+(aveqg/avegg)*(
     &                 +real(mmsq_qg_vec0(polq,polz)
     &        *conjg(Q(j)*q1+vQ(j,polq)*v2(polz)*prop)*vcouple(polz))
     &                 +real(mmsq_qg_vect(polq,polz)
     &        *conjg(Q(j)*q1+vQ(j,polq)*v2(polz)*prop)*tcouple(polz))
     &                 +real(mmsq_qg_ax(polq,polz)
     &        *conjg(Q(j)*q1+vQ(j,polq)*v2(polz)*prop)
     &          *(v2(polz)*prop)/sin2w))

c--- antiquark-gluon
      elseif ((j < 0) .and. (k == 0)) then

        msq(j,k)=msq(j,k)+(aveqg/avegg)*(
     &                 +real(mmsq_qbg_vec0(polq,polz)
     &        *conjg(Q(-j)*q1+vQ(-j,polq)*v2(polz)*prop)*vcouple(polz))
     &                 +real(mmsq_qbg_vect(polq,polz)
     &        *conjg(Q(-j)*q1+vQ(-j,polq)*v2(polz)*prop)*tcouple(polz))
     &                 +real(mmsq_qbg_ax(polq,polz)
     &        *conjg(Q(-j)*q1+vQ(-j,polq)*v2(polz)*prop)
     &          *(v2(polz)*prop)/sin2w))

c--- gluon-quark
      elseif ((j == 0) .and. (k > 0)) then
        msq(j,k)=msq(j,k)+(aveqg/avegg)*(
     &                 +real(mmsq_gq_vec0(polq,polz)
     &        *conjg(Q(k)*q1+vQ(k,polq)*v2(polz)*prop)*vcouple(polz))
     &                 +real(mmsq_gq_vect(polq,polz)
     &        *conjg(Q(k)*q1+vQ(k,polq)*v2(polz)*prop)*tcouple(polz))
     &                 +real(mmsq_gq_ax(polq,polz)
     &        *conjg(Q(k)*q1+vQ(k,polq)*v2(polz)*prop)
     &          *(v2(polz)*prop)/sin2w))

c--- gluon-antiquark
      elseif ((j == 0) .and. (k < 0)) then
        msq(j,k)=msq(j,k)+(aveqg/avegg)*(
     &                 +real(mmsq_gqb_vec0(polq,polz)
     &        *conjg(Q(-k)*q1+vQ(-k,polq)*v2(polz)*prop)*vcouple(polz))
     &                 +real(mmsq_gqb_vect(polq,polz)
     &        *conjg(Q(-k)*q1+vQ(-k,polq)*v2(polz)*prop)*tcouple(polz))
     &                 +real(mmsq_gqb_ax(polq,polz)
     &        *conjg(Q(-k)*q1+vQ(-k,polq)*v2(polz)*prop)
     &          *(v2(polz)*prop)/sin2w))

c--- gluon-gluon
      elseif ((j == 0) .and. (k == 0)) then
        msq(j,k)=msq(j,k)+real(ndo,dp)*(
     &                 +real(mmsq_gg_vec0(polq,polz)
     &        *conjg(Q(1)*q1+vQ(1,polq)*v2(polz)*prop)*vcouple(polz))
     &                 +real(mmsq_gg_vect(polq,polz)
     &        *conjg(Q(1)*q1+vQ(1,polq)*v2(polz)*prop)*tcouple(polz))
     &                 +real(mmsq_gg_ax(polq,polz)
     &        *conjg(Q(1)*q1+vQ(1,polq)*v2(polz)*prop)
     &          *(v2(polz)*prop)/sin2w))
        msq(j,k)=msq(j,k)+real(nup,dp)*(
     &                 +real(mmsq_gg_vec0(polq,polz)
     &        *conjg(Q(2)*q1+vQ(2,polq)*v2(polz)*prop)*vcouple(polz))
     &                 +real(mmsq_gg_vect(polq,polz)
     &        *conjg(Q(2)*q1+vQ(2,polq)*v2(polz)*prop)*tcouple(polz))
     &                 +real(mmsq_gg_ax(polq,polz)
     &        *conjg(Q(2)*q1+vQ(2,polq)*v2(polz)*prop)
     &          *(v2(polz)*prop)/sin2w))
      endif
      enddo
      enddo
      enddo
      enddo

      return
      end
