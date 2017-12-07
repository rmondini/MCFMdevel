      subroutine xzqqgg_floop(p,mqqb_vec0,mqqb_vect,mqqb_ax)
      implicit none
      include 'types.f'
      
************************************************************************
*     Author J.M.Campbell, February 2000                               *
*     Returns the interference of the tree and loop                    *
*     amplitudes for the process                                       *
*     0---> q(p1)+g(p2)+g(p3)+qbar(p4)+l(p5)+a(p6)                     *
*                                                                      *
*     The value of COLOURCHOICE determines which colour structures     *
*     are included in the calculation                                  *
*                                                                      *
*     mqqb(2,2) has two indices - the first for the helicity of the    *
*     quark line, the second for the helicity of the lepton line.      *
************************************************************************
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'zprods_com.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'lc.f'
      include 'masses.f'
      integer:: j,lh,h2,h3,hq,h(2:3)
      real(dp):: p(mxpart,4),pfloop(mxpart,4),fac
      complex(dp):: m(2),ml1(2),ml2(2),ml3,ml4(2),
     & mqqb_vec0(2,2),mqqb_vect(2,2),mqqb_ax(2,2)
      complex(dp):: ml_vec0(2),ml_vect(2),ml1_ax(2),ml2_ax(2)
      complex(dp):: a6treeg1,
     & a61g1lc,a61g1slc,a61g1nf,a63g1,a64v,a64ax,a65ax,
     & ampTree(2,2,2,2,2),
     & ampVec0(2,2,2,2,2),ampAx0(0:2,2,2,2,2),
     & ampVect(2,2,2,2,2),ampAxt(0:2,2,2,2,2)
      integer,parameter::i1(2)=(/1,4/),i2(2)=(/2,3/),i3(2)=(/3,2/),
     &                   i4(2)=(/4,1/),i5(2)=(/6,5/),i6(2)=(/5,6/)
      character*9,parameter:: st1(2,2)=
     & reshape((/'q+g-g-qb-','q+g-g+qb-','q+g+g-qb-','q+g+g+qb-'/)
     & ,(/2,2/))
      character*9,parameter:: st2(2,2)=
     & reshape((/'q+qb-g+g+','q+qb-g+g-','q+qb-g-g+','q+qb-g-g-'/)
     & ,(/2,2/))
      character*9,parameter:: st3(2,2)=
     & reshape((/'q+qb-g-g-','q+qb-g-g+','q+qb-g+g-','q+qb-g+g+'/)
     & ,(/2,2/))
      logical, parameter:: exact=.true.
      include 'cplx.h'

C ---final matrix element squared is needed as function of quark line helicity
C----and lepton line helicity
C----first argument is quark line helicity
C----second argument is lepton line helicity

      fac=avegg*8._dp*gsq**2*esq**2*cf*xn**3*ason2pi
c--- no extra factor here since colour algebra is already done in (2.12)

      if (exact) then
c--- calls to amplitudes with exact mt-dependence
        pfloop(1,:)=p(1,:)
        pfloop(2,:)=p(4,:)
        pfloop(3,:)=p(3,:)
        pfloop(4,:)=p(2,:)
        pfloop(5,:)=p(5,:)
        pfloop(6,:)=p(6,:)
!        call qqbZgg_floop(pfloop,zip,ampTree,ampVec0,ampAx0)
!        call qqbZgg_floop(pfloop,mt,ampTree,ampVect,ampAxt)
!        ampVect=zip
!        ampAxt=zip
! need to reset spinor products after qqbZgg_floop has spoiled them
        call spinoru(6,p,za,zb)
      endif
      
      do hq=1,2
      do lh=1,2
      mqqb_vec0(hq,lh)=czip
      mqqb_vect(hq,lh)=czip
      mqqb_ax(hq,lh)=czip

      do h2=1,2
      do h3=1,2
c--- calculate all amplitudes
!        if (exact .eqv. .false.) then
          h(2)=h2
          h(3)=h3
          do j=1,2
          if (hq == 1) then
          m(j)=  a6treeg1(st1(3-h(i2(j)),3-h(i3(j))),
     &       i1(1),i2(j),i3(j),i4(1),i6(lh),i5(lh),zb,za)
          ml_vec0(j)=a64v(st3(3-h(i2(j)),3-h(i3(j))),
     &       i1(1),i4(1),i2(j),i3(j),i6(lh),i5(lh),zb,za)
          ml_vect(j)=zip
c--- note: this symmetry relation (including minus sign) checked numerically 
          ml1_ax(j)=-a64ax(st3(3-h(i2(j)),3-h(i3(j))),
     &       i1(1),i4(1),i2(j),i3(j),i6(lh),i5(lh),zb,za)
          ml2_ax(j)=-a65ax(st3(3-h(i2(j)),3-h(i3(j))),
     &       i1(1),i4(1),i2(j),i3(j),i6(lh),i5(lh),zb,za)
          else
          m(j)=  a6treeg1(st1(h(i2(j)),h(i3(j))),
     &       i1(1),i2(j),i3(j),i4(1),i5(lh),i6(lh),za,zb)
          ml_vec0(j)=a64v(st3(h(i2(j)),h(i3(j))),
     &       i1(1),i4(1),i2(j),i3(j),i5(lh),i6(lh),za,zb)
          ml_vect(j)=zip
          ml1_ax(j)=a64ax(st3(h(i2(j)),h(i3(j))),
     &       i1(1),i4(1),i2(j),i3(j),i5(lh),i6(lh),za,zb)
          ml2_ax(j)=a65ax(st3(h(i2(j)),h(i3(j))),
     &       i1(1),i4(1),i2(j),i3(j),i5(lh),i6(lh),za,zb)
          endif
          enddo
!        endif

c--- use amplitudes with exact results, if required
c--- note that ml_vec[0/t](1)=ml_vec[0/t](2), ml2_ax(1)=ml2_ax(2)
c--- Equivalence between BDK amplitudes and explicit new calculation  
!      write(6,*) 'vc1',hq,lh,h2,h3,ml_vec0(1)/ampVec0(1,hq,h2,h3,lh)
!      write(6,*) 'vc2',hq,lh,h2,h3,ml_vec0(2)/ampVec0(2,hq,h2,h3,lh)
!      write(6,*) 'ax0',hq,lh,h2,h3,ml2_ax(1)/((ampAxt(0,hq,h2,h3,lh)-Ampax0(0,hq,h2,h3,lh))/two)
!      write(6,*) 'ax1',hq,lh,h2,h3,ml1_ax(1)/((ampAxt(2,hq,h2,h3,lh)-Ampax0(2,hq,h2,h3,lh))/two)
!      write(6,*) 'ax2',hq,lh,h2,h3,ml1_ax(2)/((ampAxt(1,hq,h2,h3,lh)-Ampax0(1,hq,h2,h3,lh))/two)
!      write(6,*) 'lo1',hq,lh,h2,h3,m(1)/ampTree(2,hq,h2,h3,lh)
!      write(6,*) 'lo2',hq,lh,h2,h3,m(2)/ampTree(1,hq,h2,h3,lh)

      if (exact) then
        m(1)=ampTree(2,hq,h2,h3,lh)
        m(2)=ampTree(1,hq,h2,h3,lh)
        ml_vec0(1)=ampVec0(1,hq,h2,h3,lh)
        ml_vec0(2)=ampVec0(2,hq,h2,h3,lh)
        ml_vect(1)=ampVect(1,hq,h2,h3,lh)
        ml_vect(2)=ampVect(2,hq,h2,h3,lh)
        ml2_ax(1)=(ampAxt(0,hq,h2,h3,lh)-Ampax0(0,hq,h2,h3,lh))/two
        ml2_ax(2)=ml2_ax(1)
        ml1_ax(1)=(ampAxt(2,hq,h2,h3,lh)-Ampax0(2,hq,h2,h3,lh))/two
        ml1_ax(2)=(ampAxt(1,hq,h2,h3,lh)-Ampax0(1,hq,h2,h3,lh))/two
      endif

c--- now interfere amplitudes
      mqqb_vec0(hq,lh)=mqqb_vec0(hq,lh)+fac/xnsq*(
     &  conjg(m(1))*(xn-4._dp/xn)*ml_vec0(1)
     & +conjg(m(2))*(xn-4._dp/xn)*ml_vec0(2))

      mqqb_vect(hq,lh)=mqqb_vect(hq,lh)+fac/xnsq*(
     &  conjg(m(1))*(xn-4._dp/xn)*ml_vect(1)
     & +conjg(m(2))*(xn-4._dp/xn)*ml_vect(2))

      mqqb_ax(hq,lh)=mqqb_ax(hq,lh)+fac/xnsq*(
     &  conjg(m(1))*(
     &    (xn-2._dp/xn)*ml1_ax(1)-2._dp/xn*ml1_ax(2)+one/xn*ml2_ax(1))
     & +conjg(m(2))*(
     &    (xn-2._dp/xn)*ml1_ax(2)-2._dp/xn*ml1_ax(1)+one/xn*ml2_ax(2)))
      
      enddo
      enddo

!      write(6,*) 'hq,lh,mqqb_ax(hq,lh)',hq,lh,mqqb_ax(hq,lh)

      enddo
      enddo
      
!      pause

      return
      end
      
