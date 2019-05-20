      subroutine lumxmsq_zh_hbb(p,xx,QB,order,xmsq)
      implicit none
      include 'types.f'
c---- Matrix element for H=> bb decay NNLO production
c---- takes in 
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'scale.f'
      include 'facscale.f'
      include 'zcouple.f'
      include 'ewcharge.f'
      include 'scet_const.f'
      include 'hbbparams.f'
      integer:: j,k,ih1,ih2,m,n,order
      real(dp),intent(out):: xmsq
      real(dp):: p(mxpart,4),fac,
     & xx(2),soft1(-1:1),soft2(-1:3),hard(2),
     & beama0(-5:5),beamb0(-5:5),
     & beama1(-5:5,-1:1),beamb1(-5:5,-1:1),
     & beama2(-5:5,-1:3),beamb2(-5:5,-1:3),
     & z1,z2,QB(2),lum0,lum1(-1:1),lum2(-1:3),bit,lobit,
     &     msq(-nf:nf,-nf:nf),assemble_dec
      real(dp) :: L12,s56,q2413,q1423,hdecay,v2(2)
      real(dp) :: jeta1(-1:1),jeta2(-1:3),jetb1(-1:1),jetb2(-1:3)
!      complex(dp):: prop,qqb,qbq
      real(dp) :: prop
      common/density/ih1,ih2
      include 'cplx.h'
      logical incpowcorr
!      common/incpowrcorr/incpowcorr
      real(dp) :: corrbit,msqpow(-5:5,-5:5)
      logical nfonly
      common/nfonly/nfonly
!$omp threadprivate(/nfonly/)

      incpowcorr=.false.
      nfonly=.false.
      
      call spinoru(6,p,za,zb)
      
      v2(1)=l1
      v2(2)=r1
c----- tree-level amplitude for ZH production and LO H decay 
c---calculate the 2 Z propagators
      prop=     ((s(1,2)-zmass**2)**2+(zmass*zwidth)**2)
      prop=prop*((s(3,4)-zmass**2)**2+(zmass*zwidth)**2)
      fac=xn*4._dp*(xw/(1._dp-xw))**2*gwsq**3*wmass**2/prop
      s56=s(5,6)
!=====LO H=>bb decay
      hdecay=xn*gwsq*mb_eff**2/(4._dp*wmass**2)*2._dp*s56
      hdecay=hdecay/((s56-hmass**2)**2+(hmass*hwidth)**2)

!---- put it all together 
      fac=fac*hdecay

      q1423=aveqq*fac*s(1,4)*s(2,3)
      q2413=aveqq*fac*s(2,4)*s(1,3)
     
      
c---------------------------------
      write(*,*) 'p(1,4)',p(1,4)
      write(*,*) 'p(1,1)',p(1,1)
      write(*,*) 'p(1,2)',p(1,2)
      write(*,*) 'p(1,3)',p(1,3)
c---------------------------------


!---- constuct log for soft function 
!      L12=log(s(5,6)/p(5,4)/p(6,4)/4._dp)

      L12=zip
      
!      call soft_qqdec(order,soft1,soft2,L12)
      call softqqbis(order,soft1,soft2)

      
      
!---- hard function 
      call hard_hbb(s(5,6),musq,hard)

c--- two jets functions for b-bar 3
!       call jetq(order,two*p(5,4),jeta1,jeta2)
!       call jetq(order,two*p(6,4),jetb1,jetb2)

      call jetq(order,sqrt(s56),jeta1,jeta2)
      call jetq(order,sqrt(s56),jetb1,jetb2)

      
      if (order >= 0) then
      call fdist_select(ih1,1,xx(1),facscale,beama0)
      call fdist_select(ih2,2,xx(2),facscale,beamb0)
      endif

!      write(6,*) jeta1,jetb1
!      write(6,*) soft1
!      write(6,*) hard
!      pause
      xmsq=zip
      do j=-nf,nf
      k=-j
      if (j == 0) cycle
      bit=assemble_dec(order,
     & beama0(j),beamb0(k),jeta1,jetb1,
     & jeta2,jetb2,soft1,soft2,hard)

      
      if ((j > 0) .and. (k < 0)) then
      lobit=
     &  +((l(j)*v2(1))**2+(r(j)*v2(2))**2)*q1423
     &        +((l(j)*v2(2))**2+(r(j)*v2(1))**2)*q2413
      bit=bit*lobit
      elseif ((j < 0) .and. (k > 0)) then
      lobit=
     &  +((l(k)*v2(1))**2+(r(k)*v2(2))**2)*q2413
     &        +((l(k)*v2(2))**2+(r(k)*v2(1))**2)*q1423
      bit=bit*lobit
      endif
      
      xmsq=xmsq+bit

c---------------------------------
      write(*,*) 'xmsq=',xmsq
      stop
c---------------------------------

      enddo

      return
      end
