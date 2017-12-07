!---  CW Sept 17
!---  this is a piece of the NLO x NLO  qqb => Z (H=>bb+g) process
!---  It corresponds to : The LO qqb => ZH  x the above cut H=>bb+g (i.e. virtual - real)

      
      subroutine qqb_zh_bb_nlopxd_ac(p,msq)
      implicit none 
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'qcdcouple.f'
      include 'scale.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'scheme.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'zcouple.f'
      include 'hbbparams.f'
      real(dp):: p(mxpart,4),msq(-nf:nf,-nf:nf)
      real(dp):: s56,q1423,q2413,prop,fac,hdecay,s567
      real(dp):: v2(2)
      integer :: j,k
      msq(:,:)=zip


      call spinoru(7,p,za,zb)


!------ LO PRODUCTION 
      v2(1)=l1
      v2(2)=r1
c----- tree-level amplitude for ZH production and LO H decay 
c---calculate the 2 Z propagators
      prop=     ((s(1,2)-zmass**2)**2+(zmass*zwidth)**2)
      prop=prop*((s(3,4)-zmass**2)**2+(zmass*zwidth)**2)
      fac=xn*4._dp*(xw/(1._dp-xw))**2*gwsq**3*wmass**2/prop
      

!---- H=>bbg decay
      call Hbbg_decay(p,5,6,7,hdecay)
!----- 3 body decay
      s567=s(5,6)+s(5,7)+s(6,7)
      
      hdecay=hdecay/((s567-hmass**2)**2+(hmass*hwidth)**2)
      
     
      q1423=aveqq*hdecay*fac*s(1,4)*s(2,3)
      q2413=aveqq*hdecay*fac*s(2,4)*s(1,3)

!      write(6,*) q1423,hdecay,fac,s(1,4),s(2,3)
!      pause
      do j=-nf,nf
      if (j == 0) go to 40
      k=-j
      if ((j > 0) .and. (k < 0)) then
      msq(j,k)=
     &  +((l(j)*v2(1))**2+(r(j)*v2(2))**2)*q1423
     &  +((l(j)*v2(2))**2+(r(j)*v2(1))**2)*q2413
      elseif ((j < 0) .and. (k > 0)) then
      msq(j,k)=
     &  +((l(k)*v2(1))**2+(r(k)*v2(2))**2)*q2413
     &  +((l(k)*v2(2))**2+(r(k)*v2(1))**2)*q1423
      endif
 40   continue
      enddo     
      return
      end
