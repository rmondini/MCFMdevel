!---  CW Sept 17
!---  this is a piece of the NLO x NLO  qqb => Z (H=>bb) process
!---  It corresponds to : The LO qqb => ZH  x the below cut H=>bb (i.e.double virtual, integrated subtractions for CS production

      
      subroutine qqb_zh_bb_nlopxd_bc(p,msq)
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
      real(dp):: s56,q1423,q2413,prop,fac,hdecay
      real(dp):: v2(2)
      integer :: j,k
      msq(:,:)=zip


      call spinoru(6,p,za,zb)


!------ LO PRODUCTION 
      v2(1)=l1
      v2(2)=r1
c----- tree-level amplitude for ZH production and LO H decay 
c---calculate the 2 Z propagators
      prop=     ((s(1,2)-zmass**2)**2+(zmass*zwidth)**2)
      prop=prop*((s(3,4)-zmass**2)**2+(zmass*zwidth)**2)
      fac=xn*4._dp*(xw/(1._dp-xw))**2*gwsq**3*wmass**2/prop
      

!------ BELOW CUT FACTORIZATION H=>bb
      s56=s(5,6)
      call Hbb_bc_bit(s56,p(5,4),p(6,4),1,hdecay)      
      hdecay=hdecay/((s56-hmass**2)**2+(hmass*hwidth)**2)
!-----COMBINATION 
      
     
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
