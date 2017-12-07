!CW Sept 17
      subroutine qqb_zh_bb_nlopxd_ac_v(p,msqv)
      implicit none
      include 'types.f'
      
      integer:: j,k
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
      real(dp):: p(mxpart,4),
     & msqv(-nf:nf,-nf:nf),msq(-nf:nf,-nf:nf),
     & dot,virt,xl12
      real(dp) :: s56,q1423,q2413,fac,hdecay,s567
      real(dp) :: v2(2),prop

      msqv(:,:)=zip




      xl12=log(two*dot(p,1,2)/musq)
      scheme='dred'
      call spinoru(7,p,za,zb)


!------ VIRTUAL PRODUCTION 
      v2(1)=l1
      v2(2)=r1
c----- tree-level amplitude for ZH production and LO H decay 
c---calculate the 2 Z propagators
      prop=     ((s(1,2)-zmass**2)**2+(zmass*zwidth)**2)
      prop=prop*((s(3,4)-zmass**2)**2+(zmass*zwidth)**2)
      fac=xn*4._dp*(xw/(1._dp-xw))**2*gwsq**3*wmass**2/prop
      
c---  calculate the multiple of the lowest order
      virt=ason2pi*cf*(-2._dp*(epinv*epinv2-epinv*xl12+half*xl12**2)
     &                 -3._dp*(epinv-xl12)
     &                 +pisq-7._dp)


!---- H=>bbg decay
      call Hbbg_decay(p,5,6,7,hdecay)
!----- 3 body decay
      s567=s(5,6)+s(5,7)+s(6,7)
      hdecay=hdecay/((s567-hmass**2)**2+(hmass*hwidth)**2)
!-----
C     OMBINATION 
      
      q1423=aveqq*virt*hdecay*fac*s(1,4)*s(2,3)
      q2413=aveqq*virt*hdecay*fac*s(2,4)*s(1,3)
         
      do j=-nf,nf
      if (j == 0) go to 40
      k=-j
      if ((j > 0) .and. (k < 0)) then
      msqv(j,k)=
     &  +((l(j)*v2(1))**2+(r(j)*v2(2))**2)*q1423
     &  +((l(j)*v2(2))**2+(r(j)*v2(1))**2)*q2413
      elseif ((j < 0) .and. (k > 0)) then
      msqv(j,k)=
     &  +((l(k)*v2(1))**2+(r(k)*v2(2))**2)*q2413
     &  +((l(k)*v2(2))**2+(r(k)*v2(1))**2)*q1423
      endif
 40   continue
      enddo
      return
      end

