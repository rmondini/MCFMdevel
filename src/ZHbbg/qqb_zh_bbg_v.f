!---- routine for calculating

!----- p  p  = > Z H ( => b b + g) virtual
!-----this routine also calls integrated dipoles and adds them together
!---- producing the finite ME. Note that since all corrections are ff type
!---- theres no need for convolution integrals with PDFs
!---- CW March 17
      
      subroutine qqb_zh_bbg_v(p,msq)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'nf.f'
      real(dp) :: p(mxpart,4),msq(-nf:nf,-nf:nf) 
      include 'cplx.h'
      include 'masses.f'
      include 'ewcouple.f'
      include 'zcouple.f'
      include 'hdecaymode.f'
      include 'hbbparams.f'
      include 'nflav.f'
      integer:: j,k
      real(dp):: s,prop,fac,q1423,q2413,s56,v2(2)
      real(dp):: hdecay,msqhgamgam
      real(dp):: hdecay_z,hdecay_v
      real(dp):: s567,test
      
      s(j,k)=2*(p(j,4)*p(k,4)-p(j,1)*p(k,1)-p(j,2)*p(k,2)-p(j,3)*p(k,3))

      v2(1)=l1
      v2(2)=r1

      msq(:,:)=zip

c---calculate the 2 Z propagators
      prop=     ((s(1,2)-zmass**2)**2+(zmass*zwidth)**2)
      prop=prop*((s(3,4)-zmass**2)**2+(zmass*zwidth)**2)

      fac=xn*4._dp*(xw/(1._dp-xw))**2*gwsq**3*wmass**2/prop
      
!---- H=>bbg decay
      call Hbbg_decay_v(p,5,6,7,hdecay_v)
!-----integrated dipoles
      call Hbbg_decay_z(p,5,6,7,hdecay_z)
      
!----- 3 body decay
      s567=s(5,6)+s(5,7)+s(6,7)
      
      hdecay=(hdecay_v+hdecay_z)/((s567-hmass**2)**2+(hmass*hwidth)**2)

      fac=fac*hdecay


!-----LO ZH production
      q1423=aveqq*fac*s(1,4)*s(2,3)
      q2413=aveqq*fac*s(2,4)*s(1,3)

    
      
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
