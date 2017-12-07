      subroutine qqb_zh_bbg_gs(p,msq)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'nf.f'
      include 'ptilde.f'

      real(dp) :: p(mxpart,4),msq(1:maxd,-nf:nf,-nf:nf) 
      include 'cplx.h'
      include 'masses.f'
      include 'ewcouple.f'
      include 'zcouple.f'
      include 'hdecaymode.f'
      include 'hbbparams.f'
      integer:: j,k,nd
      real(dp):: s,prop,fac,q1423,q2413,s56,v2(2)
      real(dp):: hdecay(maxd),msqhgamgam
      real(dp):: s5678,test

      
      s(j,k)=2*(p(j,4)*p(k,4)-p(j,1)*p(k,1)-p(j,2)*p(k,2)-p(j,3)*p(k,3))

      v2(1)=l1
      v2(2)=r1

!      ndmax=6
      msq(:,:,:)=zip

c---calculate the 2 Z propagators
      prop=     ((s(1,2)-zmass**2)**2+(zmass*zwidth)**2)
      prop=prop*((s(3,4)-zmass**2)**2+(zmass*zwidth)**2)

      fac=xn*4._dp*(xw/(1._dp-xw))**2*gwsq**3*wmass**2/prop
      
!---- H=>bbg decay subtractions
      call Hbbg_decay_gs(p,5,6,7,8,hdecay)
!---- NOTE that ALL DIPOLES SHOULD CONSERVER P5679 = ptilde 567
!----- SO can cheat and use the real-PS here to calcualte Higgs prop
!-----4 body decay
      
!      write(6,*) hdecay
!      pause
      s5678=s(5,6)+s(5,7)+s(5,8)+s(6,7)+s(6,8)+s(7,8)

!      write(6,*) test,hdecay,test/hdecay
      
      hdecay(:)=hdecay(:)/((s5678-hmass**2)**2+(hmass*hwidth)**2)

!      fac=fac*hdecay


!-----LO ZH production
      q1423=aveqq*fac*s(1,4)*s(2,3)
      q2413=aveqq*fac*s(2,4)*s(1,3)

    
      do nd=1,ndmax
         do j=-nf,nf
            if (j == 0) go to 40
            k=-j
            if ((j > 0) .and. (k < 0)) then
               msq(nd,j,k)=
     &              +((l(j)*v2(1))**2+(r(j)*v2(2))**2)*q1423*hdecay(nd)
     &              +((l(j)*v2(2))**2+(r(j)*v2(1))**2)*q2413*hdecay(nd)
            elseif ((j < 0) .and. (k > 0)) then
               msq(nd,j,k)=
     &              +((l(k)*v2(1))**2+(r(k)*v2(2))**2)*q2413*hdecay(nd)
     &              +((l(k)*v2(2))**2+(r(k)*v2(1))**2)*q1423*hdecay(nd)
            endif
 40         continue
         enddo
      enddo
      
      return
      end 
