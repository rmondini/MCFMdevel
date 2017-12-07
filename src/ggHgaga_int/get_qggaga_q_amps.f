      subroutine get_qggaga_q_amps(p,i1,i2,i5,za,zb,Amp)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'nf.f'
      include 'ewcharge.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'zprods_decl.f'
      real(dp):: p(mxpart,4)
      integer i1,i2,i5,h1,h2,h3,h4
      real(dp):: Qsumsq,fac
      complex(dp):: Amp(2,2,2,2,2),permamp(2,2,2,2)
     
      call amp_virt_nf_gmgmjt(i1,i2,3,4,i5,za,zb,permamp)

! overall factor to be applied to amplitude
      Qsumsq=Q(1)**2+Q(2)**2+Q(3)**2+Q(4)**2+Q(5)**2
      fac=Qsumsq*ason2pi*esq*sqrt(gsq)*rt2

! Need to flip sign of amplitudes with negative helicity photons.
! This is due to an extra sign when using za<->zb charge
! conjugation since these amplitudes contain a fermion line
! (I believe this has no effect elsewhere, e.g. gg2gam, since only
!  interference is sensitive to these signs)
      permamp(:,1,1,:)=-permamp(:,1,1,:)
   
c--- Need to shuffle helicities appropriately:
c--- routine fills permamp with (hq, hgam, hgam, hg)
c--- want to returm amp with (hq, hg, hgam, hgam, hq) 
      Amp(:,:,:,:,:)=czip
      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h4=1,2
      Amp(h1,h2,h3,h4,h1)=fac*permamp(h1,h3,h4,h2)
      enddo
      enddo
      enddo
      enddo
      
      return
      end
