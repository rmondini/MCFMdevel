      subroutine get_qgHgaga_q_amps(p,i1,i2,i5,za,zb,amp)
      implicit none
      include 'types.f'
c--- Returns a series of arrays representing the dressed amp[itudes
c--- for the process q+g->Higgs->ga+ga+q computed in infinite top-mass effective theory
c---
c--- Only overall factor remaining is color: delta(A,B)
c---
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'ewcouple.f'
      include 'zcouple.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      include 'scale.f'
      include 'anom_higgs.f' 
      include 'qcdcouple.f'
      integer:: i1,i2,i5,h1,h2,h3,h5
      real(dp):: p(mxpart,4)
      complex(dp):: ampprod(2,2),amp(2,2,2,2,2),Hgaga(2,2),prodfac

c--- Amplitudes for production

!===== -ve quark, +ve gluon 
      ampprod(1,2)=za(i2,i1)*zb(i2,i5)**2/s(i1,i2)
!===== +ve quark, +ve gluon by exchanging 1 and 2 
      ampprod(2,2)=za(i1,i2)*zb(i1,i5)**2/s(i1,i2)
!===== -ve quark, -ve gluon 
      ampprod(1,1)=za(i1,i5)**2*zb(i2,i1)/s(i1,i2)
!===== +ve quark, -ve gluon 
      ampprod(2,1)=za(i2,i5)**2*zb(i1,i2)/s(i1,i2)

      prodfac=as/(3._dp*pi)*sqrt(gsq/vevsq)/rt2
      
      ampprod(:,:)=prodfac*ampprod(:,:)
      
c--- Amplitudes for decay
      call fill_amp_hgamgam(s(3,4),za,zb,3,4,Hgaga)

      amp(:,:,:,:,:)=czip
c--- Assemble: insert factor of (im) from propagator here
      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h5=1,2
      amp(h1,h2,h3,h3,h1)=im*ampprod(h1,h2)*Hgaga(h3,h3)
      enddo
      enddo
      enddo
      enddo

! rescale to allow for anomalous coupling
      amp=amp*chi_higgs**2

      return
      end
      
