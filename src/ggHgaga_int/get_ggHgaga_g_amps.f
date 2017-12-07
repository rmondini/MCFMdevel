      subroutine get_ggHgaga_g_amps(p,za,zb,amp)
      implicit none
      include 'types.f'
c--- Returns a series of arrays representing the dressed amp[itudes
c--- for the process g+g->Higgs->ga+ga+g computed in infinite top-mass effective theory
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
      integer:: h1,h2,h3,h5
      real(dp):: p(mxpart,4),mc2,mb2,mt2
      complex(dp):: ampprod(2,2,2),amp(2,2,2,2,2),
     & Mloop_cquark(2,2,2,2),Mloop_bquark(2,2,2,2),Mloop_tquark(2,2,2,2),
     & ggHmt(2,2),ggHmb(2,2),ggHmc(2,2),qlI3,C0mt,C0mb,C0mc,prop34,
     & Hgaga(2,2),sinthw,higgsprop,prodfac

c--- propagator factors
      prop34=higgsprop(s(3,4))

c--- Amplitudes for production
      ampprod(1,1,1)=im*(za(3,4)*zb(4,3))**2/(zb(1,2)*zb(2,5)*zb(5,1)) 
      ampprod(1,1,2)=im*za(1,2)**3/(za(2,5)*za(5,1))
      ampprod(1,2,1)=im*za(5,1)**3/(za(1,2)*za(2,5))
      ampprod(2,1,1)=im*za(2,5)**3/(za(5,1)*za(1,2))

      ampprod(2,2,2)=im*(zb(3,4)*za(4,3))**2/(za(1,2)*za(2,5)*za(5,1)) 
      ampprod(2,2,1)=im*zb(1,2)**3/(zb(2,5)*zb(5,1))
      ampprod(2,1,2)=im*zb(5,1)**3/(zb(1,2)*zb(2,5))
      ampprod(1,2,2)=im*zb(2,5)**3/(zb(5,1)*zb(1,2))

      prodfac=as/(3._dp*pi)*sqrt(gsq/vevsq)/rt2
      
      ampprod(:,:,:)=prodfac*ampprod(:,:,:)
      
c--- Amplitudes for decay
      call fill_amp_hgamgam(s(3,4),za,zb,3,4,Hgaga)

      amp(:,:,:,:,:)=czip
c--- Assemble: insert factor of (im) from propagator here
      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h5=1,2
      amp(h1,h2,h3,h3,h5)=im*ampprod(h1,h2,h5)*Hgaga(h3,h3)!*prop34
      enddo
      enddo
      enddo
      enddo

! rescale to allow for anomalous coupling
      amp=amp*chi_higgs**2

      return
      end
      
