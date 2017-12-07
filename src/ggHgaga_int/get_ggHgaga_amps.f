      subroutine get_ggHgaga_amps(p,doexact,Amp_heft,Amp_exact)
      implicit none
      include 'types.f'
c--- Returns a series of arrays representing the dressed amp[itudes
c--- for the process g+g->Higgs->ga+ga; there are:
c---         Amp_heft(h1,h2,h3,h4)   in the effective theory (mt->infty)
c---        Amp_exact(h1,h2,h3,h4)   including t,b,c mass effects
c---
c--- Note that Amp_exact returns zero if doexact=.false.
c---
c--- Only overall factor remaining is color: delta(A,B)
c---
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'msbarmasses.f'
      include 'ewcouple.f'
      include 'zcouple.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      include 'scale.f'
      include 'anom_higgs.f'
      include 'couple.f' 
      include 'qcdcouple.f'
      include 'gamgamintflags.f'
      include 'kpart.f'
      include 'first.f'
      logical doexact
      integer:: h1,h3
      real(dp):: p(mxpart,4),mc2,mb2,mt2,massfrun
      complex(dp)::
     & Mloop_cquark(2,2,2,2),Mloop_bquark(2,2,2,2),Mloop_tquark(2,2,2,2),
     & ggHmt(2,2),ggHmb(2,2),ggHmc(2,2),qlI3,C0mt,C0mb,C0mc,Hgaga(2,2),prodfac,
     & Amp_heft(2,2,2,2),Amp_exact(2,2,2,2),ggH_heft(2,2),ggH_exact(2,2)
      save mt2,mb2,mc2
!$omp threadprivate(mt2,mb2,mc2)

      Amp_heft(:,:,:,:)=czip
      Amp_exact(:,:,:,:)=czip
     
      call spinoru(4,p,za,zb)

      if (first) then
c--- run masses to appropriate scale (Higgs mass)
        if (kpart==klord) then
          mc2=massfrun(mc_msbar,hmass,amz,1)**2
          mb2=massfrun(mb_msbar,hmass,amz,1)**2
          mt2=massfrun(mt_msbar,hmass,amz,1)**2
        else
          mc2=massfrun(mc_msbar,hmass,amz,2)**2
          mb2=massfrun(mb_msbar,hmass,amz,2)**2
          mt2=massfrun(mt_msbar,hmass,amz,2)**2
!          write(6,*) 'mt(mh)=',sqrt(mt2)
!          write(6,*) 'mb(mh)=',sqrt(mb2)
!          write(6,*) 'mc(mh)=',sqrt(mc2)
        endif
        first=.false.
      endif

! S. Martin values for mq(mh)**2
!      mt2=168.2_dp**2
!      mb2=2.78_dp**2
!      mc2=0.72_dp**2

      ggH_heft(2,2)=s(1,2)/3._dp
      ggH_heft(1,1)=ggH_heft(2,2)*za(1,2)/zb(1,2)
      ggH_heft(2,2)=ggH_heft(2,2)*zb(1,2)/za(1,2)

      if (doexact) then

c--- Amplitudes for production 
        C0mt=qlI3(zip,zip,s(1,2),mt2,mt2,mt2,musq,0)
        C0mb=qlI3(zip,zip,s(1,2),mb2,mb2,mb2,musq,0)
        C0mc=qlI3(zip,zip,s(1,2),mc2,mc2,mc2,musq,0)
   
c------ top quark in the loop
        ggHmt(2,2)=mt2*(two-s(1,2)*C0mt*(1._dp-4._dp*mt2/s(1,2)))
        ggHmt(1,1)=ggHmt(2,2)*za(1,2)/zb(1,2)
        ggHmt(2,2)=ggHmt(2,2)*zb(1,2)/za(1,2)

c------ bottom quark in the loop
        ggHmb(2,2)=mb2*(two-s(1,2)*C0mb*(1._dp-4._dp*mb2/s(1,2)))
        ggHmb(1,1)=ggHmb(2,2)*za(1,2)/zb(1,2)
        ggHmb(2,2)=ggHmb(2,2)*zb(1,2)/za(1,2)

c------ charm quark in the loop
        ggHmc(2,2)=mc2*(two-s(1,2)*C0mc*(1._dp-4._dp*mc2/s(1,2)))
        ggHmc(1,1)=ggHmc(2,2)*za(1,2)/zb(1,2)
        ggHmc(2,2)=ggHmc(2,2)*zb(1,2)/za(1,2)
        
        ggH_exact(:,:)=ggHmt(:,:)+ggHmb(:,:)+ggHmc(:,:)
        
      else
      
        ggH_exact(:,:)=czip

      endif

! overall factor on production amplitudes: delta(A,B) remains
      prodfac=im*gw/wmass*gsq/(16d0*pisq)
      ggH_heft(:,:)=prodfac*ggH_heft(:,:)
      ggH_exact(:,:)=prodfac*ggH_exact(:,:)
      
c--- Amplitudes for decay
      call fill_amp_hgamgam(s(1,2),za,zb,3,4,Hgaga)

c--- Assemble: insert factor of (im) from propagator here
      do h1=1,2
      do h3=1,2
      Amp_heft(h1,h1,h3,h3)=im*ggH_heft(h1,h1)*Hgaga(h3,h3)
      Amp_exact(h1,h1,h3,h3)=im*ggH_exact(h1,h1)*Hgaga(h3,h3)
!      write(6,*) h1,h3,ggHmt(h1,h1),Hgaga(h3,h3)
      enddo
      enddo
      
! rescale to allow for anomalous couplings
      Amp_heft=Amp_heft*chi_higgs**2
      Amp_exact=Amp_exact*chi_higgs**2

      return
      end
      
