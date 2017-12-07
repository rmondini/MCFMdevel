      subroutine gg_Hgaga_Hpi_gvec(p,n,in,msq)
      implicit none
      include 'types.f'
c--- Author: J. M. Campbell, March 2017
c--- Effect of the Higgs boson in gamma, gamma final states:
c--- includes both gg -> H -> ga ga signal process and its interference
c--- with gg -> ga ga NNLO contribution to continuum background
c---
c--- It is based on the implementation of the continuum (squared) process
c--- in the routine gg_2gam_gvec.f
c--- For a description of the procedure implemented in this routine,
c--- see Bern, Dixon and Schmidt, hep-ph/0206194, Eqs.(9)-(11).
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'ewcharge.f'
      include 'zprods_com.f'
      include 'gamgamintflags.f'
C  in is the label of the momentum contracted with n
      integer:: j,k,in,h1,h2,h3,h4,i1,i2
      real(dp):: msq(-nf:nf,-nf:nf)
      real(dp):: n(4),p(mxpart,4),gg,hpropR,hpropI,dot
      complex(dp):: amp(2,2,2,2),ampCont1(2,2,2,2),ampCont2(2,2,2,2),phase,phasec,
     & AmpHiggs1(2,2,2,2),ggH_cquark(2,2,2,2),ggH_bquark(2,2,2,2),ggH_tquark(2,2,2,2),
     & hprop,higgsprop,amp1(2,2,2,2),amp2(2,2,2,2),csum,AmpHiggs1_exact(2,2,2,2)
      real(dp),parameter::statfac=0.5_dp

      msq(:,:)=0._dp

c--- this amplitude doesn't contain any n dependence
c--- (it has been explicitly written out), and the emitted momentum
c--- (label 5) has been passed in from dipolesub via "n"
      do j=1,4
      p(5,j)=n(j)
      enddo

      call spinoru(5,p,za,zb)

      if     (in == 1) then
        i1=1
        i2=2
      elseif (in == 2) then
        i1=2
        i2=1
      else
        write(6,*) 'Unexpected value of in in gg_Hgaga_Hpi_gvec.f: in=',in
      endif
      phase=za(i1,5)*zb(5,i2)*za(i2,i1)
     &    /(zb(i1,5)*za(5,i2)*zb(i2,i1))
      phasec=conjg(phase)
      
!      phasec=zb(i1,5)*za(5,i2)*zb(i2,i1)
!     &    /(za(i1,5)*zb(5,i2)*za(i2,i1))
!      phasec=phase

! Get Higgs amplitudes
      call get_ggHgaga_amps(p,.false.,AmpHiggs1,AmpHiggs1_exact)

! Get 1-loop continuum amplitudes
!      call gg_gaga_amps(1,2,3,4,amp)
      call get_gggaga_amps_twoloop(p,AmpCont1,AmpCont2,.true.)
      
!      do h1=1,2
!      do h2=1,2
!      do h3=1,2
!      do h4=1,2
!      write(6,*) h1,h2,h3,h4,amp(h1,h2,h3,h4),ampCont1(h1,h2,h3,h4)/amp(h1,h2,h3,h4)
!      enddo
!      enddo
!      enddo
!      enddo
!      pause
      
! get Higgs propagator and real, imaginary parts
      hprop=higgsprop(two*dot(p,3,4))
      hpropR=real(hprop)
      hpropI=aimag(hprop)
      
      amp1(:,:,:,:)=AmpHiggs1(:,:,:,:)
      
      gg=0._dp
      do h3=1,2
      do h4=1,2

      if (checkHsq) then
        amp2(:,:,h3,h4)=AmpHiggs1(:,:,h3,h4)
      else
        amp2(:,:,h3,h4)=AmpCont1(:,:,h3,h4)
      endif

      do h1=1,2
      do h2=1,2
      
!      if (in==1) then
!        if ((h1/=1).or.(h2/=1).or.(h3/=2).or.(h4/=2)) cycle
!      else
!        if ((h1/=1).or.(h2/=1).or.(h3/=2).or.(h4/=2)) cycle
!      endif

c--- first, LO term
!      gg=gg+abs(amp(h1,h2,h3,h4))**2
      csum=amp1(h1,h2,h3,h4)*conjg(amp2(h1,h2,h3,h4))
c--- now subtract correlated term
      if     (in == 1) then
        if (h1 == 2) then
!          gg=gg-real(phase*amp(h1,h2,h3,h4)*conjg(amp(3-h1,h2,h3,h4)))
          csum=csum-phase*amp1(h1,h2,h3,h4)*conjg(amp2(3-h1,h2,h3,h4))
        endif
        if (h1 == 1) then
!          gg=gg-real(phasec*amp(h1,h2,h3,h4)*conjg(amp(3-h1,h2,h3,h4)))
          csum=csum-phasec*amp1(h1,h2,h3,h4)*conjg(amp2(3-h1,h2,h3,h4))
        endif
      elseif (in == 2) then
        if (h2 == 2) then
!          gg=gg-real(phase*amp(h1,h2,h3,h4)*conjg(amp(h1,3-h2,h3,h4)))
          csum=csum-phase*amp1(h1,h2,h3,h4)*conjg(amp2(h1,3-h2,h3,h4))
        endif
        if (h2 == 1) then
!          gg=gg-real(phasec*amp(h1,h2,h3,h4)*conjg(amp(h1,3-h2,h3,h4)))
          csum=csum-phasec*amp1(h1,h2,h3,h4)*conjg(amp2(h1,3-h2,h3,h4))
        endif
      endif

      if (checkHsq) then
        gg=gg+abs(hprop)**2*real(csum)
      else
        if (propimag) then
          gg=gg-two*hpropI*aimag(csum)
        endif
        if (propreal) then
          gg=gg+two*hpropR*real(csum)
        endif
      endif
      
      enddo
      enddo
      enddo
      enddo
      
      msq(0,0)=avegg*V*statfac*gg

c--- divide by extra factor of two to compensate for normalization
c--- of "subv" in dipolesub.f
      msq(0,0)=msq(0,0)/2._dp

      return
      end

