      subroutine gg_Hgaga_Hpi(p,msq)
      implicit none
      include 'types.f'      
c--- Author: J. M. Campbell, March 2017
c--- Effect of the Higgs boson in gamma, gamma final states:
c--- includes both gg -> H -> ga ga signal process and its interference
c--- with gg -> ga ga NNLO contribution to continuum background
c--- The effect of massive bottom and top quark loops is included
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'noglue.f'
      include 'qlfirst.f'
      include 'gamgamintflags.f'
      integer:: h1,h2,h3,h4
      real(dp):: p(mxpart,4),msq(fn:nf,fn:nf),msqgg,fac,hpropR,hpropI,dot
      complex(dp):: AmpCont(2,2,2,2),AmpCont1(2,2,2,2),AmpCont2(2,2,2,2),
     & AmpHiggs1(2,2,2,2),AmpHiggs2(2,2,2,2),Amp_compare(2,2,2,2),
     & ggH_cquark(2,2,2,2),ggH_bquark(2,2,2,2),ggH_tquark(2,2,2,2),
     & Acont,Ahiggs,hprop,higgsprop,AmpHiggs1_exact(2,2,2,2)
      logical:: includegens1and2,includebottom,includetop,masslesscont
      double precision, parameter:: statfac=0.5_dp

      if (qlfirst) then
        qlfirst=.false. 
        call qlinit
      endif

c--- if noglue print warning message and stop
      if (noglue) then
         write(6,*) 'Please set noglue .false. in input file'
         stop
      endif
            
      msq(:,:)=zero

! Higgs amplitudes
      if (useHEFT) then
        call get_ggHgaga_amps(p,.false.,AmpHiggs1,AmpHiggs1_exact)
      else
        call get_ggHgaga_amps(p,.true.,AmpHiggs1,AmpHiggs1_exact)
        AmpHiggs1(:,:,:,:)=AmpHiggs1_exact(:,:,:,:)
      endif

! Continuum amplitudes
! masslesscont=.false. : return AmpCont1 that contains mass effects
! masslesscont=.true.  : massless loops only, required for 2-loop poles
!      masslesscont=.true. ! useful for checking I1ggtogamgam cancels
!      call get_gggaga_amps(p,masslesscont,AmpCont1)
!      Amp_compare=AmpCont1

! this is a quicker routine since it contains only massless result!
      call get_gggaga_amps_twoloop(p,AmpCont1,AmpCont2,.true.)

!      do h1=1,2
!      do h2=1,2
!      do h3=1,2
!      do h4=1,2
!      if (abs(Amp_compare(h1,h2,h3,h4)) /= zip) then
!      write(6,*) h1,h2,h3,h4,AmpCont1(h1,h2,h3,h4)/Amp_compare(h1,h2,h3,h4)
!      endif
!      enddo
!      enddo
!      enddo
!      enddo
!      pause
      
! get Higgs propagator and real, imaginary parts
      hprop=higgsprop(two*dot(p,1,2))
      hpropR=real(hprop)
      hpropI=aimag(hprop)
      
      msqgg=zero
      do h1=1,2
      do h2=1,2
      do h3=1,2

!      if ((h1/=1).or.(h2/=1).or.(h3/=2)) cycle

      if (checkHsq) then
! Accumulate |AmpHiggs|^2 as a check ...
        msqgg=msqgg+(
     &          +AmpHiggs1(h1,h2,h3,h3)*hprop
     &    *conjg(AmpHiggs1(h1,h2,h3,h3)*hprop)
     &    )
     
      else     
! Accumulate Higgs/continuum interference
        if (propimag) then
          msqgg=msqgg-two*hpropI*(
     &      +aimag(AmpHiggs1(h1,h2,h3,h3)*conjg(AmpCont1(h1,h2,h3,h3))))
        endif
        if (propreal) then
          msqgg=msqgg+two*hpropR*(
     &      +real(AmpHiggs1(h1,h2,h3,h3)*conjg(AmpCont1(h1,h2,h3,h3))))
        endif
      endif
     
      enddo
      enddo
      enddo

c--- overall factor extracted (c.f. getggZZamps.f and getggHZZamps.f )
      fac=avegg*V*statfac
      
      msq(0,0)=msqgg*fac

      return
      end
