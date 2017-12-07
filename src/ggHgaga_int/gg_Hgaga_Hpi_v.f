      subroutine gg_Hgaga_Hpi_v(p,msq)
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
      include 'scheme.f'
      include 'gamgamintflags.f'
      integer:: h1,h2,h3,h4
      real(dp):: p(mxpart,4),msq(fn:nf,fn:nf),msqgg,fac,hpropR,hpropI,dot
      complex(dp):: AmpCont(2,2,2,2),AmpCont1(2,2,2,2),AmpCont2(2,2,2,2),
     & AmpHiggs1(2,2,2,2),AmpHiggs2(2,2,2,2),
     & ggH_cquark(2,2,2,2),ggH_bquark(2,2,2,2),ggH_tquark(2,2,2,2),
     & Acont,Ahiggs,hprop,higgsprop,AmpHiggs1_exact(2,2,2,2)
      logical:: includegens1and2,includebottom,includetop,masslesscont
      double precision, parameter:: statfac=0.5_dp

      scheme='dred'

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
      call get_ggHgaga_amps(p,.true.,AmpHiggs1,AmpHiggs1_exact)
      call get_ggHgaga_amps_twoloop(AmpHiggs1,AmpHiggs2)

! Continuum amplitudes for massless quarks
      call get_gggaga_amps_twoloop(p,AmpCont1,AmpCont2,.false.)
      
! get Higgs propagator and real, imaginary parts
      hprop=higgsprop(two*dot(p,1,2))
      hpropR=real(hprop)
      hpropI=aimag(hprop)
      
      msqgg=zero
      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h4=1,2

      if (checkHsq) then
! Accumulate |AmpHiggs|^2 as a check ...
        msqgg=msqgg+two*real(
     &          +AmpHiggs1(h1,h2,h3,h4)*hprop
     &    *conjg(AmpHiggs2(h1,h2,h3,h4)*hprop)
     &    )

! Replace LO (1-loop) HEFT contribution with exact treatment of masses 
        msqgg=msqgg+abs(hprop)**2*(
     &     AmpHiggs1_exact(h1,h2,h3,h4)*conjg(AmpHiggs1_exact(h1,h2,h3,h4))
     &    -AmpHiggs1(h1,h2,h3,h4)*conjg(AmpHiggs1(h1,h2,h3,h4)))

      else

        if (propimag) then
          msqgg=msqgg-two*hpropI*(
! Two-loop continuum
     &      +aimag(AmpHiggs1(h1,h2,h3,h4)*conjg(AmpCont2(h1,h2,h3,h4)))
! Two-loop Higgs
     &      +aimag(AmpHiggs2(h1,h2,h3,h4)*conjg(AmpCont1(h1,h2,h3,h4))))
        endif
        if (propreal) then
          msqgg=msqgg+two*hpropR*(
! Two-loop continuum
     &      +real(AmpHiggs1(h1,h2,h3,h4)*conjg(AmpCont2(h1,h2,h3,h4)))
! Two-loop Higgs
     &      +real(AmpHiggs2(h1,h2,h3,h4)*conjg(AmpCont1(h1,h2,h3,h4))))
        endif

! Replace LO (1-loop) HEFT contribution with exact treatment of masses 
        if (propimag) then
          msqgg=msqgg-two*hpropI*(
     &      +aimag((AmpHiggs1_exact(h1,h2,h3,h4)-AmpHiggs1(h1,h2,h3,h4))
     &       *conjg(AmpCont1(h1,h2,h3,h3))))
        endif
        if (propreal) then
          msqgg=msqgg+two*hpropR*(
     &      +real((AmpHiggs1_exact(h1,h2,h3,h4)-AmpHiggs1(h1,h2,h3,h4))
     &       *conjg(AmpCont1(h1,h2,h3,h3))))
        endif

      endif

      enddo
      enddo
      enddo
      enddo

c--- overall factor extracted (c.f. getggZZamps.f and getggHZZamps.f )
      fac=avegg*V*statfac
      
      msq(0,0)=msqgg*fac

      return
      end
