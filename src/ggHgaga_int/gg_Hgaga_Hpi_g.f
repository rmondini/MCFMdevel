      subroutine gg_Hgaga_Hpi_g(p,msq)
      implicit none
      include 'types.f'      
c--- Author: J. M. Campbell, March 2017
c--- Effect of the Higgs boson in gamma, gamma final states:
c--- includes both gg -> H -> ga ga + g signal process and its interference
c--- with gg -> ga ga + g NNLO contribution to continuum background
c---
c--- Higgs amplitude computed in effective theory
c--- Continuum amplitude computed with light quark loops only
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'noglue.f'
      include 'qlfirst.f'
      include 'epinv.f'
      include 'gamgamintflags.f'
      include 'zprods_decl.f'
      integer:: h1,h2,h3,h4,h5,j
      real(dp):: p(mxpart,4),msq(fn:nf,fn:nf),fac,facqg,hpropR,hpropI,dot,
     & msq_gg,msq_qg,msq_gq,msq_qqb!,msq_qbg,msq_gqb
      complex(dp):: AmpCont(2,2,2,2,2),AmpCont1(2,2,2,2,2),AmpCont2(2,2,2,2,2),
     & AmpHiggs(2,2,2,2,2),Amplo(2,2,2,2),Amp_compare(2,2,2,2),AmpHlo(2,2,2,2),
     & hprop,higgsprop,AmpHiggs1_exact(2,2,2,2),
     & AmpHiggs_qg(2,2,2,2,2),AmpHiggs_qbg(2,2,2,2,2),
     & AmpHiggs_gq(2,2,2,2,2),AmpHiggs_gqb(2,2,2,2,2),AmpHiggs_qqb(2,2,2,2,2),
     & AmpCont_qg(2,2,2,2,2),AmpCont_qbg(2,2,2,2,2),
     & AmpCont_gq(2,2,2,2,2),AmpCont_gqb(2,2,2,2,2),AmpCont_qqb(2,2,2,2,2)
      double precision, parameter:: statfac=0.5_dp

      if (qlfirst) then
        qlfirst=.false. 
        call qlinit
        epinv=0d0
      endif

c--- if noglue print warning message and stop
      if (noglue) then
         write(6,*) 'Please set noglue .false. in input file'
         stop
      endif
            
      msq(:,:)=zero

      call spinoru(5,p,za,zb)
      
! Higgs amplitudes
      call get_ggHgaga_g_amps(p,za,zb,AmpHiggs)

! Higgs amplitudes involving quarks
      call get_qgHgaga_q_amps(p,5,1,2,za,zb,AmpHiggs_qg)
      call get_qgHgaga_q_amps(p,5,2,1,za,zb,AmpHiggs_gq)
      call get_qgHgaga_q_amps(p,2,1,5,za,zb,AmpHiggs_qqb)
!      call get_qgHgaga_q_amps(p,1,5,2,za,zb,AmpHiggs_qbg)
!      call get_qgHgaga_q_amps(p,2,5,1,za,zb,AmpHiggs_gqb)

! Continuum amplitudes
      if (checkHsq) then
        continue
      else
        call get_gggaga_g_amps(p,za,zb,AmpCont)
        call get_qggaga_q_amps(p,5,1,2,za,zb,AmpCont_qg)
        call get_qggaga_q_amps(p,5,2,1,za,zb,AmpCont_gq)
        call get_qggaga_q_amps(p,2,1,5,za,zb,AmpCont_qqb)
!        call get_qggaga_q_amps(p,1,5,2,za,zb,AmpCont_qbg)
!        call get_qggaga_q_amps(p,2,5,1,za,zb,AmpCont_gqb)
      endif

!! This block of code checks the soft limit of Higgs and continuum amplitudes
!      if (1==2) then
!! LO amplitudes
!      call get_gggaga_amps_twoloop(p,Amplo,Amp_compare,.true.)
!!      call get_gggaga_amps(p,.true.,Amplo)
!      call get_ggHgaga_amps(p,.false.,AmpHlo,AmpHiggs1_exact)
!      call spinoru(5,p,za,zb)
!! check soft limit ...
!      do h1=1,2
!      do h2=1,2
!      do h3=1,2
!      do h4=1,2
!!      write(6,*) h1,h2,h3,h4,1,AmpCont(h1,h2,h3,h4,1)**2,AmpCont(h2,h1,h3,h4,1)**2
!!      write(6,*) h1,h2,h3,h4,2,AmpCont(h1,h2,h3,h4,2)**2,AmpCont(h2,h1,h3,h4,2)**2
!      write(6,*) h1,h2,h3,h4,1,(AmpCont(h1,h2,h3,h4,1))/
!     &           (sqrt(2d0*gsq)*Amplo(h1,h2,h3,h4)*zb(2,1)/(zb(2,5)*zb(5,1)))
!     &           ,(AmpHiggs(h1,h2,h3,h4,1))/
!     &           (sqrt(2d0*gsq)*AmpHlo(h1,h2,h3,h4)*zb(2,1)/(zb(2,5)*zb(5,1)))
!      write(6,*) h1,h2,h3,h4,2,(AmpCont(h1,h2,h3,h4,2))/
!     &           (sqrt(2d0*gsq)*Amplo(h1,h2,h3,h4)*za(2,1)/(za(2,5)*za(5,1)))
!     &           ,(AmpHiggs(h1,h2,h3,h4,2))/
!     &           (sqrt(2d0*gsq)*AmpHlo(h1,h2,h3,h4)*za(2,1)/(za(2,5)*za(5,1)))
!      enddo
!      enddo
!      enddo
!      enddo
!      pause
!      endif

! get Higgs propagator and real, imaginary parts
      hprop=higgsprop(two*dot(p,3,4))
      hpropR=real(hprop)
      hpropI=aimag(hprop)
      
      msq_gg=zero
      msq_qg=zero
      msq_gq=zero
      msq_qqb=zero
!      msq_qbg=zero
!      msq_gqb=zero
      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h4=1,2
      do h5=1,2
      
! DEBUG: continuum only
!        msq_gg=msq_gg+(
!     &          +AmpCont(h1,h2,h3,h4,h5)
!     &    *conjg(AmpCont(h1,h2,h3,h4,h5)))
!        msq_qg=msq_qg+(
!     &          +AmpCont_qg(h1,h2,h3,h4,h5)
!     &    *conjg(AmpCont_qg(h1,h2,h3,h4,h5)))
!        msq_gq=msq_gq+(
!     &          +AmpCont_gq(h1,h2,h3,h4,h5)
!     &    *conjg(AmpCont_gq(h1,h2,h3,h4,h5)))
!        msq_qbg=msq_qbg+(
!     &          +AmpCont_qbg(h1,h2,h3,h4,h5)
!     &    *conjg(AmpCont_qbg(h1,h2,h3,h4,h5)))
!        msq_gqb=msq_gqb+(
!     &          +AmpCont_gqb(h1,h2,h3,h4,h5)
!     &    *conjg(AmpCont_gqb(h1,h2,h3,h4,h5)))      
!      cycle

      if (checkHsq) then
! Accumulate |AmpHiggs|^2 as a check ...
        msq_gg=msq_gg+(
     &          +AmpHiggs(h1,h2,h3,h4,h5)*hprop
     &    *conjg(AmpHiggs(h1,h2,h3,h4,h5)*hprop))

        msq_qg=msq_qg+(
     &          +AmpHiggs_qg(h1,h2,h3,h4,h5)*hprop
     &    *conjg(AmpHiggs_qg(h1,h2,h3,h4,h5)*hprop))

        msq_gq=msq_gq+(
     &          +AmpHiggs_gq(h1,h2,h3,h4,h5)*hprop
     &    *conjg(AmpHiggs_gq(h1,h2,h3,h4,h5)*hprop))

        msq_qqb=msq_qqb+(
     &          +AmpHiggs_qqb(h1,h2,h3,h4,h5)*hprop
     &    *conjg(AmpHiggs_qqb(h1,h2,h3,h4,h5)*hprop))

!        msq_qbg=msq_qbg+(
!     &          +AmpHiggs_qbg(h1,h2,h3,h4,h5)*hprop
!     &    *conjg(AmpHiggs_qbg(h1,h2,h3,h4,h5)*hprop))

!        msq_gqb=msq_gqb+(
!     &          +AmpHiggs_gqb(h1,h2,h3,h4,h5)*hprop
!     &    *conjg(AmpHiggs_gqb(h1,h2,h3,h4,h5)*hprop))

      else
! Accumulate Higgs/continuum interference
        if (propimag) then
          msq_gg=msq_gg-two*hpropI*(
     &      +aimag(AmpHiggs(h1,h2,h3,h4,h5)*conjg(AmpCont(h1,h2,h3,h4,h5))))

          msq_qg=msq_qg-two*hpropI*(
     &      +aimag(AmpHiggs_qg(h1,h2,h3,h4,h5)*conjg(AmpCont_qg(h1,h2,h3,h4,h5))))

          msq_gq=msq_gq-two*hpropI*(
     &      +aimag(AmpHiggs_gq(h1,h2,h3,h4,h5)*conjg(AmpCont_gq(h1,h2,h3,h4,h5))))

          msq_qqb=msq_qqb-two*hpropI*(
     &      +aimag(AmpHiggs_qqb(h1,h2,h3,h4,h5)*conjg(AmpCont_qqb(h1,h2,h3,h4,h5))))

!          msq_qbg=msq_qbg-two*hpropI*(
!     &      +aimag(AmpHiggs_qbg(h1,h2,h3,h4,h5)*conjg(AmpCont_qbg(h1,h2,h3,h4,h5))))

!          msq_gqb=msq_gqb-two*hpropI*(
!     &      +aimag(AmpHiggs_gqb(h1,h2,h3,h4,h5)*conjg(AmpCont_gqb(h1,h2,h3,h4,h5))))

        endif
        if (propreal) then
          msq_gg=msq_gg+two*hpropR*(
     &      +real(AmpHiggs(h1,h2,h3,h4,h5)*conjg(AmpCont(h1,h2,h3,h4,h5))))

          msq_qg=msq_qg+two*hpropR*(
     &      +real(AmpHiggs_qg(h1,h2,h3,h4,h5)*conjg(AmpCont_qg(h1,h2,h3,h4,h5))))

          msq_gq=msq_gq+two*hpropR*(
     &      +real(AmpHiggs_gq(h1,h2,h3,h4,h5)*conjg(AmpCont_gq(h1,h2,h3,h4,h5))))

          msq_qqb=msq_qqb+two*hpropR*(
     &      +real(AmpHiggs_qqb(h1,h2,h3,h4,h5)*conjg(AmpCont_qqb(h1,h2,h3,h4,h5))))

!          msq_qbg=msq_qbg+two*hpropR*(
!     &      +real(AmpHiggs_qbg(h1,h2,h3,h4,h5)*conjg(AmpCont_qbg(h1,h2,h3,h4,h5))))

!          msq_gqb=msq_gqb+two*hpropR*(
!     &      +real(AmpHiggs_gqb(h1,h2,h3,h4,h5)*conjg(AmpCont_gqb(h1,h2,h3,h4,h5))))

        endif

      endif
     
      enddo
      enddo
      enddo
      enddo
      enddo

c--- overall color factors to be applied to amplitudes
      fac=avegg*xn*V*statfac
      facqg=aveqg*TR*V*statfac
      
      msq(0,0)=msq_gg*fac
      
      do j=1,nf
      msq(j,0)=msq_qg*facqg
      msq(0,j)=msq_gq*facqg
      msq(j,-j)=msq_qqb*facqg*aveqq/aveqg
      msq(-j,0)=msq(j,0)
      msq(0,-j)=msq(0,j)
      msq(-j,j)=msq(j,-j)
      enddo

      return
      end


