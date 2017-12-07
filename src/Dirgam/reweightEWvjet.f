      function reweightEWvjet(p)
      implicit none
c--- Relative EW corrections to V+jet production as a function of transverse
c--- momentum of V and sqrts, according to Becher et al., 1509.01961
      include 'types.f'
      include 'mxpart.f'
      include 'energy.f'
      include 'kprocess.f'
      real(dp):: reweightEWvjet,p(mxpart,4),pt,pttwo,ptV,a0V,a1V,a2V,a3V,
     & Delta,x,y,fV,rts0,s1,s2
      
c--- routine can handle both photons and vector bosons;
c--- determine correct branch by examining kcase
      if ((kcase == kdirgam) .or. (kcase==kgamjet)) then
        ptV=pt(3,p)
c--- parameters for photon
        a0V=1.713_dp
        a1V=-21.68_dp
        a2V=12.16_dp
        a3V=-3.050_dp
        s1=-2.3355e-2_dp
        s2=1.2310e-3_dp
      elseif (kcase == kZ_1jet) then
        ptV=pttwo(3,4,p)
c--- parameters for Z-boson
        a0V=2.763_dp
        a1V=-40.76_dp
        a2V=21.95_dp
        a3V=-5.356_dp
        s1=-1.0376e-2_dp
        s2=5.4060e-4_dp
      elseif (kcase == kW_1jet) then
        ptV=pttwo(3,4,p)
        write(6,*) 'EW corrections not yet coded for W+jet'
        stop
      else
        write(6,*) 'EW corrections not available for this process'
        stop
      endif
      
c--- no application outside [100, 1500] GeV range
      if ((ptV < 100._dp) .or. (ptV > 1500._dp)) then
        Delta=0._dp
        fV=1._dp
        goto 99
      endif
      
c--- corrections at sqrts=7 TeV according to Eq. (2)
      x=ptV/1.e3_dp ! convert pt to TeV
      Delta=a0V+a1V*x+a2V*x**2+a3V*x**3
      
c--- scaling to general sqrts according to Eqs. (3), (5)
      rts0=7000._dp
      y=(sqrts-rts0)/rts0
      fV=1._dp+s1*y+s2*y**2

   99 continue

! Divide by 100 in next line to convert Delta from %
!      reweightEWvjet=1._dp+Delta/fV/100._dp ! Result including EW corrections
      reweightEWvjet=Delta/fV/100._dp ! EW corrections alone
      
      return
      end
      
