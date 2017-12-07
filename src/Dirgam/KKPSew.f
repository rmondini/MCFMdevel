! A collection of routines to implement the EW corrections to direct photon
! production as calculated in:
! J. Kuhn, A. Kulesza, S. Pozzorini, M. Schulze,
! hep-ph/0507178 and hep-ph/0508253

! Non-Abelian terms

! hep-ph/0508253, Eq. (18)
! Note: this function differs in the photon+jet case from the Z+jet case
!  that is given in hep-ph/0507178, Eq. (64) line one
      function g0N(shat,that,uhat,mw2)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'epinv.f'
      include 'masses.f'
      real(dp):: g0N,shat,that,uhat,mw2,
     & ashat,athat,auhat,DeltaUVbar,theta
C--- define Heaviside theta function (=1 for x>0) and (0 for x < 0)
      theta(shat)=half+half*sign(one,shat)
      
      ashat=abs(shat)
      athat=abs(that)
      auhat=abs(uhat)

c--- Eq. (14) - only include 1/epsilon since it should cancel
      DeltaUVbar=epinv
      
      g0N=two*(DeltaUVbar-log(mw2/zmass**2))
     & +log(ashat/mw2)**2-log(athat/mw2)**2-log(auhat/mw2)**2
     & +log(athat/auhat)**2-1.5_dp*(log(athat/ashat)**2+log(auhat/ashat)**2)
     & -pisq/two*theta(-shat)
     
      return
      end

      
! hep-ph/0507178, Eq. (64) line two
      function g1N(shat,that,uhat)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'epinv.f'
      include 'masses.f'
      real(dp):: g1N,shat,that,uhat,
     & ashat,athat,auhat,theta
C--- define Heaviside theta function (=1 for x>0) and (0 for x < 0)
      theta(shat)=half+half*sign(one,shat)
      
      ashat=abs(shat)
      athat=abs(that)
      auhat=abs(uhat)
      
      g1N=half*(log(auhat/ashat)**2-log(athat/ashat)**2)
     & +pisq/two*(theta(that)-theta(uhat))
     
      return
      end

      
! hep-ph/0507178, Eq. (64) line three
      function g2N(shat,that,uhat)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'epinv.f'
      include 'masses.f'
      real(dp):: g2N,shat,that,uhat,
     & ashat,athat,auhat,theta
C--- define Heaviside theta function (=1 for x>0) and (0 for x < 0)
      theta(shat)=half+half*sign(one,shat)
      
      ashat=abs(shat)
      athat=abs(that)
      auhat=abs(uhat)
      
      g2N=-two*(log(athat/ashat)**2+log(auhat/ashat)**2
     & +log(athat/ashat)+log(auhat/ashat))-two*pisq*theta(-shat)
     
      return
      end


! Abelian terms

! hep-ph/0507178, Eq. (65) line one
      function g0A(shat,that,uhat,mv2)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'epinv.f'
      include 'masses.f'
      real(dp):: g0A,shat,that,uhat,mv2,
     & ashat,athat,auhat,theta
C--- define Heaviside theta function (=1 for x>0) and (0 for x < 0)
      theta(shat)=half+half*sign(one,shat)
      
      ashat=abs(shat)
      athat=abs(that)
      auhat=abs(uhat)
      
      g0A=
     & -log(ashat/mv2)**2+three*log(ashat/mv2)
     & +1.5_dp*(log(athat/ashat)**2+log(auhat/ashat)**2
     &         +log(athat/ashat)+log(auhat/ashat))
     & +pisq/six*(two+three*theta(-shat))-five/two
     
      return
      end

      
! hep-ph/0507178, Eq. (65) line two
      function g1A(shat,that,uhat)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'epinv.f'
      include 'masses.f'
      real(dp):: g1A,g1N,shat,that,uhat,
     & ashat,athat,auhat
      
      ashat=abs(shat)
      athat=abs(that)
      auhat=abs(uhat)
      
      g1A=-g1N(shat,that,uhat)
     & +three/two*(log(auhat/ashat)-log(athat/ashat))
     
      return
      end

      
! hep-ph/0507178, Eq. (65) line three
      function g2A(shat,that,uhat)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'epinv.f'
      include 'masses.f'
      real(dp):: g2A,g2N,shat,that,uhat
            
      g2A=-g2N(shat,that,uhat)
     
      return
      end
      
