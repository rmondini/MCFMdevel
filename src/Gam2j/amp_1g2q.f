      subroutine amp_1g2q(h1234,j1,j2,j3,j4,j5,za,zb,A10,A01,B10,B01)
      implicit none
      include 'types.f'
      complex(dp):: A10,A01,B10,B01,denom,fac
c--- -i*Amplitude for the process
c---     0  -->  qb(j1) + q(j2) + qbar(j3) + q(j4) + g(j5) 
c--- helicities:  h1      h2       h3        h4       +     
c---
c--- Taken from "Multi-Photon Amplitudes for Next-to-Leading Order QCD"
c---  V. Del Duca, W. Kilgore and F. Maltoni, hep-ph/9910253
c---
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      integer:: h1234,j1,j2,j3,j4,j5
      integer, parameter:: hmpmp=1,hmppm=2,hpmmp=3,hpmpm=4  
          
      denom=cone/(za(j1,j2)*za(j3,j4))
      
      if     (h1234 == hmpmp) then
        fac=-za(j1,j3)**2
      elseif (h1234 == hmppm) then
        fac=za(j1,j4)**2
      elseif (h1234 == hpmmp) then
        fac=za(j2,j3)**2
      elseif (h1234 == hpmpm) then
        fac=-za(j2,j4)**2
      else
        write(6,*) 'Invalid h123 in amp_1g2q, h1234 =',h1234
	stop
      endif
      
      A10=fac*denom*za(j2,j3)/(za(j2,j5)*za(j5,j3))
      A01=fac*denom*za(j4,j1)/(za(j4,j5)*za(j5,j1))
      B10=fac*denom*za(j2,j1)/(za(j2,j5)*za(j5,j1))
      B01=fac*denom*za(j4,j3)/(za(j4,j5)*za(j5,j3))
      
      return
      end
      

      
      
