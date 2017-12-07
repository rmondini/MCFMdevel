      function amp_1gam3g_mppppm(j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: amp_1gam3g_mppppm
c--- -i*Amplitude for the process
c---       0  -->  qb(j5) + q(j6) + g(j1) + g(j2) + g(j3) + gam(j4)
c--- helicities:     +       -      -       +       +        +
c---
c--- Taken from "Multi-Photon Amplitudes for Next-to-Leading Order QCD"
c---  V. Del Duca, W. Kilgore and F. Maltoni, hep-ph/9910253
c---
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      integer:: j1,j2,j3,j4,j5,j6

c--- Equation (2.9) with r=3, m=1 and labels permuted appropriately
c--- negative helicity gluon i=1
      amp_1gam3g_mppppm=      
     & za(j5,j1)*za(j6,j1)**3/(za(j5,j6)*za(j6,j1)*za(j1,j2)*za(j2,j3)*za(j3,j5))
     & *za(j6,j5)/(za(j6,j4)*za(j4,j5))
     
      return
      end

      
      function amp_1gam3g_pmpppm(j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: amp_1gam3g_pmpppm
c--- -i*Amplitude for the process
c---       0  -->  qb(j5) + q(j6) + g(j1) + g(j2) + g(j3) + gam(j4)
c--- helicities:     +       -      +       -       +        +
c---
c--- Taken from "Multi-Photon Amplitudes for Next-to-Leading Order QCD"
c---  V. Del Duca, W. Kilgore and F. Maltoni, hep-ph/9910253
c---
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      integer:: j1,j2,j3,j4,j5,j6

c--- Equation (2.9) with r=3, m=1 and labels permuted appropriately
c--- negative helicity gluon i=2
      amp_1gam3g_pmpppm=      
     & za(j5,j2)*za(j6,j2)**3/(za(j5,j6)*za(j6,j1)*za(j1,j2)*za(j2,j3)*za(j3,j5))
     & *za(j6,j5)/(za(j6,j4)*za(j4,j5))
     
      return
      end

      
      function amp_1gam3g_ppmppm(j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: amp_1gam3g_ppmppm
c--- -i*Amplitude for the process
c---       0  -->  qb(j5) + q(j6) + g(j1) + g(j2) + g(j3) + gam(j4)
c--- helicities:     +       -      +       +       -        +
c---
c--- Taken from "Multi-Photon Amplitudes for Next-to-Leading Order QCD"
c---  V. Del Duca, W. Kilgore and F. Maltoni, hep-ph/9910253
c---
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      integer:: j1,j2,j3,j4,j5,j6

c--- Equation (2.9) with r=3, m=1 and labels permuted appropriately
c--- negative helicity gluon i=3
      amp_1gam3g_ppmppm=      
     & za(j5,j3)*za(j6,j3)**3/(za(j5,j6)*za(j6,j1)*za(j1,j2)*za(j2,j3)*za(j3,j5))
     & *za(j6,j5)/(za(j6,j4)*za(j4,j5))
     
      return
      end

      
      function amp_1gam3g_pppmpm(j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: amp_1gam3g_pppmpm
c--- -i*Amplitude for the process
c---       0  -->  qb(j5) + q(j6) + g(j1) + g(j2) + g(j3) + gam(j4)
c--- helicities:     +       -      +       +       +        -
c---
c--- Taken from "Multi-Photon Amplitudes for Next-to-Leading Order QCD"
c---  V. Del Duca, W. Kilgore and F. Maltoni, hep-ph/9910253
c---
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      integer:: j1,j2,j3,j4,j5,j6

c--- Equation (2.9) with r=3, m=1 and labels permuted appropriately
c--- negative helicity photon i=4
      amp_1gam3g_pppmpm=      
     & za(j5,j4)*za(j6,j4)**3/(za(j5,j6)*za(j6,j1)*za(j1,j2)*za(j2,j3)*za(j3,j5))
     & *za(j6,j5)/(za(j6,j4)*za(j4,j5))
     
      return
      end

      
      function amp_1gam3g_mmpppm(j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: amp_1gam3g_mmpppm
c--- -i*Amplitude for the process
c---       0  -->  qb(j5) + q(j6) + g(j1) + g(j2) + g(j3) + gam(j4)
c--- helicities:     +       -      -       -       +        +
c---
c--- Taken from "Multi-Photon Amplitudes for Next-to-Leading Order QCD"
c---  V. Del Duca, W. Kilgore and F. Maltoni, hep-ph/9910253
c---
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      integer:: j1,j2,j3,j4,j5,j6
      real(dp):: t
      complex(dp):: zba
c--- statement function
      zba(j1,j2,j3,j4)=zb(j1,j2)*za(j2,j4)+zb(j1,j3)*za(j3,j4)
      
c--- Equation (3.1)
c--- note: I have interpreted, e.g. <5|1+2|6> as [5|1+2|6> in order
c--- to get the correct <| and |] required by the helicities
      amp_1gam3g_mmpppm=
     &+zba(j3,j1,j2,j5)*zba(j4,j3,j5,j2)*t(j1,j2,j6)
     & /(zb(j1,j2)*za(j4,j5)*zb(j1,j6)*za(j3,j5)*s(j2,j3)*s(j4,j6))
     &+zba(j3,j2,j5,j1)*zba(j4,j3,j5,j2)**2
     & /(za(j3,j5)*zb(j1,j6)*s(j2,j3)*s(j4,j6)*t(j1,j4,j6))
     &-za(j1,j2)*zb(j4,j5)*zba(j3,j1,j2,j6)*zba(j3,j1,j2,j5)
     & /(zb(j1,j2)*za(j4,j5)*s(j2,j3)*s(j4,j6)*t(j1,j2,j3))
      return
      end

      
      function amp_1gam3g_pmmppm(j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: amp_1gam3g_pmmppm
c--- -i*Amplitude for the process
c---       0  -->  qb(j5) + q(j6) + g(j1) + g(j2) + gam(j3) + gam(j4)
c--- helicities:     +       -      +       -         -         +
c---
c--- Taken from "Multi-Photon Amplitudes for Next-to-Leading Order QCD"
c---  V. Del Duca, W. Kilgore and F. Maltoni, hep-ph/9910253
c---
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      integer:: j1,j2,j3,j4,j5,j6
      real(dp):: t
      complex(dp):: zba
c--- statement function
      zba(j1,j2,j3,j4)=zb(j1,j2)*za(j2,j4)+zb(j1,j3)*za(j3,j4)
      
c--- Equation (3.2)
c--- note: I have interpreted, e.g. <4|3+6|2> as [4|3+6|2> in order
c--- to get the correct <| and |] required by the helicities
      amp_1gam3g_pmmppm=      
     &+za(j2,j3)*zba(j1,j2,j3,j6)**2
     & /(zb(j2,j3)*za(j5,j4)*za(j4,j6)*s(j1,j2)*t(j1,j2,j3))
     &+za(j2,j6)*zba(j1,j2,j3,j6)*zba(j1,j2,j6,j3)*t(j1,j4,j6)
     & /(zb(j2,j3)*za(j5,j4)*za(j4,j6)*s(j1,j2)*s(j1,j6)*s(j3,j5))
     &+za(j2,j3)*za(j2,j6)*zb(j4,j1)*zba(j5,j1,j4,j6)
     & /(zb(j2,j3)*za(j1,j2)*za(j4,j6)*s(j1,j6)*s(j3,j5))
     &+za(j2,j6)*zba(j1,j2,j6,j3)**2*zba(j4,j3,j5,j2)
     & /(za(j5,j4)*s(j1,j2)*s(j1,j6)*s(j3,j5)*t(j1,j2,j6))
     
      return
      end
    
      
      function amp_1gam3g_mpmppm(j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: amp_1gam3g_mpmppm
c--- -i*Amplitude for the process
c---       0  -->  qb(j5) + q(j6) + g(j1) + g(j2) + gam(j3) + gam(j4)
c--- helicities:     +       -      -       +         -         +
c---
c--- Taken from "Multi-Photon Amplitudes for Next-to-Leading Order QCD"
c---  V. Del Duca, W. Kilgore and F. Maltoni, hep-ph/9910253
c---
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      integer:: j1,j2,j3,j4,j5,j6
      real(dp):: t
      complex(dp):: zba
c--- statement function
      zba(j1,j2,j3,j4)=zb(j1,j2)*za(j2,j4)+zb(j1,j3)*za(j3,j4)
      
c--- Equation (3.3)
c--- note: I have interpreted, e.g. <4|2+5|3> as [4|2+5|3> in order
c--- to get the correct <| and |] required by the helicities
      amp_1gam3g_mpmppm= 
     &-za(j1,j3)**2*zba(j2,j1,j3,j6)**2
     & /(za(j5,j4)*za(j4,j6)*s(j1,j2)*s(j2,j3)*t(j1,j2,j3))
     &+zb(j2,j6)*zb(j4,j5)*za(j6,j1)*za(j5,j3)*zba(j2,j4,j5,j3)
     & /(zb(j1,j6)*za(j5,j4)*s(j1,j2)*s(j3,j5)*t(j1,j2,j6))
     &+zb(j2,j5)**2*za(j6,j1)*zba(j4,j2,j5,j3)
     & /(zb(j1,j6)*zb(j5,j3)*za(j4,j6)*s(j2,j3)*t(j1,j4,j6))
     &+za(j1,j6)*zb(j2,j5)*zba(j2,j1,j3,j5)*zba(j2,j1,j6,j3)
     & /(zb(j1,j6)*zb(j2,j3)*za(j5,j4)*za(j4,j6)*s(j1,j2)*s(j3,j5))
     &+za(j1,j3)*zba(j2,j1,j3,j6)*zba(j2,j1,j6,j3)
     & /(zb(j1,j6)*za(j5,j4)*za(j4,j6)*s(j1,j2)*s(j2,j3))
     &+za(j1,j6)*zb(j2,j5)*za(j1,j3)*zba(j4,j2,j5,j3)
     & /(zb(j1,j6)*za(j4,j6)*za(j2,j1)*s(j2,j3)*s(j3,j5))

      return
      end
      
