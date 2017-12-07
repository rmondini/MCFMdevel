      subroutine amp_1gam1g2q_mpmppp(j1,j2,j3,j4,j5,j6,za,zb,Q12,Q34,
     & A10,A01,B10,B01)
      implicit none
      include 'types.f'
      complex(dp):: A10,A01,B10,B01,fac
c--- -i*Amplitude for the process
c---     0  -->  qb(j1) + q(j2) + qbar(j3) + q(j4) + g(j5) + gam(j6)
c--- helicities:  -       +        -         +       +        +
c---  quark charges q(j2) -> Q12,  q(j4) -> Q34
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
      real(dp):: Q12,Q34
      integer, parameter:: hmpmp=1,hmppm=2,hpmmp=3,hpmpm=4  

c--- Equation (2.15) with m=1, r=1, (k,l)=(0,1) and (1,0) (n=6) 
      fac=      
     & (Q12*za(j2,j1)/(za(j2,j6)*za(j6,j1))
     & +Q34*za(j4,j3)/(za(j4,j6)*za(j6,j3)))
      
      call amp_1g2q(hmpmp,j1,j2,j3,j4,j5,za,zb,A10,A01,B10,B01)
      A10=fac*A10
      A01=fac*A01
      B10=fac*B10
      B01=fac*B01
          
      return
      end

      
      subroutine amp_1gam1g2q_mppmpp(j1,j2,j3,j4,j5,j6,za,zb,Q12,Q34,
     & A10,A01,B10,B01)
      implicit none
      include 'types.f'
      complex(dp):: A10,A01,B10,B01,fac
c--- -i*Amplitude for the process
c---     0  -->  qb(j1) + q(j2) + qbar(j3) + q(j4) + g(j5) + gam(j6)
c--- helicities:  -       +        +         -       +        +
c---  quark charges q(j2) -> Q12,  q(j4) -> Q34
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
      real(dp):: Q12,Q34
      integer, parameter:: hmpmp=1,hmppm=2,hpmmp=3,hpmpm=4  

c--- Equation (2.15) with m=2, r=0 (n=6) 
      fac=      
     & (Q12*za(j2,j1)/(za(j2,j6)*za(j6,j1))
     & +Q34*za(j4,j3)/(za(j4,j6)*za(j6,j3)))
      
      call amp_1g2q(hmppm,j1,j2,j3,j4,j5,za,zb,A10,A01,B10,B01)
      A10=fac*A10
      A01=fac*A01
      B10=fac*B10
      B01=fac*B01
          
      return
      end

      
      subroutine amp_1gam1g2q_pmmppp(j1,j2,j3,j4,j5,j6,za,zb,Q12,Q34,
     & A10,A01,B10,B01)
      implicit none
      include 'types.f'
      complex(dp):: A10,A01,B10,B01,fac
c--- -i*Amplitude for the process
c---     0  -->  qb(j1) + q(j2) + qbar(j3) + q(j4) + g(j5) + gam(j6)
c--- helicities:  +       -        -         +       +        +
c---  quark charges q(j2) -> Q12,  q(j4) -> Q34
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
      real(dp):: Q12,Q34
      integer, parameter:: hmpmp=1,hmppm=2,hpmmp=3,hpmpm=4  

c--- Equation (2.15) with m=2, r=0 (n=6) 
      fac=      
     & (Q12*za(j2,j1)/(za(j2,j6)*za(j6,j1))
     & +Q34*za(j4,j3)/(za(j4,j6)*za(j6,j3)))
      
      call amp_1g2q(hpmmp,j1,j2,j3,j4,j5,za,zb,A10,A01,B10,B01)
      A10=fac*A10
      A01=fac*A01
      B10=fac*B10
      B01=fac*B01
          
      return
      end

      
      subroutine amp_1gam1g2q_pmpmpp(j1,j2,j3,j4,j5,j6,za,zb,Q12,Q34,
     & A10,A01,B10,B01)
      implicit none
      include 'types.f'
      complex(dp):: A10,A01,B10,B01,fac
c--- -i*Amplitude for the process
c---     0  -->  qb(j1) + q(j2) + qbar(j3) + q(j4) + g(j5) + gam(j6)
c--- helicities:  +       -        +         -       +        +
c---  quark charges q(j2) -> Q12,  q(j4) -> Q34
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
      real(dp):: Q12,Q34
      integer, parameter:: hmpmp=1,hmppm=2,hpmmp=3,hpmpm=4  

c--- Equation (2.15) with m=2, r=0 (n=6) 
      fac=      
     & (Q12*za(j2,j1)/(za(j2,j6)*za(j6,j1))
     & +Q34*za(j4,j3)/(za(j4,j6)*za(j6,j3)))
      
      call amp_1g2q(hpmpm,j1,j2,j3,j4,j5,za,zb,A10,A01,B10,B01)
      A10=fac*A10
      A01=fac*A01
      B10=fac*B10
      B01=fac*B01
          
      return
      end

      
      subroutine amp_1gam1g2q_pmpmmp(j1,j2,j3,j4,j5,j6,za,zb,Q12,Q34,
     & A10,A01,B10,B01)
      implicit none
      include 'types.f'
      complex(dp):: A10,A01,B10,B01
c--- -i*Amplitude for the process
c---     0  -->  qb(j1) + q(j2) + qbar(j3) + q(j4) + gam(j5) + gam(j6)
c--- helicities:  +       -        +         -        -         +
c---  quark charges q(j2) -> Q12,  q(j4) -> Q34
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
      real(dp):: Q12,Q34
      complex(dp):: f1_1gam1g2q,f2_1gam1g2q,g1_2gam2q,g2_2gam2q

c--- Equation (3.5)
      A01=      
     & +Q12*f1_1gam1g2q(j1,j2,j3,j4,j5,j6,za,zb)
     & +Q34*f2_1gam1g2q(j3,j4,j1,j2,j5,j6,za,zb)
      A10=      
     & +Q12*f2_1gam1g2q(j1,j2,j3,j4,j5,j6,za,zb)
     & +Q34*f1_1gam1g2q(j3,j4,j1,j2,j5,j6,za,zb)
      B10=      
     & +Q12*g1_2gam2q(j1,j2,j3,j4,j5,j6,za,zb)
     & +Q34*g2_2gam2q(j3,j4,j1,j2,j5,j6,za,zb)
      B01=      
     & +Q12*g2_2gam2q(j1,j2,j3,j4,j5,j6,za,zb)
     & +Q34*g1_2gam2q(j3,j4,j1,j2,j5,j6,za,zb)
         
      return
      end


      subroutine amp_1gam1g2q_mpmpmp(j1,j2,j3,j4,j5,j6,za,zb,Q12,Q34,
     & A10,A01,B10,B01)
      implicit none
      include 'types.f'
      complex(dp):: A10,A01,B10,B01
c--- -i*Amplitude for the process
c---     0  -->  qb(j1) + q(j2) + qbar(j3) + q(j4) + gam(j5) + gam(j6)
c--- helicities:  -       +        -         +        -         +
c---  quark charges q(j2) -> Q12,  q(j4) -> Q34
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
      real(dp):: Q12,Q34
      complex(dp):: f1_1gam1g2q,f2_1gam1g2q,g1_2gam2q,g2_2gam2q

c--- Equation (3.13)
      A01=
     & +Q12*f2_1gam1g2q(j2,j1,j4,j3,j5,j6,za,zb)
     & +Q34*f1_1gam1g2q(j4,j3,j2,j1,j5,j6,za,zb)
      A10=
     & +Q12*f1_1gam1g2q(j2,j1,j4,j3,j5,j6,za,zb)
     & +Q34*f2_1gam1g2q(j4,j3,j2,j1,j5,j6,za,zb)
      B10=
     & +Q12*g1_2gam2q(j2,j1,j4,j3,j5,j6,za,zb)
     & +Q34*g2_2gam2q(j4,j3,j2,j1,j5,j6,za,zb)
      B01=
     & +Q12*g2_2gam2q(j2,j1,j4,j3,j5,j6,za,zb)
     & +Q34*g1_2gam2q(j4,j3,j2,j1,j5,j6,za,zb)

      return
      end


      subroutine amp_1gam1g2q_pmmpmp(j1,j2,j3,j4,j5,j6,za,zb,Q12,Q34,
     & A10,A01,B10,B01)
      implicit none
      include 'types.f'
      complex(dp):: A10,A01,B10,B01
c--- -i*Amplitude for the process
c---     0  -->  qb(j1) + q(j2) + qbar(j3) + q(j4) + gam(j5) + gam(j6)
c--- helicities:  +       -        -         +        -         +
c---  quark charges q(j2) -> Q12,  q(j4) -> Q34
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
      real(dp):: Q12,Q34
      complex(dp):: f3_1gam1g2q,f4_1gam1g2q,g1_2gam2q,g2_2gam2q

c--- Equation (3.6)
      A01=      
     & +Q12*f3_1gam1g2q(j1,j2,j3,j4,j5,j6,za,zb)
     & +Q34*f3_1gam1g2q(j4,j3,j2,j1,j5,j6,za,zb)
      A10=      
     & +Q12*f4_1gam1g2q(j1,j2,j3,j4,j5,j6,za,zb)
     & +Q34*f4_1gam1g2q(j4,j3,j2,j1,j5,j6,za,zb)
      B10=      
     & +Q12*g1_2gam2q(j1,j2,j4,j3,j5,j6,za,zb)
     & -Q34*g2_2gam2q(j4,j3,j1,j2,j5,j6,za,zb)
      B01=      
     & -Q12*g2_2gam2q(j1,j2,j4,j3,j5,j6,za,zb)
     & +Q34*g1_2gam2q(j4,j3,j1,j2,j5,j6,za,zb)
         
      return
      end


      
      subroutine amp_1gam1g2q_mppmmp(j1,j2,j3,j4,j5,j6,za,zb,Q12,Q34,
     & A10,A01,B10,B01)
      implicit none
      include 'types.f'
      complex(dp):: A10,A01,B10,B01
c--- -i*Amplitude for the process
c---     0  -->  qb(j1) + q(j2) + qbar(j3) + q(j4) + gam(j5) + gam(j6)
c--- helicities:  -       +        +         -        -         +
c---  quark charges q(j2) -> Q12,  q(j4) -> Q34
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
      real(dp):: Q12,Q34
      complex(dp):: f3_1gam1g2q,f4_1gam1g2q,g1_2gam2q,g2_2gam2q

c--- Equation (3.6) + cc, parity like Eq. (3.13)
      A01=      
     & +Q12*f3_1gam1g2q(j2,j1,j4,j3,j5,j6,za,zb)
     & +Q34*f3_1gam1g2q(j3,j4,j1,j2,j5,j6,za,zb)
      A10=      
     & +Q12*f4_1gam1g2q(j2,j1,j4,j3,j5,j6,za,zb)
     & +Q34*f4_1gam1g2q(j3,j4,j1,j2,j5,j6,za,zb)
      B10=      
     & +Q12*g1_2gam2q(j2,j1,j3,j4,j5,j6,za,zb)
     & -Q34*g2_2gam2q(j3,j4,j2,j1,j5,j6,za,zb)
      B01=      
     & -Q12*g2_2gam2q(j2,j1,j3,j4,j5,j6,za,zb)
     & +Q34*g1_2gam2q(j3,j4,j2,j1,j5,j6,za,zb)
         
      return
      end


      
      function f1_1gam1g2q(j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: f1_1gam1g2q
c--- -i*f1(j1,j2,j3,j4,j5,j6): auxiliary function for the process
c---     0  -->  qb(j1) + q(j2) + qbar(j3) + q(j4) + g(j5) + gam(j6)
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

c--- Equation (3.7)
      f1_1gam1g2q=
     &+zba(j3,j1,j6,j2)*zba(j3,j1,j6,j4)
     & /(za(j1,j6)*za(j6,j2)*zb(j3,j5)*zb(j5,j2)*s(j3,j4))
     &+za(j5,j2)*zb(j3,j1)*zba(j6,j2,j5,j4)
     & /(zb(j5,j2)*za(j6,j2)*s(j3,j4)*t(j1,j3,j4))
     &+za(j4,j5)*zba(j3,j1,j6,j2)**2
     & /(za(j1,j6)*za(j6,j2)*zb(j3,j5)*s(j3,j4)*t(j1,j2,j6))
     
      return
      end
            
            
      function f2_1gam1g2q(j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: f2_1gam1g2q
c--- -i*f1(j1,j2,j3,j4,j5,j6): auxiliary function for the process
c---     0  -->  qb(j1) + q(j2) + qbar(j3) + q(j4) + g(j5) + gam(j6)
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

c--- Equation (3.8)
      f2_1gam1g2q=
     &+zba(j1,j4,j5,j2)*zba(j3,j1,j6,j2)
     & /(za(j1,j6)*za(j6,j2)*zb(j1,j5)*zb(j5,j4)*s(j3,j4))
     &+zb(j1,j6)*za(j2,j4)*zba(j3,j1,j6,j5)
     & /(za(j1,j6)*zb(j1,j5)*s(j3,j4)*t(j1,j5,j6))
     &+za(j5,j3)*zba(j3,j1,j6,j2)**2
     & /(za(j1,j6)*za(j6,j2)*zb(j5,j4)*s(j3,j4)*t(j1,j2,j6))
     
      return
      end

            
      function f3_1gam1g2q(j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: f3_1gam1g2q
c--- -i*f1(j1,j2,j3,j4,j5,j6): auxiliary function for the process
c---     0  -->  qb(j1) + q(j2) + qbar(j3) + q(j4) + g(j5) + gam(j6)
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

c--- Equation (3.11)
      f3_1gam1g2q=
     &-t(j1,j4,j6)*zba(j4,j1,j6,j2)
     & /(za(j1,j6)*za(j6,j2)*zb(j3,j5)*zb(j5,j2)*s(j3,j4))
     &+za(j5,j4)*zba(j4,j1,j6,j2)**2
     & /(za(j1,j6)*za(j6,j2)*zb(j3,j5)*s(j3,j4)*t(j1,j2,j6))
     &+za(j5,j2)*zb(j4,j1)*zba(j6,j2,j5,j3)
     & /(zb(j5,j2)*za(j6,j2)*s(j3,j4)*t(j1,j3,j4))
     
      return
      end
            

      function f4_1gam1g2q(j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: f4_1gam1g2q
c--- -i*f1(j1,j2,j3,j4,j5,j6): auxiliary function for the process
c---     0  -->  qb(j1) + q(j2) + qbar(j3) + q(j4) + g(j5) + gam(j6)
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

c--- Equation (3.12)
      f4_1gam1g2q=
     &+zb(j4,j1)*za(j2,j3)*zba(j4,j1,j6,j2)
     & /(za(j1,j6)*za(j6,j2)*zb(j1,j5)*zb(j5,j4)*s(j3,j4))
     &+zb(j1,j6)*za(j2,j3)*zba(j4,j1,j6,j5)
     & /(za(j1,j6)*zb(j1,j5)*s(j3,j4)*t(j1,j5,j6))
     &+za(j3,j5)*zba(j4,j1,j6,j2)**2
     & /(za(j1,j6)*za(j6,j2)*zb(j5,j4)*s(j3,j4)*t(j1,j2,j6))
     
      return
      end
            
            
            
            
