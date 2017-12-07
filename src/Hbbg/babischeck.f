!---  collection of results from Babis et al for H=>bb @ NNLO used for checking
!---  MCFM amplitude implementaition
  
      
      function babisHbbg_msq(i1,i2,i3)
      implicit none
      include 'types.f' 
      real(dp):: babisHbbg_msq
      integer i1,i2,i3
      include 'constants.f'
      include 'qcdcouple.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'ewcouple.f'
      include 'hbbparams.f'
   
      real(dp) :: t,mhsq,babis_NLOfac      
      t(i1,i2,i3)=s(i1,i2)+s(i2,i3)+s(i1,i3)

      babis_NLOfac=mb_eff**2/vevsq*ason2pi*2._dp
      mhsq=s(i1,i2)+s(i1,i3)+s(i2,i3)
      
      babisHbbg_msq=  (16*Cf*mhsq*Pi**2*xn*(one+ s(i1,i2)**2/mhsq**2)*
     -    (one/s(i1,i3) + one/s(i2,i3)))/(one - s(i1,i2)/mhsq)

      babisHbbg_msq=babisHbbg_msq*babis_NLOfac
      return
      end

      function babisHbbgg_msq(i1,i2,i3,i4)
      implicit none
      include 'types.f' 
      real(dp):: babisHbbgg_msq
      integer i1,i2,i3,i4
      include 'constants.f'
      include 'mxpart.f'
      include 'qcdcouple.f'
      include 'sprods_com.f'
      include 'ewcouple.f' 
      include 'hbbparams.f'
      real(dp) :: t,mhsq,babis_NNLOfac
      real(dp) :: babisA0,babisB0,babisC0
      
      t(i1,i2,i3)=s(i1,i2)+s(i2,i3)+s(i1,i3)
      babis_NNLOfac=mb_eff**2/vevsq*ason2pi**2*4._dp
      
      mhsq=s(i1,i2)+s(i1,i3)+s(i1,i4)+s(i2,i3)+s(i2,i4)+s(i3,i4)

      babisHbbgg_msq=16._dp*pisq**2*(
     &     cf*babisA0(i1,i2,i3,i4)+ xn*cf**2*babisB0(i1,i2,i3,i4)
     &     +cf*(one+2*xn*cf)*babisC0(i1,i2,i3,i4))

      babisHbbgg_msq= babisHbbgg_msq*babis_NNLOfac
      return
      end
      
      function babisHbbQQ_msq(i1,i2,i3,i4)
      implicit none
      include 'types.f' 
      real(dp):: babisHbbQQ_msq
      integer i1,i2,i3,i4
      include 'constants.f'
      include 'mxpart.f'
      include 'qcdcouple.f'
      include 'sprods_com.f'
      include 'ewcouple.f' 
      include 'hbbparams.f'
      real(dp) :: t,mhsq,babis_NNLOfac
      
      t(i1,i2,i3)=s(i1,i2)+s(i2,i3)+s(i1,i3)
      babis_NNLOfac=mb_eff**2/vevsq*ason2pi**2*4._dp
      
      mhsq=s(i1,i2)+s(i1,i3)+s(i1,i4)+s(i2,i3)+s(i2,i4)+s(i3,i4)

      babisHbbQQ_msq=  64*cf*pisq**2*xn*(2/s(i3,i4) - 1/t(i1,i3,i4) + 
     -    (-2*(mhsq + s(i1,i4)) + s(i2,i3) - s(i2,i4))/
     -     (s(i3,i4)*t(i1,i3,i4)) - 
     -    mhsq*(t(i1,i3,i4)**(-2) + t(i2,i3,i4)**(-2)) - 
     -    (2*mhsq*(s(i1,i4)/t(i1,i3,i4)**2 + 
     -         s(i2,i4)/t(i2,i3,i4)**2))/s(i3,i4) - 
     -    (2*mhsq*(s(i1,i4)/t(i1,i3,i4) - s(i2,i4)/t(i2,i3,i4))**
     -        2)/s(i3,i4)**2 - 1/t(i2,i3,i4) + 
     -    (s(i1,i3) - s(i1,i4) - 2*(mhsq + s(i2,i4)))/
     -     (s(i3,i4)*t(i2,i3,i4)) + 
     -    (2*(mhsq**2 + mhsq*(s(i1,i4) + s(i2,i4)) + 
     -         (s(i1,i4) + s(i2,i4))**2))/
     -     (s(i3,i4)*t(i1,i3,i4)*t(i2,i3,i4)) + 
     -    (2*(mhsq + 2*(s(i1,i4) + s(i2,i4)) + s(i3,i4)))/
     -     (t(i1,i3,i4)*t(i2,i3,i4)))

      babisHbbQQ_msq= babisHbbQQ_msq*babis_NNLOfac
      return
      end


      function babisC0(i1,i2,i3,i4)
      implicit none
      include 'types.f' 
      real(dp):: babisC0
      integer i1,i2,i3,i4
      include 'constants.f'
      include 'mxpart.f'
      include 'qcdcouple.f'
      include 'sprods_com.f'
      include 'ewcouple.f' 
      include 'hbbparams.f'
      real(dp) :: t,mhsq,babis_NNLOfac
      
      t(i1,i2,i3)=s(i1,i2)+s(i2,i3)+s(i1,i3)
      mhsq=s(i1,i2)+s(i1,i3)+s(i1,i4)+s(i2,i3)+s(i2,i4)+s(i3,i4)

      babisC0=(-4*(2*mhsq - 2*s(i1,i3) - s(i1,i4) - 2*s(i2,i4)))/
     -   (s(i2,i3)*s(i3,i4)) - 
     -  (4*(4*mhsq - 4*s(i1,i3) - 2*s(i2,i3) - s(i2,i4)))/
     -   (s(i1,i4)*s(i3,i4)) - 
     -  (4*(2*mhsq - 2*s(i1,i4) - 2*s(i2,i3) - s(i2,i4)))/
     -   (s(i1,i3)*s(i3,i4)) - 
     -  (4*(4*mhsq - 2*s(i1,i3) - s(i1,i4) - 4*s(i2,i3)))/
     -   (s(i2,i4)*s(i3,i4)) - 
     -  (4*(-2*mhsq**2 + 2*mhsq*s(i1,i4) - s(i1,i4)**2 + 
     -       2*mhsq*s(i2,i3) - 2*s(i1,i4)*s(i2,i3) - s(i2,i3)**2))
     -    /(s(i1,i3)*s(i2,i4)*s(i3,i4)) - 
     -  (4*(-2*mhsq**2 + 2*mhsq*s(i1,i3) - s(i1,i3)**2 + 
     -       2*mhsq*s(i2,i4) - 2*s(i1,i3)*s(i2,i4) - s(i2,i4)**2))
     -    /(s(i1,i4)*s(i2,i3)*s(i3,i4)) + 
     -  (16*mhsq*s(i1,i4)**2)/(s(i3,i4)**2*t(i1,i3,i4)**2) + 
     -  (16*mhsq*s(i1,i4))/(s(i3,i4)*t(i1,i3,i4)**2) - 
     -  (8*(-4*mhsq - 2*s(i1,i4) + s(i2,i3) - s(i2,i4)))/
     -   (s(i3,i4)*t(i1,i3,i4)) - 
     -  (4*(-2*mhsq**2 - 2*mhsq*s(i1,i4) - s(i1,i4)**2 + 
     -       2*mhsq*s(i2,i3) + 2*s(i1,i4)*s(i2,i3) - s(i2,i3)**2))
     -    /(s(i2,i4)*s(i3,i4)*t(i1,i3,i4)) - 
     -  (4*(-2*mhsq**2 + 2*mhsq*s(i1,i4) - s(i1,i4)**2 + 
     -       2*mhsq*s(i2,i4) - 2*s(i1,i4)*s(i2,i4) - s(i2,i4)**2))
     -    /(s(i2,i3)*s(i3,i4)*t(i1,i3,i4)) + 
     -  (16*mhsq*s(i2,i4)**2)/(s(i3,i4)**2*t(i2,i3,i4)**2) + 
     -  (16*mhsq*s(i2,i4))/(s(i3,i4)*t(i2,i3,i4)**2) + 
     -  (8*(4*mhsq - s(i1,i3) + s(i1,i4) + 2*s(i2,i4)))/
     -   (s(i3,i4)*t(i2,i3,i4)) - 
     -  (4*(-2*mhsq**2 + 2*mhsq*s(i1,i3) - s(i1,i3)**2 - 
     -       2*mhsq*s(i2,i4) + 2*s(i1,i3)*s(i2,i4) - s(i2,i4)**2))
     -    /(s(i1,i4)*s(i3,i4)*t(i2,i3,i4)) - 
     -  (4*(-2*mhsq**2 + 2*mhsq*s(i1,i4) - s(i1,i4)**2 + 
     -       2*mhsq*s(i2,i4) - 2*s(i1,i4)*s(i2,i4) - s(i2,i4)**2))
     -    /(s(i1,i3)*s(i3,i4)*t(i2,i3,i4)) - 
     -  (32*mhsq*s(i1,i4)*s(i2,i4))/
     -   (s(i3,i4)**2*t(i1,i3,i4)*t(i2,i3,i4)) + 
     -  (16*(-4*mhsq**2 - mhsq*s(i1,i4) - s(i1,i4)**2 - 
     -       mhsq*s(i2,i4) - 2*s(i1,i4)*s(i2,i4) - s(i2,i4)**2))/
     -     (s(i3,i4)*t(i1,i3,i4)*t(i2,i3,i4))

      return
      end
      
      function babisB0(i1,i2,i3,i4)
      implicit none
      include 'types.f' 
      real(dp):: babisB0
      integer i1,i2,i3,i4
      include 'constants.f'
      include 'mxpart.f'
      include 'qcdcouple.f'
      include 'sprods_com.f'
      include 'ewcouple.f' 
      include 'hbbparams.f'
      real(dp) :: t,mhsq,babis_NNLOfac
      
      t(i1,i2,i3)=s(i1,i2)+s(i2,i3)+s(i1,i3)
      mhsq=s(i1,i2)+s(i1,i3)+s(i1,i4)+s(i2,i3)+s(i2,i4)+s(i3,i4)

      babisB0=     16/s(i1,i3) + 8/s(i1,i4) + 16/s(i2,i3) + 8/s(i2,i4) + 
     -  (8*(-mhsq + 2*s(i1,i4) + 2*s(i2,i3) + s(i3,i4)))/
     -   (s(i1,i3)*s(i2,i4)) + 
     -  (8*(-mhsq + 2*s(i1,i3) + 2*s(i2,i4) + s(i3,i4)))/
     -   (s(i1,i4)*s(i2,i3)) + (24*mhsq)/t(i1,i3,i4)**2 - 
     -  (8*mhsq*s(i1,i4))/(s(i1,i3)*t(i1,i3,i4)**2) + 
     -  (8*mhsq*s(i3,i4))/(s(i1,i4)*t(i1,i3,i4)**2) + 
     -  40/t(i1,i3,i4) + (8*mhsq)/(s(i1,i3)*t(i1,i3,i4)) - 
     -  (8*mhsq*(2*mhsq + s(i1,i4) - s(i2,i3)))/
     -   (s(i1,i3)*s(i2,i4)*t(i1,i3,i4)) + 
     -  (8*(2*mhsq + s(i2,i4)))/(s(i2,i3)*t(i1,i3,i4)) - 
     -  (8*(-3*mhsq - 2*s(i1,i4) + s(i2,i3) - s(i3,i4)))/
     -   (s(i2,i4)*t(i1,i3,i4)) - 
     -  (8*(-3*mhsq + s(i2,i3) - s(i2,i4) - s(i3,i4)))/
     -   (s(i1,i4)*t(i1,i3,i4)) + 
     -  (8*mhsq*(-2*mhsq + s(i2,i4) + s(i3,i4)))/
     -   (s(i1,i4)*s(i2,i3)*t(i1,i3,i4)) + 
     -  (24*mhsq)/t(i2,i3,i4)**2 - 
     -  (8*mhsq*s(i2,i4))/(s(i2,i3)*t(i2,i3,i4)**2) + 
     -  (8*mhsq*s(i3,i4))/(s(i2,i4)*t(i2,i3,i4)**2) + 
     -  40/t(i2,i3,i4) + (8*(2*mhsq + s(i1,i4)))/
     -   (s(i1,i3)*t(i2,i3,i4)) + 
     -  (8*mhsq)/(s(i2,i3)*t(i2,i3,i4)) - 
     -  (8*mhsq*(2*mhsq - s(i1,i3) + s(i2,i4)))/
     -   (s(i1,i4)*s(i2,i3)*t(i2,i3,i4)) - 
     -  (8*(-3*mhsq + s(i1,i3) - s(i1,i4) - s(i3,i4)))/
     -   (s(i2,i4)*t(i2,i3,i4)) - 
     -  (8*(-3*mhsq + s(i1,i3) - 2*s(i2,i4) - s(i3,i4)))/
     -   (s(i1,i4)*t(i2,i3,i4)) + 
     -  (8*mhsq*(-2*mhsq + s(i1,i4) + s(i3,i4)))/
     -   (s(i1,i3)*s(i2,i4)*t(i2,i3,i4)) + 
     -  (16*mhsq**3)/
     -   (s(i1,i4)*s(i2,i3)*t(i1,i3,i4)*t(i2,i3,i4)) + 
     -  (8*(-4*mhsq**2 + 3*mhsq*s(i1,i4) - s(i1,i4)**2))/
     -   (s(i2,i3)*t(i1,i3,i4)*t(i2,i3,i4)) + 
     -  (16*mhsq**3)/
     -   (s(i1,i3)*s(i2,i4)*t(i1,i3,i4)*t(i2,i3,i4)) + 
     -  (8*(-4*mhsq**2 + 3*mhsq*s(i2,i4) - s(i2,i4)**2))/
     -   (s(i1,i3)*t(i1,i3,i4)*t(i2,i3,i4)) + 
     -  (16*(-6*mhsq - 4*s(i1,i4) - 4*s(i2,i4) - 3*s(i3,i4)))/
     -   (t(i1,i3,i4)*t(i2,i3,i4)) + 
     -  (8*(-4*mhsq**2 - 3*mhsq*s(i1,i4) - s(i1,i4)**2 - 
     -       3*mhsq*s(i3,i4) - 2*s(i1,i4)*s(i3,i4) - s(i3,i4)**2))
     -    /(s(i2,i4)*t(i1,i3,i4)*t(i2,i3,i4)) + 
     -  (8*(-4*mhsq**2 - 3*mhsq*s(i2,i4) - s(i2,i4)**2 - 
     -       3*mhsq*s(i3,i4) - 2*s(i2,i4)*s(i3,i4) - s(i3,i4)**2))
     -    /(s(i1,i4)*t(i1,i3,i4)*t(i2,i3,i4))
      return
      end
      
      function babisA0(i1,i2,i3,i4)
      implicit none
      include 'types.f' 
      real(dp):: babisA0
      integer i1,i2,i3,i4
      include 'constants.f'
      include 'mxpart.f'
      include 'qcdcouple.f'
      include 'sprods_com.f'
      include 'ewcouple.f' 
      include 'hbbparams.f'
      real(dp) :: t,mhsq,babis_NNLOfac
      
      t(i1,i2,i3)=s(i1,i2)+s(i2,i3)+s(i1,i3)
      mhsq=s(i1,i2)+s(i1,i3)+s(i1,i4)+s(i2,i3)+s(i2,i4)+s(i3,i4)

      babisA0= (8*mhsq)/t(i1,i3,i4)**2 + 32/s(i1,i3) + 24/s(i1,i4) + 
     -  32/s(i2,i3) + 24/s(i2,i4) - 
     -  (4*(8*mhsq - 5*s(i1,i4) - 5*s(i2,i3) - 7*s(i3,i4)))/
     -   (s(i1,i3)*s(i2,i4)) - 
     -  (4*(8*mhsq - 5*s(i1,i3) - 5*s(i2,i4) - 7*s(i3,i4)))/
     -   (s(i1,i4)*s(i2,i3)) + 
     -  (8*(-6*mhsq + 3*s(i1,i3) + 3*s(i2,i3) + 3*s(i3,i4)))/
     -   (s(i1,i4)*s(i2,i4)) + 
     -  (8*(-6*mhsq + 3*s(i1,i4) + 3*s(i2,i4) + 3*s(i3,i4)))/
     -   (s(i1,i3)*s(i2,i3)) + 
     -  (4*(-6*mhsq + 3*s(i1,i3) + 3*s(i1,i4) + 6*s(i3,i4)))/
     -   (s(i2,i3)*s(i2,i4)) + 
     -  (4*(-6*mhsq + 3*s(i2,i3) + 3*s(i2,i4) + 6*s(i3,i4)))/
     -   (s(i1,i3)*s(i1,i4)) + 
     -  (4*(mhsq - s(i3,i4))*
     -     (-2*mhsq**2 + 2*mhsq*s(i3,i4) - s(i3,i4)**2))/
     -   (s(i1,i3)*s(i1,i4)*s(i2,i3)*s(i2,i4)) + 
     -  (4*(4*mhsq**2 - 3*mhsq*s(i1,i3) + s(i1,i3)**2 - 
     -       6*mhsq*s(i3,i4) + 3*s(i1,i3)*s(i3,i4) + 3*s(i3,i4)**2
     -       ))/(s(i1,i4)*s(i2,i3)*s(i2,i4)) + 
     -  (4*(4*mhsq**2 - 3*mhsq*s(i1,i4) + s(i1,i4)**2 - 
     -       6*mhsq*s(i3,i4) + 3*s(i1,i4)*s(i3,i4) + 3*s(i3,i4)**2
     -       ))/(s(i1,i3)*s(i2,i3)*s(i2,i4)) + 
     -  (4*(4*mhsq**2 - 3*mhsq*s(i2,i3) + s(i2,i3)**2 - 
     -       6*mhsq*s(i3,i4) + 3*s(i2,i3)*s(i3,i4) + 3*s(i3,i4)**2
     -       ))/(s(i1,i3)*s(i1,i4)*s(i2,i4)) + 
     -  (4*(4*mhsq**2 - 3*mhsq*s(i2,i4) + s(i2,i4)**2 - 
     -       6*mhsq*s(i3,i4) + 3*s(i2,i4)*s(i3,i4) + 3*s(i3,i4)**2
     -       ))/(s(i1,i3)*s(i1,i4)*s(i2,i3)) + 20/t(i1,i3,i4) - 
     -  (4*(-6*mhsq - 3*s(i1,i4) + 3*s(i2,i3) - 3*s(i3,i4)))/
     -   (s(i2,i4)*t(i1,i3,i4)) - 
     -  (4*(-3*mhsq + s(i2,i3) - s(i2,i4) - 2*s(i3,i4)))/
     -   (s(i1,i4)*t(i1,i3,i4)) + 
     -  (4*(4*mhsq - s(i1,i4) - s(i2,i4) + s(i3,i4)))/
     -   (s(i2,i3)*t(i1,i3,i4)) - 
     -  (4*(-4*mhsq**2 + 3*mhsq*s(i2,i3) - s(i2,i3)**2 - 
     -       3*mhsq*s(i3,i4) + s(i2,i3)*s(i3,i4) - s(i3,i4)**2))/
     -   (s(i1,i4)*s(i2,i4)*t(i1,i3,i4)) - 
     -  (4*(-4*mhsq**2 + 3*mhsq*s(i2,i4) - s(i2,i4)**2 - 
     -       3*mhsq*s(i3,i4) + s(i2,i4)*s(i3,i4) - s(i3,i4)**2))/
     -   (s(i1,i3)*s(i2,i3)*t(i1,i3,i4)) + 
     -  (8*mhsq)/t(i2,i3,i4)**2 + 20/t(i2,i3,i4) - 
     -  (4*(-3*mhsq + s(i1,i3) - s(i1,i4) - 2*s(i3,i4)))/
     -   (s(i2,i4)*t(i2,i3,i4)) + 
     -  (4*(4*mhsq - s(i1,i4) - s(i2,i4) + s(i3,i4)))/
     -   (s(i1,i3)*t(i2,i3,i4)) + 
     -  (4*(6*mhsq - 3*s(i1,i3) + 3*s(i2,i4) + 3*s(i3,i4)))/
     -   (s(i1,i4)*t(i2,i3,i4)) - 
     -  (4*(-4*mhsq**2 + 3*mhsq*s(i1,i3) - s(i1,i3)**2 - 
     -       3*mhsq*s(i3,i4) + s(i1,i3)*s(i3,i4) - s(i3,i4)**2))/
     -   (s(i1,i4)*s(i2,i4)*t(i2,i3,i4)) - 
     -  (4*(-4*mhsq**2 + 3*mhsq*s(i1,i4) - s(i1,i4)**2 - 
     -       3*mhsq*s(i3,i4) + s(i1,i4)*s(i3,i4) - s(i3,i4)**2))/
     -   (s(i1,i3)*s(i2,i3)*t(i2,i3,i4)) + 
     -  (8*(-8*mhsq - 4*s(i1,i4) - 4*s(i2,i4) - 5*s(i3,i4)))/
     -   (t(i1,i3,i4)*t(i2,i3,i4)) + 
     -  (4*(-4*mhsq**2 - 3*mhsq*s(i1,i4) - s(i1,i4)**2 - 
     -       6*mhsq*s(i3,i4) - 3*s(i1,i4)*s(i3,i4) - 3*s(i3,i4)**2
     -       ))/(s(i2,i4)*t(i1,i3,i4)*t(i2,i3,i4)) + 
     -  (4*(-4*mhsq**2 - 3*mhsq*s(i2,i4) - s(i2,i4)**2 - 
     -       6*mhsq*s(i3,i4) - 3*s(i2,i4)*s(i3,i4) - 3*s(i3,i4)**2
     -       ))/(s(i1,i4)*t(i1,i3,i4)*t(i2,i3,i4)) + 
     -  (4*(mhsq + s(i3,i4))*
     -     (-2*mhsq**2 - 2*mhsq*s(i3,i4) - s(i3,i4)**2))/
     -   (s(i1,i3)*s(i2,i3)*t(i1,i3,i4)*t(i2,i3,i4)) + 
     -  (4*(mhsq + s(i3,i4))*
     -     (-2*mhsq**2 - 2*mhsq*s(i3,i4) - s(i3,i4)**2))/
     -   (s(i1,i4)*s(i2,i4)*t(i1,i3,i4)*t(i2,i3,i4)) + 
     -  (4*(-4*mhsq**2 + 3*mhsq*s(i1,i4) - s(i1,i4)**2 - 
     -       3*mhsq*s(i3,i4) + s(i1,i4)*s(i3,i4) - s(i3,i4)**2))/
     -   (s(i2,i3)*t(i1,i3,i4)*t(i2,i3,i4)) + 
     -  (4*(-4*mhsq**2 + 3*mhsq*s(i2,i4) - s(i2,i4)**2 - 
     -       3*mhsq*s(i3,i4) + s(i2,i4)*s(i3,i4) - s(i3,i4)**2))/
     -     (s(i1,i3)*t(i1,i3,i4)*t(i2,i3,i4))

      return
      end
