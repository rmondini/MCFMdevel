!---- real amplitudes for Hbgg (i.e. Hbbgg and HbbQQ) 
      subroutine Hbbg_realamps(i1,i2,i3,i4,za,zb,ampsgg,ampsQQ)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f' 
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer i1,i2,i3,i4
      integer icol
      complex(dp) :: ampsgg(0:2,2,2,2),ampsQQ(1:2,2,2)
      complex(dp) amps(2,2),s1234
      complex(dp) Hbbgg_allm,Hbbgg_mmpp,Hbbgg_mmmp,HbbQQ_mmmp
      integer h1,h2,h3,h4
      integer iint
           
      s1234=s(i1,i2)+s(i1,i3)+s(i1,i4)+s(i2,i3)+s(i2,i4)+s(i3,i4)

      icol=1
!--- icol=1, i3,i4 ordering
      ampsgg(icol,1,1,1)=Hbbgg_allm(i1,i2,i3,i4,za,zb)
      ampsgg(icol,1,1,2)=-Hbbgg_mmmp(i1,i2,i3,i4,za,zb)
      ampsgg(icol,1,2,2)=-Hbbgg_mmpp(i1,i2,i3,i4,za,zb)
      ampsgg(icol,1,2,1)=Hbbgg_mmmp(i2,i1,i4,i3,za,zb)

      ampsgg(icol,2,2,2)=Hbbgg_allm(i1,i2,i3,i4,zb,za)
      ampsgg(icol,2,2,1)=-Hbbgg_mmmp(i1,i2,i3,i4,zb,za)
      ampsgg(icol,2,1,1)=-Hbbgg_mmpp(i1,i2,i3,i4,zb,za)
      ampsgg(icol,2,1,2)=Hbbgg_mmmp(i2,i1,i4,i3,zb,za)
 
    
      
      icol=2
!---- icol=2 i4,i3 ordering
      ampsgg(icol,1,1,1)=Hbbgg_allm(i1,i2,i4,i3,za,zb)
      ampsgg(icol,1,2,1)=-Hbbgg_mmmp(i1,i2,i4,i3,za,zb)
      ampsgg(icol,1,2,2)=-Hbbgg_mmpp(i1,i2,i4,i3,za,zb)
      ampsgg(icol,1,1,2)=Hbbgg_mmmp(i2,i1,i3,i4,za,zb)

      ampsgg(icol,2,2,2)=Hbbgg_allm(i1,i2,i4,i3,zb,za)
      ampsgg(icol,2,1,2)=-Hbbgg_mmmp(i1,i2,i4,i3,zb,za)
      ampsgg(icol,2,1,1)=-Hbbgg_mmpp(i1,i2,i4,i3,zb,za)
      ampsgg(icol,2,2,1)=Hbbgg_mmmp(i2,i1,i3,i4,zb,za)

      
      icol=0
      ampsgg(icol,:,:,:)=ampsgg(1,:,:,:)+ampsgg(2,:,:,:)

      iint=1
!----- int =1, regular ordering, int=2  i2 <-> i4 for interference 

!---- i1,i2 i3,i4 QQ pair (i1,i2 couple to Higgs) 
      
      ampsQQ(iint,1,1)=HbbQQ_mmmp(i1,i2,i3,i4,za,zb)
      ampsQQ(iint,2,2)=HbbQQ_mmmp(i1,i2,i3,i4,zb,za)
      ampsQQ(iint,1,2)=HbbQQ_mmmp(i1,i2,i4,i3,za,zb)
      ampsQQ(iint,2,1)=HbbQQ_mmmp(i1,i2,i4,i3,zb,za)

      iint=2
      ampsQQ(iint,1,1)=HbbQQ_mmmp(i1,i4,i3,i2,za,zb)
      ampsQQ(iint,2,2)=HbbQQ_mmmp(i1,i4,i3,i2,zb,za)
      ampsQQ(iint,1,2)=HbbQQ_mmmp(i1,i4,i2,i3,za,zb)
      ampsQQ(iint,2,1)=HbbQQ_mmmp(i1,i4,i2,i3,zb,za)

      
!----- debug, check amplitudes
      
!      do h1=1,2
!         do h2=1,2
!            do h3=1,2
!               write(6,*) h1,h2,h3,ampsgg(1,h1,h2,h3),ampsgg(2,h1,h2,h3)
!            enddo
!         enddo
!      enddo
!      stop
      return
      end

      function HbbQQ_mmmp(i1,i2,i3,i4,za,zb)
      implicit none
      include 'types.f'
      complex(dp) :: HbbQQ_mmmp 
      include 'constants.f' 
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer i1,i2,i3,i4
      complex(dp) :: zab2
      real(dp) :: t
      
      zab2(i1,i2,i3,i4) = za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)
      t(i1,i2,i3)=s(i1,i2)+s(i1,i3)+s(i2,i3)

      HbbQQ_mmmp= za(i1,i3)*zab2(i2,i1,i3,i4)/(t(i1,i3,i4)*s(i3,i4))
     &     +za(i2,i3)*zab2(i1,i2,i3,i4)/(t(i2,i3,i4)*s(i3,i4))
      
      
      return
      end

      function Hbbgg_allm(i1,i2,i3,i4,za,zb)
      implicit none
      include 'types.f'
      complex(dp) :: Hbbgg_allm 
      include 'constants.f' 
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer i1,i2,i3,i4
      real(dp) s1234
      s1234=s(i1,i2)+s(i1,i3)+s(i1,i4)+s(i2,i3)+s(i2,i4)+s(i3,i4)

      Hbbgg_allm = s1234/(zb(i1,i4)*zb(i2,i3)*zb(i3,i4))
      
      return
      end

      
      function Hbbgg_mmpp(i1,i2,i3,i4,za,zb)
      implicit none
      include 'types.f'
      complex(dp) :: Hbbgg_mmpp 
      include 'constants.f' 
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer i1,i2,i3,i4
   
      Hbbgg_mmpp = za(i1,i2)**2/(za(i1,i4)*za(i2,i3)*za(i3,i4))
      return
      end

      function Hbbgg_mmmp(i1,i2,i3,i4,za,zb)
      implicit none
      include 'types.f'
      complex(dp) :: Hbbgg_mmmp 
      include 'constants.f' 
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer i1,i2,i3,i4
      complex(dp) :: zab2
      real(dp) t
      zab2(i1,i2,i3,i4) = za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)
      t(i1,i2,i3)=s(i1,i2)+s(i2,i3)+s(i1,i3)
      
      Hbbgg_mmmp = za(i2,i3)*zb(i4,i2)*zab2(i1,i2,i3,i4)
     &     /(t(i2,i3,i4)*za(i3,i4)*zb(i2,i3)*zb(i3,i4))

      Hbbgg_mmmp=Hbbgg_mmmp+za(i1,i3)**2*zab2(i2,i1,i3,i4)
     &     /(t(i1,i3,i4)*za(i1,i4)*za(i3,i4)*zb(i3,i4))

      Hbbgg_mmmp=Hbbgg_mmmp-za(i1,i3)*zab2(i1,i2,i3,i4)
     &     /(za(i1,i4)*za(i3,i4)*zb(i2,i3)*zb(i3,i4))
      
      return
      end
