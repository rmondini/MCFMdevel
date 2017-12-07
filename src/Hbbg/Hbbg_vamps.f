!---- lc amplitudes
!----- amplitudes are for H=>(b(i1)+bb(i2)+g(i3))
      subroutine Hbbg_vamp_lc(i1,i2,i3,za,zb,vamps)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f' 
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'epinv.f'
      include 'scale.f' 
      integer i1,i2,i3
      complex(dp) vamps(2,2),loamps(2,2),s123
      complex(dp) lsm1,vpole,rat,box,ratc,lnrat
      complex(dp) s12,s13,s23
      
      s123=s(i1,i2)+s(i2,i3)+s(i1,i3)
      s12=s(i1,i2)
      s23=s(i2,i3)
      s13=s(i1,i3)

      
      call Hbbg_loamp(i1,i2,i3,za,zb,loamps)
!---- vpole is helicity universal
      Vpole= 2*epinv**2 + lnrat(musq,-s13)**2/2. +
     -  epinv*(+lnrat(musq,-s13) + lnrat(musq,-s23)) +
     -  lnrat(musq,-s23)**2/2._dp
      box=Lsm1(-s13,-s123,-s23,-s123)

      vamps(:,:)=loamps(:,:)*(Vpole+box)

!---- rational terms for diagonal helicity
!--- this is rat for all minus amplitude 
      rat=-(s13+s23)/(2._dp*zb(i2,i3)*zb(i1,i3))
!---- rat for all plus amplitude
      ratc = (s13+s23)/(2._dp*za(i2,i3)*za(i1,i3))

!---- add in rational term
      vamps(1,1)=vamps(1,1)+rat*1._dp
      vamps(2,2)=vamps(2,2)+ratc*1._dp
      
      return
      end

!---- subleading amplitudes

      
      subroutine Hbbg_vamp_slc(i1,i2,i3,za,zb,vamps)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f' 
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'epinv.f'
      include 'scale.f' 
      integer i1,i2,i3
      complex(dp) vamps(2,2),loamps(2,2),s123
      complex(dp) lsm1,vpole,rat,box1,box2,ratc,lnrat
      complex(dp) s12,s13,s23
      
      s123=s(i1,i2)+s(i2,i3)+s(i1,i3)
      s12=s(i1,i2)
      s23=s(i2,i3)
      s13=s(i1,i3)

      
      call Hbbg_loamp(i1,i2,i3,za,zb,loamps)
!---- vpole is helicity universal
      Vpole=  -epinv**2 - epinv*lnrat(musq,-s12) - 
     -  lnrat(musq,-s12)**2/2._dp


!----- box function 
      box1=-Lsm1(-s13,-s123,-s12,-s123)
      box2=-Lsm1(-s23,-s123,-s12,-s123)

      vamps(:,:)=loamps(:,:)*(Vpole+box1+box2)

!---- rational terms for diagonal helicity
!--- this is rat for all minus amplitude 
      rat=-(s13+s23)/(2._dp*zb(i2,i3)*zb(i1,i3))
!---- rat for all plus amplitude
      ratc = (s13+s23)/(2._dp*za(i2,i3)*za(i1,i3))
    
!---- add in rational term
      vamps(1,1)=vamps(1,1)+rat*1._dp
      vamps(2,2)=vamps(2,2)+ratc*1._dp
      
      return
      end
