!     =--- nf loops -+-+
      
      function amp_qqbQQbga_del1_nf(i1,i2,i3,i4,i5,za,zb)
      implicit none
      include 'types.f'
      complex(dp) :: amp_qqbQQbga_del1_nf
      include 'constants.f'
      include 'mxpart.f' 
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'scale.f' 
      integer i1,i2,i3,i4,i5
      integer Nbox,Ntri,Nbub
      parameter(Nbox=5,Ntri=6,Nbub=5)
      complex(dp) ::  d0(Nbox),c0(Ntri),b0(Nbub),rat
      complex(dp) :: ALO,Vpole
      complex(dp) :: L0,L1,Bubs,Boxes,Lsm1,lnrat
      complex(dp) :: Boxint(4),boxcoeff(4)
      integer i 

      ALO=-(za(i1,i3)**2/(za(i1,i5)*za(i2,i5)*za(i3,i4)))

!      write(6,*) ALO 
!---- total amplitude is pole piece

      Vpole = (2*epinv) + (2*(2 + lnrat(musq,-s(i3,i4)))
     &     -2._dp/3._dp)

      amp_qqbQQbga_del1_nf=Vpole*ALO

!      write(6,*) 'dell 1 -+-+'
!      write(6,*) 'LO',ALO
!      write(6,*) 'total = ',amp_qqbQQbga_del1_nf*im
      return
      end
     

      function amp_qqbQQbga_del2_nf(i1,i2,i3,i4,i5,za,zb)
      implicit none
      include 'types.f'
      complex(dp) :: amp_qqbQQbga_del2_nf
      complex(dp) :: amp_qqbQQbga_del1_nf
      include 'constants.f'
      include 'mxpart.f' 
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'scale.f' 
      integer i1,i2,i3,i4,i5
      integer Nbox,Ntri,Nbub
      parameter(Nbox=5,Ntri=6,Nbub=5)
      complex(dp) ::  d0(Nbox),c0(Ntri),b0(Nbub),rat
      complex(dp) :: ALO,Vpole
      complex(dp) :: L0,L1,Bubs,Boxes,Lsm1,lnrat
      complex(dp) :: Boxint(4),boxcoeff(4)
      integer i 

      ALO=-(za(i1,i3)**2/(za(i1,i5)*za(i2,i5)*za(i3,i4)))
    
!---- total amplitude is pole piece

      Vpole = (2*epinv) + (2*(2 + lnrat(musq,-s(i3,i4)))
     &     -2._dp/3._dp)
      amp_qqbQQbga_del2_nf=-Vpole*ALO

!      write(6,*) 'dell 2 -+-+'
!      write(6,*) 'LO',ALO
!      write(6,*) 'total = ',amp_qqbQQbga_del2_nf*im
!      write(6,*)
      return
      end
!=====nf loops -++-

      function amp_qqbQQbga_mhvalt_del1_nf(i1,i2,i3,i4,i5,za,zb)
      implicit none
      include 'types.f'
      complex(dp) :: amp_qqbQQbga_mhvalt_del1_nf
      include 'constants.f'
      include 'mxpart.f' 
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'scale.f' 
      integer i1,i2,i3,i4,i5
      integer Nbox,Ntri,Nbub
      parameter(Nbox=5,Ntri=6,Nbub=5)
      complex(dp) ::  d0(Nbox),c0(Ntri),b0(Nbub),rat
      complex(dp) :: ALO,Vpole
      complex(dp) :: L0,L1,Bubs,Boxes,Lsm1,lnrat
      complex(dp) :: Boxint(4),boxcoeff(4)
      integer i 

      ALO=za(i1,i4)**2/(za(i1,i5)*za(i2,i5)*za(i3,i4))

!---- total amplitude is pole piece

      Vpole = (2*epinv) + (2*(2 + lnrat(musq,-s(i3,i4)))
     &     -2._dp/3._dp)
      
      amp_qqbQQbga_mhvalt_del1_nf=Vpole*ALO

!      write(6,*) 'dell 1 --++'
 
!      write(6,*) 'LO',ALO
!      write(6,*) 'total = ',amp_qqbQQbga_mhvalt_del1_nf*im
!      write(6,*)
      return
      end
     

      function amp_qqbQQbga_mhvalt_del2_nf(i1,i2,i3,i4,i5,za,zb)
      implicit none
      include 'types.f'
      complex(dp) :: amp_qqbQQbga_mhvalt_del2_nf
      complex(dp) :: amp_qqbQQbga_mhvalt_del1_nf
      include 'constants.f'
      include 'mxpart.f' 
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'scale.f' 
      integer i1,i2,i3,i4,i5
      integer Nbox,Ntri,Nbub
      parameter(Nbox=5,Ntri=6,Nbub=5)
      complex(dp) ::  d0(Nbox),c0(Ntri),b0(Nbub),rat
      complex(dp) :: ALO,Vpole
      complex(dp) :: L0,L1,Bubs,Boxes,Lsm1,lnrat
      complex(dp) :: Boxint(4),boxcoeff(4)
      integer i 


      ALO=za(i1,i4)**2/(za(i1,i5)*za(i2,i5)*za(i3,i4))

!---- total amplitude is pole piece

      Vpole = (2*epinv) + (2*(2 + lnrat(musq,-s(i3,i4)))
     &     -2._dp/3._dp)
        
      amp_qqbQQbga_mhvalt_del2_nf=-Vpole*ALO
!      write(6,*) 'dell 2 -++-'
! 
!      write(6,*) 'LO',ALO
!      write(6,*) 'total = ',amp_qqbQQbga_mhvalt_del2_nf*im
!      write(6,*)
      return
      end
     
