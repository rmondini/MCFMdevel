!---  CW May 16
!---- collection of QQB qqb ga amplitudes
!==== del1 is leading color delta combination
!==== del2 is sub leading color combination (i.e. del-1/N_c del2)
!==== at NLO both have detailed color structure

!------ this is the -++- configuration 
      
!------photon couples to i1,i2 line (qqb)
      function amp_qqbQQbga_mhvalt_del1_lc(i1,i2,i3,i4,i5,za,zb)
      implicit none
      include 'types.f'
      complex(dp) :: amp_qqbQQbga_mhvalt_del1_lc
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

      
      Vpole= 2*epinv*epinv2 + 
     -  epinv*(-2._dp/3._dp + 
     -     lnrat(musq,-s(i1,i4)) + 
     -     lnrat(musq,-s(i2,i3))) + 
     -  (-8 + 3*lnrat(musq,-s(i1,i4))**2 + 
     -     3*lnrat(musq,-s(i2,i3))**2 - 
     -     4*lnrat(musq,-s(i3,i4)))/6.

      vpole=ALO*vpole

      boxint(1)=Lsm1(-s(i2,i3),-s(i1,i5),-s(i3,i4),-s(i1,i5))
      boxint(2)=Lsm1(-s(i1,i4),-s(i2,i5),-s(i3,i4),-s(i2,i5))
      boxint(3)=Lsm1(-s(i2,i5),-s(i1,i4),-s(i2,i3),-s(i1,i4))
      boxint(4)=Lsm1(-s(i1,i4),-s(i2,i3),-s(i1,i5),-s(i2,i3))

      boxcoeff(1)=ALO
      boxcoeff(2)=ALO
      boxcoeff(3)=ALO
      boxcoeff(4)=ALO

      boxes=czip
      do i=1,4
         boxes=boxes+boxcoeff(i)*boxint(i)
      enddo

      bubs=(L0(-s(i1,i5),-s(i3,i4))*za(i1,i4)*za(i2,i4)*zb(i5,i2))/
     -   (s(i3,i4)*za(i2,i5)*za(i3,i4)) + 
     -  (L1(-s(i1,i5),-s(i3,i4))*za(i1,i5)*za(i2,i4)**2*zb(i5,i2)**2)/
     -     (2.*s(i3,i4)**2*za(i2,i5)*za(i3,i4))

      Rat= -((za(i1,i2)*za(i1,i4)*zb(i3,i2))/
     -     (za(i1,i5)*za(i2,i5)*za(i3,i4)*zb(i4,i3))) - 
     -  (za(i1,i4)**2*zb(i5,i1))/
     -     (2.*za(i2,i5)*za(i3,i4)**2*zb(i4,i3))
!======missing bubble piece
      rat=rat+2._dp/9._dp*ALO

      amp_qqbQQbga_mhvalt_del1_lc=vpole+boxes-bubs+Rat

 !     write(6,*) 'LO mv =',LO
 !     write(6,*) 'total ', amp_qqbQQbga_mhvalt_del1_lc*im
 !     write(6,*) 'pole pieces',Vpole*im
 !     write(6,*) 'boxes ',Boxes*im
 !     write(6,*) 'bubs ',Bubs*im
 !     write(6,*) 'completed rat',rat*im
 !     write(6,*) 'total/LO',amp_qqbQQbga_mhvalt_del1_lc*im/ALO
      return
      end


      function amp_qqbQQbga_mhvalt_del1_slc(i1,i2,i3,i4,i5,za,zb)
      implicit none
      include 'types.f'
      complex(dp) :: amp_qqbQQbga_mhvalt_del1_slc
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
      complex(dp) :: Boxint(10),boxcoeff(10)
      integer i 

      ALO=za(i1,i4)**2/(za(i1,i5)*za(i2,i5)*za(i3,i4))

      Vpole= -6 - 2*epinv*epinv2 - 
     -  lnrat(musq,-s(i1,i2))**2/2. + 
     -  lnrat(musq,-s(i1,i3))**2 - 
     -  lnrat(musq,-s(i1,i4))**2 - 
     -  lnrat(musq,-s(i2,i3))**2 + 
     -  lnrat(musq,-s(i2,i4))**2 + 
     -  epinv*(-3 - lnrat(musq,-s(i1,i2)) + 
     -     2*lnrat(musq,-s(i1,i3)) - 
     -     2*lnrat(musq,-s(i1,i4)) - 
     -     2*lnrat(musq,-s(i2,i3)) + 
     -     2*lnrat(musq,-s(i2,i4)) - 
     -     lnrat(musq,-s(i3,i4))) - 
     -  3*lnrat(musq,-s(i3,i4)) - 
     -  lnrat(musq,-s(i3,i4))**2/2.
            
      Vpole=Vpole*ALO

      boxint(1)=Lsm1(-s(i1,i5),-s(i3,i4),-s(i1,i2),-s(i3,i4))
      boxint(2)=Lsm1(-s(i2,i5),-s(i3,i4),-s(i1,i2),-s(i3,i4))
      boxint(3)=Lsm1(-s(i2,i3),-s(i1,i5),-s(i3,i4),-s(i1,i5))
      boxint(4)=Lsm1(-s(i2,i4),-s(i1,i5),-s(i3,i4),-s(i1,i5))
      boxint(5)=Lsm1(-s(i1,i3),-s(i2,i5),-s(i3,i4),-s(i2,i5))
      boxint(6)=Lsm1(-s(i1,i4),-s(i2,i5),-s(i3,i4),-s(i2,i5))
      boxint(7)=Lsm1(-s(i1,i5),-s(i2,i4),-s(i1,i3),-s(i2,i4))
      boxint(8)=Lsm1(-s(i2,i4),-s(i1,i3),-s(i2,i5),-s(i1,i3))
      boxint(9)=Lsm1(-s(i2,i5),-s(i1,i4),-s(i2,i3),-s(i1,i4))
      boxint(10)=Lsm1(-s(i1,i4),-s(i2,i3),-s(i1,i5),-s(i2,i3))

!------ box coefficients 
      boxcoeff(1)=  -(za(i1,i2)**2*za(i4,i5)**2)/
     -  (za(i1,i5)*za(i2,i5)**3*za(i3,i4))
      boxcoeff(2)=-ALO
      boxcoeff(3)= -2.*ALO
      boxcoeff(4)=(2*za(i1,i3)**2*za(i2,i4)**2)/
     -  (za(i1,i5)*za(i2,i3)**2*za(i2,i5)*
     -     za(i3,i4))
      boxcoeff(5)=2.*ALO
      boxcoeff(6)=-2.*ALO
      boxcoeff(7)=(2*za(i1,i3)**2*za(i4,i5)**2)/
     -  (za(i1,i5)*za(i2,i5)*za(i3,i4)*
     -     za(i3,i5)**2)
      boxcoeff(8)=2*ALO
      boxcoeff(9) = -2*ALO
      boxcoeff(10) = -2.*ALO


!      write(6,*) boxcoeff(1)/2._dp*(s(i1,i5)*s(i1,i2))*im
!      write(6,*) boxcoeff(2)/2._dp*(s(i2,i5)*s(i1,i2))*im
!      write(6,*) boxcoeff(3)/2._dp*(s(i3,i4)*s(i2,i3))*im
!      write(6,*) boxcoeff(4)/2._dp*(s(i3,i4)*s(i2,i4))*im
!      write(6,*) boxcoeff(5)/2._dp*(s(i3,i4)*s(i1,i3))*im
!      write(6,*) boxcoeff(6)/2._dp*(s(i3,i4)*s(i1,i4))*im
!      write(6,*) boxcoeff(7)/2._dp*(s(i1,i5)*s(i1,i3))*im
!      write(6,*) boxcoeff(8)/2._dp*(s(i2,i5)*s(i2,i4))*im
!      write(6,*) boxcoeff(9)/2._dp*(s(i2,i5)*s(i2,i3))*im
!      write(6,*) boxcoeff(10)/2._dp*(s(i1,i5)*s(i1,i4))*im

      Boxes=czip
      do i=1,10
         Boxes=Boxes+boxint(i)*boxcoeff(i)
      enddo
      
      bubs=(-2*lnrat(-s(i1,i5),-s(i3,i4))*za(i1,i3)*za(i1,i4))/
     -   (za(i1,i5)*za(i2,i3)*za(i3,i5)) + 
     -  (2*lnrat(-s(i1,i3),-s(i3,i4))*za(i1,i4))/
     -   (za(i2,i5)*za(i3,i5)) - 
     -  (2*L0(-s(i2,i4),-s(i1,i5))*za(i1,i3)**2*za(i2,i4)*
     -     zb(i3,i2))/(s(i1,i5)*za(i1,i5)*za(i2,i3)*za(i3,i5)) + 
     -  (2*L0(-s(i2,i4),-s(i1,i3))*za(i1,i5)*za(i2,i4)*zb(i5,i2))/
     -   (s(i1,i3)*za(i2,i5)*za(i3,i5)) - 
     -  (2*L0(-s(i1,i2),-s(i3,i4))*za(i1,i2)*za(i2,i4)*za(i4,i5)*
     -     zb(i5,i2))/(s(i3,i4)*za(i2,i5)**2*za(i3,i4)) - 
     -  (L1(-s(i1,i5),-s(i3,i4))*za(i1,i2)*za(i2,i4)*zb(i3,i2)*
     -     zb(i5,i2))/(2.*s(i3,i4)**2*za(i2,i5)) - 
     -  (L1(-s(i1,i2),-s(i3,i4))*za(i1,i2)**2*za(i4,i5)**2*
     -     zb(i5,i1)*zb(i5,i2))/
     -   (s(i3,i4)**2*za(i2,i5)**2*za(i3,i4)) - 
     -  (L0(-s(i1,i5),-s(i3,i4))*
     -     ((-2*za(i1,i2)*za(i2,i4)*zb(i5,i2))/
     -        (za(i2,i3)*za(i2,i5)) + 
     -       (za(i1,i4)*za(i2,i4)*zb(i5,i2))/
     -        (2.*za(i2,i5)*za(i3,i4)) + 
     -       (za(i1,i2)*za(i2,i4)*za(i4,i5)*zb(i5,i2))/
     -     (za(i2,i5)**2*za(i3,i4))))/s(i3,i4)


      rat= -(za(i1,i4)**2/(za(i1,i5)*za(i2,i5)*za(i3,i4))) + 
     -  (za(i1,i2)*za(i2,i4)*za(i4,i5)*zb(i5,i2))/
     -   (za(i2,i5)**2*za(i3,i4)**2*zb(i4,i3)) + 
     -  (za(i1,i4)*zb(i5,i3))/
     -   (2.*za(i2,i5)*za(i3,i4)*zb(i4,i3))
      amp_qqbQQbga_mhvalt_del1_slc=Vpole+boxes-bubs+rat

!      write(6,*) 'total ', amp_qqbQQbga_mhvalt_del1_slc*im
!      write(6,*) 'pole pieces',Vpole*im
!      write(6,*) 'boxes ',Boxes*im
!      write(6,*) 'bubs ',Bubs*im
!      write(6,*) 'completed rat',rat*im
!      write(6,*) 'total/LO',amp_qqbQQbga_mhvalt_del1_slc*im/ALO


      return
      end


!=====SECOND DELTA PIECS LEADING COLOR

      
      function amp_qqbQQbga_mhvalt_del2_lc(i1,i2,i3,i4,i5,za,zb)
      implicit none
      include 'types.f'
      complex(dp) :: amp_qqbQQbga_mhvalt_del2_lc
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
      vpole=   -2*epinv*epinv2 + 
     -  epinv*(2._dp/3._dp - 
     -     lnrat(musq,-s(i1,i3)) - 
     -     lnrat(musq,-s(i2,i4))) + 
     -  (8 - 3*lnrat(musq,-s(i1,i3))**2 - 
     -     3*lnrat(musq,-s(i2,i4))**2 + 
     -     4*lnrat(musq,-s(i3,i4)))/6.

      vpole=ALO*vpole
      
      boxint(1)=Lsm1(-s(i2,i4),-s(i1,i5),-s(i3,i4),-s(i1,i5))
      boxint(2)=Lsm1(-s(i1,i3),-s(i2,i5),-s(i3,i4),-s(i2,i5))
      boxint(3)=Lsm1(-s(i1,i5),-s(i2,i4),-s(i1,i3),-s(i2,i4))
      boxint(4)=Lsm1(-s(i2,i4),-s(i1,i3),-s(i2,i5),-s(i1,i3))

      boxcoeff(1)=-((za(i1,i3)**2*za(i2,i4)**2)/
     -     (za(i1,i5)*za(i2,i3)**2*za(i2,i5)*za(i3,i4)))
      boxcoeff(2)=-ALO
      boxcoeff(3)=-((za(i1,i3)**2*za(i4,i5)**2)/
     -     (za(i1,i5)*za(i2,i5)*za(i3,i4)*za(i3,i5)**2))
      boxcoeff(4)=-ALO

!      write(6,*) boxcoeff(1)/2._dp*(s(i2,i4)*s(i3,i4))*im
!      write(6,*) boxcoeff(2)/2._dp*(s(i1,i3)*s(i3,i4))*im
!      write(6,*) boxcoeff(3)/2._dp*(s(i1,i5)*s(i1,i3))*im
!      write(6,*) boxcoeff(4)/2._dp*(s(i2,i5)*s(i2,i4))*im

      
      boxes=czip
      do i=1,4
         boxes=boxes+boxcoeff(i)*boxint(i)
      enddo

      bubs= -((lnrat(-s(i2,i4),-s(i3,i4))*za(i1,i4)*za(i1,i3))/
     -     (za(i1,i5)*za(i2,i3)*za(i3,i5))) + 
     -  (lnrat(-s(i2,i4),-s(i3,i4))*za(i1,i4))/(za(i2,i5)*za(i3,i5)) - 
     -  (L0(-s(i1,i5),-s(i3,i4))*za(i1,i3)*za(i2,i4)**2*zb(i5,i2))/
     -   (s(i3,i4)*za(i2,i3)*za(i2,i5)*za(i4,i3)) - 
     -  (L1(-s(i1,i5),-s(i3,i4))*za(i1,i5)*za(i2,i4)**2*zb(i5,i2)**2)/
     -   (2.*s(i3,i4)**2*za(i2,i5)*za(i4,i3)) + 
     -  (L0(-s(i1,i5),-s(i2,i4))*za(i1,i3)*za(i4,i3)*zb(i5,i3))/
     -   (s(i2,i4)*za(i2,i3)*za(i3,i5)) + 
     -  (L0(-s(i1,i3),-s(i2,i4))*za(i1,i3)*za(i4,i5)*zb(i5,i3))/
     -     (s(i2,i4)*za(i2,i5)*za(i3,i5))

      
      Rat=    -((za(i1,i2)*za(i1,i4)*zb(i3,i2))/
     -     (za(i1,i5)*za(i2,i5)*za(i3,i4)*
     -       zb(i4,i3))) - 
     -  (za(i1,i4)**2*zb(i5,i1))/
     -     (2.*za(i2,i5)*za(i3,i4)**2*zb(i4,i3))

!======missing bubble piece
      rat=rat+2._dp/9._dp*ALO

      amp_qqbQQbga_mhvalt_del2_lc=vpole+boxes+bubs-Rat
    
!      write(6,*) 'total ', amp_qqbQQbga_mhvalt_del2_lc*im
!      write(6,*) 'pole pieces',Vpole*im
!      write(6,*) 'boxes ',Boxes*im
!      write(6,*) 'bubs ',Bubs*im
!      write(6,*) 'completed rat',rat*im
!      write(6,*) 'total/LO',amp_qqbQQbga_mhvalt_del2_lc*im/ALO

      return
      end


      
!=====SECOND DELTA PIECS SUBLEADING COLOR

      
      function amp_qqbQQbga_mhvalt_del2_slc(i1,i2,i3,i4,i5,za,zb)
      implicit none
      include 'types.f'
      complex(dp) :: amp_qqbQQbga_mhvalt_del2_slc
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
      complex(dp) :: Boxint(10),boxcoeff(10)
      integer i 
      
      ALO=za(i1,i4)**2/(za(i1,i5)*za(i2,i5)*za(i3,i4))
      vpole= 2*epinv*epinv2 + 
     -  epinv*(3 + lnrat(musq,-s(i1,i2)) - 
     -     lnrat(musq,-s(i1,i3)) + 
     -     lnrat(musq,-s(i1,i4)) + 
     -     lnrat(musq,-s(i2,i3)) - 
     -     lnrat(musq,-s(i2,i4)) + 
     -     lnrat(musq,-s(i3,i4))) + 
     -  (12 + lnrat(musq,-s(i1,i2))**2 - 
     -     lnrat(musq,-s(i1,i3))**2 + 
     -     lnrat(musq,-s(i1,i4))**2 + 
     -     lnrat(musq,-s(i2,i3))**2 - 
     -     lnrat(musq,-s(i2,i4))**2 + 
     -     6*lnrat(musq,-s(i3,i4)) + 
     -     lnrat(musq,-s(i3,i4))**2)/2.

      vpole=vpole*ALO

      
           boxint(1)=Lsm1(-s(i1,i5),-s(i3,i4),-s(i1,i2),-s(i3,i4))
      boxint(2)=Lsm1(-s(i2,i5),-s(i3,i4),-s(i1,i2),-s(i3,i4))
      boxint(3)=Lsm1(-s(i2,i3),-s(i1,i5),-s(i3,i4),-s(i1,i5))
      boxint(4)=Lsm1(-s(i2,i4),-s(i1,i5),-s(i3,i4),-s(i1,i5))
      boxint(5)=Lsm1(-s(i1,i3),-s(i2,i5),-s(i3,i4),-s(i2,i5))
      boxint(6)=Lsm1(-s(i1,i4),-s(i2,i5),-s(i3,i4),-s(i2,i5))
      boxint(7)=Lsm1(-s(i1,i5),-s(i2,i4),-s(i1,i3),-s(i2,i4))
      boxint(8)=Lsm1(-s(i2,i4),-s(i1,i3),-s(i2,i5),-s(i1,i3))
      boxint(9)=Lsm1(-s(i2,i5),-s(i1,i4),-s(i2,i3),-s(i1,i4))
      boxint(10)=Lsm1(-s(i1,i4),-s(i2,i3),-s(i1,i5),-s(i2,i3))

!------ box coefficients 
      boxcoeff(1)=  (za(i1,i2)**2*za(i4,i5)**2)/
     -  (za(i1,i5)*za(i2,i5)**3*za(i3,i4))
      boxcoeff(2)=ALO
      boxcoeff(3)= ALO
      boxcoeff(4)=-(za(i1,i3)**2*za(i2,i4)**2)/
     -  (za(i1,i5)*za(i2,i3)**2*za(i2,i5)*
     -     za(i3,i4))
      boxcoeff(5)=-ALO
      boxcoeff(6)=ALO
      boxcoeff(7)=-(za(i1,i3)**2*za(i4,i5)**2)/
     -  (za(i1,i5)*za(i2,i5)*za(i3,i4)*
     -     za(i3,i5)**2)
      boxcoeff(8)=-ALO
      boxcoeff(9) = ALO
      boxcoeff(10) = ALO

!      write(6,*) boxcoeff(1)/2._dp*(s(i1,i5)*s(i1,i2))*im
!      write(6,*) boxcoeff(2)/2._dp*(s(i2,i5)*s(i1,i2))*im
!      write(6,*) boxcoeff(3)/2._dp*(s(i3,i4)*s(i2,i3))*im
!      write(6,*) boxcoeff(4)/2._dp*(s(i3,i4)*s(i2,i4))*im
!      write(6,*) boxcoeff(5)/2._dp*(s(i3,i4)*s(i1,i3))*im
!      write(6,*) boxcoeff(6)/2._dp*(s(i3,i4)*s(i1,i4))*im
!      write(6,*) boxcoeff(7)/2._dp*(s(i1,i5)*s(i1,i3))*im
!      write(6,*) boxcoeff(8)/2._dp*(s(i2,i5)*s(i2,i4))*im
!      write(6,*) boxcoeff(9)/2._dp*(s(i2,i5)*s(i2,i3))*im
!      write(6,*) boxcoeff(10)/2._dp*(s(i1,i5)*s(i1,i4))*im

      Boxes=czip
      do i=1,10
         Boxes=Boxes+boxint(i)*boxcoeff(i)
      enddo
      
      bubs = (lnrat(-s(i1,i5),-s(i3,i4))*za(i1,i5)*
     -     za(i2,i4)**2)/(2.*za(i2,i5)**3*za(i3,i4)) + 
     -  (lnrat(-s(i1,i5),-s(i3,i4))*za(i1,i3)*za(i1,i4))/
     -   (za(i1,i5)*za(i2,i3)*za(i3,i5)) - 
     -  (lnrat(-s(i1,i3),-s(i3,i4))*za(i1,i4))/
     -   (za(i2,i5)*za(i3,i5)) - 
     -  (L1(-s(i1,i5),-s(i3,i4))*za(i1,i2)**2*za(i1,i5)*
     -     za(i2,i4)**2*zb(i2,i1)**2)/
     -   (2.*s(i3,i4)**2*za(i2,i5)**3*za(i3,i4)) + 
     -  (L0(-s(i2,i4),-s(i1,i5))*za(i1,i3)**2*za(i2,i4)*
     -     zb(i3,i2))/
     -   (s(i1,i5)*za(i1,i5)*za(i2,i3)*za(i3,i5)) - 
     -  (L0(-s(i1,i5),-s(i3,i4))*za(i1,i2)*za(i2,i4)*
     -     zb(i5,i2))/(s(i3,i4)*za(i2,i3)*za(i2,i5)) - 
     -  (L0(-s(i2,i4),-s(i1,i3))*za(i1,i5)*za(i2,i4)*
     -     zb(i5,i2))/(s(i1,i3)*za(i2,i5)*za(i3,i5)) + 
     -  (2*L0(-s(i1,i2),-s(i3,i4))*za(i1,i2)*za(i2,i4)*
     -     za(i4,i5)*zb(i5,i2))/
     -   (s(i3,i4)*za(i2,i5)**2*za(i3,i4)) + 
     -  (2*L0(-s(i1,i5),-s(i3,i4))*za(i1,i2)*za(i2,i4)*
     -     za(i4,i5)*zb(i5,i2))/
     -   (s(i3,i4)*za(i2,i5)**2*za(i3,i4)) + 
     -  (L1(-s(i1,i2),-s(i3,i4))*za(i1,i2)**2*
     -     za(i4,i5)**2*zb(i5,i1)*zb(i5,i2))/
     -     (s(i3,i4)**2*za(i2,i5)**2*za(i3,i4))


      Rat=(-(za(i1,i4)**2/(za(i1,i5)*za(i2,i5)*za(i3,i4))) + 
     -  (za(i1,i2)*za(i2,i4)*za(i4,i5)*zb(i5,i2))/
     -   (za(i2,i5)**2*za(i3,i4)**2*zb(i4,i3)) + 
     -  (za(i1,i4)*zb(i5,i3))/
     -     (2.*za(i2,i5)*za(i3,i4)*zb(i4,i3)))
      
!===  plus correction for deviation from del1 slc
      Rat=Rat-(-(za(i1,i2)*za(i2,i4)*
     -     (za(i1,i5)*za(i2,i4)*zb(i2,i1) + 
     -       za(i2,i5)*za(i4,i5)*zb(i5,i2)))/
     -  (2.*za(i2,i5)**3*za(i3,i4)**2*zb(i4,i3)))
      
      amp_qqbQQbga_mhvalt_del2_slc=Vpole+boxes-bubs-Rat
     
!      write(6,*) 'total ', amp_qqbQQbga_mhvalt_del2_slc*im
!      write(6,*) 'pole pieces',Vpole*im
!      write(6,*) 'boxes ',Boxes*im
!      write(6,*) 'bubs ',Bubs*im
!      write(6,*) 'completed rat',rat*im
!      write(6,*) 'total/LO',amp_qqbQQbga_mhvalt_del2_slc*im/ALO


      return
      end



      
