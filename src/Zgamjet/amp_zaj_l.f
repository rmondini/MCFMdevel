!---- CW Jan 17
!---- reimplementaiton of Amplitudes for qqb + ell +ell + ga + g where ga is
!---- radiated from lepton line
!---- two basic amplitudes corresponding to two helicity choices for photon



      function amp_zaj_ql_lc_MHV(i1,i2,i3,i4,i5,i6,za,zb)
      implicit none
      include 'types.f'
      complex(dp) :: amp_zaj_ql_lc_MHV
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f' 
      include 'scale.f'
      include 'epinv.f' 
      integer i1,i2,i3,i4,i5,i6
      complex(dp) Atree, Vpole,Box(1),Boxint(1),rat,bubs,boxes
      complex(dp) lnrat,L0,L1,Lsm1
      real(dp) t
      
      t(i1,i2,i3)=s(i1,i2)+s(i2,i3)+s(i1,i3)
      
!---  lc ampltidue for q^-(i1) + q^+(i2) + ell^-(i3)+ell^+(i4)+ga^+(i5)+g^+(i6)

      
      Atree =za(i2,i3)**2/(za(i1,i6)*za(i2,i6)*za(i3,i5)*za(i4,i5))

      Vpole =-2*epinv**2 + epinv*(-1.5 - lnrat(musq,-s(i1,i6)) - 
     -     lnrat(musq,-s(i2,i6))) + 
     -  (-6 - lnrat(musq,-s(i1,i6))**2 - 3*lnrat(musq,-s(i2,i6)) - 
     -     lnrat(musq,-s(i2,i6))**2)/2.
      
      Vpole=Vpole*Atree
!===== boxes
      Boxint(1)=Lsm1(-s(i1,i6),-t(i1,i2,i6),-s(i2,i6),-t(i1,i2,i6))
      Box(1)=-Atree

!      write(6,*) Box(1)*(s(i2,i6)*s(i1,i6)/2._dp)

      
      boxes=Boxint(1)*Box(1)         
      
!---  completed bubbles
      Bubs=  -((L0(-s(i2,i6),-t(i3,i4,i5))*za(i1,i2)*za(i2,i3)*
     -       (za(i3,i4)*zb(i4,i1) + za(i3,i5)*zb(i5,i1)))/
     -     (t(i3,i4,i5)*za(i1,i6)*za(i2,i6)*za(i3,i5)*za(i4,i5))) + 
     -  (L1(-s(i2,i6),-t(i3,i4,i5))*za(i1,i2)**2*
     -     (za(i3,i4)*zb(i4,i1) + za(i3,i5)*zb(i5,i1))*
     -     (-(za(i2,i3)*zb(i2,i1)) + za(i3,i6)*zb(i6,i1)))/
     -     (2.*t(i3,i4,i5)**2*za(i1,i6)*za(i2,i6)*za(i3,i5)*za(i4,i5))
      
!----- additional boxes      
      Rat=czip

      
      amp_zaj_ql_lc_MHV = Vpole+boxes+bubs+rat

!---- check
!      write(6,*) 'Atree* im  = ',Atree*im
!      write(6,*) 'Boxes  = ',Boxes
!      write(6,*) 'Completed Bubs   = ',bubs
!      write(6,*) 'Rat   = ',Rat
!      write(6,*) 'Total  = ',amp_zaj_ql_lc_MHV
      
 
!      stop
      return
      end
      

      

      function amp_zaj_ql_slc_MHV(i1,i2,i3,i4,i5,i6,za,zb)
      implicit none
      include 'types.f'
      complex(dp) :: amp_zaj_ql_slc_MHV
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f' 
      include 'scale.f'
      include 'epinv.f' 
         integer i1,i2,i3,i4,i5,i6
      complex(dp) Atree, Vpole,Box(2),Boxint(2),rat,bubs,boxes
      complex(dp) lnrat,L0,L1,Lsm1
      real(dp) t
      integer i
      
      t(i1,i2,i3)=s(i1,i2)+s(i2,i3)+s(i1,i3)
      
!---  slc ampltidue for q^-(i1) + q^+(i2) + ell^-(i3)+ell^+(i4)+ga^+(i5)+g^+(i6)

      
      Atree =za(i2,i3)**2/(za(i1,i6)*za(i2,i6)*za(i3,i5)*za(i4,i5))

      Vpole=  (  epinv**2 + epinv*(1.5 + lnrat(musq,-s(i1,i2))) + 
     -     (6 + lnrat(musq,-s(i1,i2))**2 + 3*lnrat(musq,-t(i3,i4,i5)))/2)
      
      Vpole=Vpole*Atree
      
      Boxint(1)=Lsm1(-s(i1,i2),-t(i1,i2,i6),-s(i2,i6),-t(i1,i2,i6))
      Boxint(2)=Lsm1(-s(i1,i6),-t(i1,i2,i6),-s(i1,i2),-t(i1,i2,i6))

      Box(1)=(za(i1,i2)**2*za(i3,i6)**2)
     &     /(za(i1,i6)**3*za(i2,i6)*za(i3,i5)*za(i4,i5))
      Box(2)=Atree

 !     write(6,*) Box(1)*(s(i1,i2)*s(i2,i6))/2._dp
 !     write(6,*) Box(2)*(s(i1,i2)*s(i1,i6))/2._dp
      
      Boxes=czip
      do i=1,2
         Boxes=Boxes+Box(i)*Boxint(i)
      enddo

      bubs= (L0(-t(i3,i4,i5),-s(i1,i2))*za(i1,i2)*za(i1,i3)*za(i3,i6)*
     -     zb(i6,i1))/(s(i1,i2)*za(i1,i6)**2*za(i3,i5)*za(i4,i5)) - 
     -  (L1(-t(i3,i4,i5),-s(i2,i6))*za(i1,i3)**2*za(i2,i6)*zb(i6,i1)**2)/
     -   (2.*s(i2,i6)**2*za(i1,i6)*za(i3,i5)*za(i4,i5)) - 
     -  (L0(-t(i3,i4,i5),-s(i2,i6))*
     -     (2*s(i1,i6)*za(i1,i3)*za(i2,i3) - 
     -       za(i1,i3)**2*za(i2,i6)*zb(i6,i1)))/
     -   (s(i2,i6)*za(i1,i6)**2*za(i3,i5)*za(i4,i5)) + 
     -  (L1(-t(i3,i4,i5),-s(i1,i2))*za(i1,i2)*za(i3,i6)*zb(i6,i1)*
     -     (-(za(i3,i4)*zb(i6,i4)) - za(i3,i5)*zb(i6,i5)))/
     -   (s(i1,i2)**2*za(i1,i6)*za(i3,i5)*za(i4,i5))

      Rat=   (za(i1,i3)*(za(i3,i4)*zb(i4,i2) + za(i3,i5)*zb(i5,i2))*zb(i6,i1)**2)/
     -   (2.*t(i3,i4,i5)*za(i1,i6)*za(i3,i5)*za(i4,i5)*zb(i2,i1)*zb(i6,i2)) + 
     -  ((za(i3,i4)*zb(i4,i1) + za(i3,i5)*zb(i5,i1))*
     -     (za(i1,i3)*zb(i6,i1) + za(i2,i3)*zb(i6,i2)))/
     -   (2.*t(i3,i4,i5)*za(i1,i6)*za(i3,i5)*za(i4,i5)*zb(i2,i1)) + 
     -  (za(i2,i3)*(za(i3,i4)*zb(i4,i1) + za(i3,i5)*zb(i5,i1))*
     -     (za(i1,i2)*zb(i2,i1) + za(i2,i6)*zb(i6,i2)))/
     -   (2.*t(i3,i4,i5)*za(i1,i6)*za(i2,i6)*za(i3,i5)*za(i4,i5)*zb(i2,i1))

      
      amp_zaj_ql_slc_MHV = Vpole+boxes+bubs+rat
      
      
!---- check
!      write(6,*) 'Atree* im  = ',Atree*im
!      write(6,*) 'Vpole* im  = ',Vpole
!      write(6,*) 'Boxes* im  = ',Boxes
!      write(6,*) 'Completed Bubs* im  = ',bubs
!      write(6,*) 'Rat*im  = ',Rat
!     write(6,*) 'Total  = ',amp_zaj_ql_slc_MHV
  
      
!      stop
      return
      end


      function amp_zaj_ql_lc_NMHV(i1,i2,i3,i4,i5,i6,za,zb)
      implicit none
      include 'types.f'
      complex(dp) :: amp_zaj_ql_lc_NMHV
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f' 
      include 'scale.f'
      include 'epinv.f' 
      integer i1,i2,i3,i4,i5,i6
      complex(dp) Atree, Vpole,Box(1),Boxint(1),rat,bubs,boxes
      complex(dp) lnrat,L0,L1,Lsm1
      real(dp) t
      
      t(i1,i2,i3)=s(i1,i2)+s(i2,i3)+s(i1,i3)
      
!---  lc ampltidue for q^-(i1) + q^+(i2) + ell^-(i3)+ell^+(i4)+ga^-(i5)+g^+(i6)

      
      Atree = (za(i1,i2)*zb(i4,i1) + za(i2,i6)*zb(i6,i4))**2/
     -  (t(i1,i2,i6)*za(i1,i6)*za(i2,i6)*zb(i5,i3)*zb(i5,i4))

      Vpole = -2*epinv**2 + epinv*
     -   (-1.5 - lnrat(musq,-s(i1,i6)) - lnrat(musq,-s(i2,i6))) + 
     -  (-6 - lnrat(musq,-s(i1,i6))**2 - 
     -     lnrat(musq,-s(i2,i6))**2 - 3*lnrat(musq,-t(i3,i4,i5)))/
     -     2.
      
      Vpole=Vpole*Atree
!===== boxes
      Boxint(1)=Lsm1(-s(i1,i6),-t(i1,i2,i6),-s(i2,i6),-t(i1,i2,i6))
      Box(1)=-Atree

!      write(6,*) Box(1)*(s(i2,i6)*s(i1,i6)/2._dp)

      
      boxes=Boxint(1)*Box(1)         
      
!---  completed bubbles
      Bubs=  (L1(-s(i2,i6),-t(i3,i4,i5))*za(i1,i2)*zb(i4,i1)*zb(i6,i1)*
     -     (-(za(i1,i2)*zb(i4,i2)) + za(i1,i6)*zb(i6,i4)))/
     -   (2.*t(i3,i4,i5)**2*za(i1,i6)*zb(i5,i3)*zb(i5,i4)) + 
     -     (3*L0(-s(i2,i6),-t(i3,i4,i5))*zb(i6,i1)*(-(za(i1,i2)*zb(i4,i2))
     & + za(i1,i6)*zb(i6,i4))*
     -     (za(i1,i2)*zb(i4,i1) + za(i2,i6)*zb(i6,i4)))/
     -   (2.*t(i3,i4,i5)**2*za(i1,i6)*zb(i5,i3)*zb(i5,i4))
      
!----- additional boxes      
      Rat=
     & -(za(i1,i2)*zb(i4,i1)*(za(i1,i2)*zb(i4,i1) + za(i2,i6)*zb(i6,i4)))/
     -  (2.*t(i3,i4,i5)*za(i1,i6)*za(i2,i6)*zb(i5,i3)*zb(i5,i4))

      
      amp_zaj_ql_lc_NMHV = Vpole+boxes+bubs+rat

!---- check
!      write(6,*) 'Atree* im  = ',Atree*im
!      write(6,*) 'Vpole  = ',Vpole
!      write(6,*) 'Boxes  = ',Boxes
!      write(6,*) 'Completed Bubs   = ',bubs
!      write(6,*) 'Rat   = ',Rat
!      write(6,*) 'Total  = ',amp_zaj_ql_lc_NMHV
      
 
!      stop
      return
      end
      


      
      function amp_zaj_ql_slc_NMHV(i1,i2,i3,i4,i5,i6,za,zb)
      implicit none
      include 'types.f'
      complex(dp) :: amp_zaj_ql_slc_NMHV
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f' 
      include 'scale.f'
      include 'epinv.f' 
         integer i1,i2,i3,i4,i5,i6
      complex(dp) Atree, Vpole,Box(2),Boxint(2),rat,bubs,boxes
      complex(dp) lnrat,L0,L1,Lsm1
      real(dp) t
      integer i
      
      t(i1,i2,i3)=s(i1,i2)+s(i2,i3)+s(i1,i3)
      
!---  slc ampltidue for q^-(i1) + q^+(i2) + ell^-(i3)+ell^+(i4)+ga^-(i5)+g^+(i6)

      
      Atree =(za(i1,i2)*zb(i4,i1) + za(i2,i6)*zb(i6,i4))**2/
     -  (t(i1,i2,i6)*za(i1,i6)*za(i2,i6)*zb(i5,i3)*zb(i5,i4))

      
      Vpole=  (  epinv**2 + epinv*(1.5 + lnrat(musq,-s(i1,i2))) + 
     -     (6 + lnrat(musq,-s(i1,i2))**2 + 3*lnrat(musq,-t(i3,i4,i5)))/2)
      
      Vpole=Vpole*Atree
      
      Boxint(1)=Lsm1(-s(i1,i2),-t(i1,i2,i6),-s(i2,i6),-t(i1,i2,i6))
      Boxint(2)=Lsm1(-s(i1,i6),-t(i1,i2,i6),-s(i1,i2),-t(i1,i2,i6))

      Box(1)= -((za(i1,i2)**2*(za(i1,i6)*zb(i4,i1) + 
     -        za(i2,i6)*zb(i4,i2))*
     -      (za(i3,i6)*zb(i4,i3) - za(i5,i6)*zb(i5,i4)))/
     -    (t(i3,i4,i5)*za(i1,i6)**3*za(i2,i6)*zb(i5,i3)*zb(i5,i4)))
      Box(2)=Atree

!      write(6,*) Box(1)*(s(i1,i2)*s(i2,i6))/2._dp
!      write(6,*) Box(2)*(s(i1,i2)*s(i1,i6))/2._dp
      
      Boxes=czip
      do i=1,2
         Boxes=Boxes+Box(i)*Boxint(i)
      enddo

      bubs=  -((L1(-s(i1,i2),-t(i3,i4,i5))*za(i1,i2)**2*
     -       (za(i1,i6)*zb(i4,i1) + za(i2,i6)*zb(i4,i2))**2*zb(i6,i1)*
     -       zb(i6,i2))/
     -     (t(i3,i4,i5)**3*za(i1,i6)**2*zb(i5,i3)*zb(i5,i4))) - 
     -  (2*L0(-s(i1,i2),-t(i3,i4,i5))*za(i1,i2)*
     -     (za(i1,i6)*zb(i4,i1) + za(i2,i6)*zb(i4,i2))*zb(i6,i1)*
     -     (-(za(i1,i2)*zb(i4,i2)) + za(i1,i6)*zb(i6,i4)))/
     -   (t(i3,i4,i5)**2*za(i1,i6)**2*zb(i5,i3)*zb(i5,i4)) + 
     -  (L1(-t(i3,i4,i5),-s(i2,i6))*za(i2,i6)*
     -     (-(za(i1,i3)*zb(i4,i3)) + za(i1,i5)*zb(i5,i4))*zb(i6,i1)**2*
     -     (-(za(i1,i2)*zb(i4,i2)) + za(i1,i6)*zb(i6,i4)))/
     -   (2.*s(i2,i6)**2*t(i3,i4,i5)*za(i1,i6)*zb(i5,i3)*zb(i5,i4)) + 
     -  (L0(-t(i3,i4,i5),-s(i2,i6))*
     -     (2*s(i1,i6)*(-(za(i2,i3)*zb(i4,i3)) + 
     -          za(i2,i5)*zb(i5,i4)) - 
     -       za(i2,i6)*(-(za(i1,i3)*zb(i4,i3)) + za(i1,i5)*zb(i5,i4))*
     -        zb(i6,i1))*(-(za(i1,i2)*zb(i4,i2)) + za(i1,i6)*zb(i6,i4))
     -     )/(s(i2,i6)*t(i3,i4,i5)*za(i1,i6)**2*zb(i5,i3)*zb(i5,i4))

      Rat=(za(i1,i2)**2*zb(i4,i2)*(za(i1,i6)*zb(i4,i1) + za(i2,i6)*zb(i4,i2))*zb(i6,i1))/
     -   (2.*t(i3,i4,i5)**2*za(i1,i6)**2*zb(i5,i3)*zb(i5,i4)) + 
     -  (za(i1,i2)*(za(i1,i6)*zb(i4,i1) + za(i2,i6)*zb(i4,i2))*
     -     (-(za(i1,i3)*zb(i4,i3)) + za(i1,i5)*zb(i5,i4))*zb(i6,i1))/
     -   (2.*t(i3,i4,i5)**2*za(i1,i6)**2*zb(i5,i3)*zb(i5,i4)) + 
     -  (za(i1,i2)*zb(i2,i1)*(-(za(i1,i3)*zb(i4,i3)) + za(i1,i5)*zb(i5,i4))*zb(i6,i1)*
     -     zb(i6,i4))/(2.*t(i3,i4,i5)**2*za(i1,i6)*zb(i5,i3)*zb(i5,i4)*zb(i6,i2)) - 
     -  (zb(i6,i1)*zb(i6,i4)**2)/(2.*t(i3,i4,i5)*zb(i5,i3)*zb(i5,i4)*zb(i6,i2)) + 
     -  (zb(i6,i4)*(za(i1,i2)*zb(i4,i1) + za(i2,i6)*zb(i6,i4)))/
     -   (2.*t(i3,i4,i5)*za(i1,i6)*zb(i5,i3)*zb(i5,i4)) + 
     -  (za(i1,i2)*zb(i4,i1) + za(i2,i6)*zb(i6,i4))**2/
     -     (2.*t(i3,i4,i5)*za(i1,i6)*za(i2,i6)*zb(i5,i3)*zb(i5,i4))
      
      amp_zaj_ql_slc_NMHV = Vpole+boxes+bubs+rat
      
      
!---- check
!      write(6,*) 'Atree* im  = ',Atree*im
!      write(6,*) 'Vpole* im  = ',Vpole
!      write(6,*) 'Boxes* im  = ',Boxes
!      write(6,*) 'Completed Bubs* im  = ',bubs
!      write(6,*) 'Rat*im  = ',Rat
!      write(6,*) 'Total  = ',amp_zaj_ql_slc_NMHV
  
      
!      stop
      return
      end
