!---- CW Sept 16
!---- this is the amplitude for the following processes

!     in = u[i1],D[i2],N[i3],e[i4],ga[i5],g[i6];
!=====this is the prp34 piece
      
!---- with helicity configuration
!-------+  -  + - ++

!======pieces are broken into Qu = Q(i1) and Qd = Q(i2)
      
      function Wgajprp34_Qu_mhv_lc(i1,i2,i3,i4,i5,i6,za,zb)
      implicit none
      include 'types.f'
      complex(dp) Wgajprp34_Qu_mhv_lc
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f' 
      include 'scale.f'
      include 'epinv.f'
      include 'epinv2.f'
      integer :: i1,i2,i3,i4,i5,i6
      complex(dp) :: Atree,Boxes,Bubs,Rat,Vpole
      complex(dp) :: Boxcoeff(3),boxint(3)
      integer i
      real(dp)::t 
      complex(dp) :: zab2
      complex(dp) :: Lsm1_2me,Lsm1
      complex(dp) :: L0,L1,L2,lnrat
      zab2(i1,i2,i3,i4)=za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)
      t(i1,i2,i3)=s(i1,i2)+s(i2,i3)+s(i1,i3)

      
      Wgajprp34_Qu_mhv_lc=czip
      rat=czip
      bubs=czip
      boxes=czip
!===== poles 
      Atree=  -(za(i2,i4)**2*zab2(i1,i3,i4,i5))/
     -  (za(i1,i5)*za(i1,i6)*za(i2,i6)*za(i3,i4)*zab2(i5,i3,i4,i5))

      write(6,*) Atree
      pause
      
      Vpole=    2*epinv**2 + epinv*
     &  (1.5_dp + lnrat(musq,-s(i1,i6)) + 
     &     lnrat(musq,-s(i2,i6))) + 
     &  (6._dp + lnrat(musq,-s(i1,i6))**2 + 
     &     lnrat(musq,-s(i2,i6))**2 + 
     &     3*lnrat(musq,-t(i1,i2,i6)))/2.
      

      Vpole=Atree*Vpole
!=====boxes

      boxint(1)=Lsm1_2me(t(i1,i5,i6),t(i1,i2,i6),s(i1,i6),s(i3,i4))
      boxint(2)=Lsm1(-s(i1,i5),-t(i1,i5,i6),-s(i1,i6),-t(i1,i5,i6))
      boxint(3)=Lsm1(-s(i1,i6),-t(i1,i2,i6),-s(i2,i6),-t(i1,i2,i6))

      boxcoeff(1)=(za(i1,i2)*za(i2,i4)**2)
     &     /(za(i1,i5)*za(i1,i6)*za(i2,i5)*za(i2,i6)*za(i3,i4))
      boxcoeff(2)=-(za(i2,i4)**2
     &     /(za(i1,i5)*za(i2,i6)*za(i3,i4)*za(i5,i6)))
      boxcoeff(3)=-(za(i2,i4)**2*zab2(i1,i3,i4,i5))/
     -  (za(i1,i5)*za(i1,i6)*za(i2,i6)*za(i3,i4)*zab2(i5,i3,i4,i5))
      
      boxes=czip
      do i=1,3
         boxes=boxes+boxcoeff(i)*boxint(i)
      enddo

!      write(6,*) boxcoeff(1)*zab2(i5,i3,i4,i2)*zab2(i2,i3,i4,i5)/2._dp
!      write(6,*) boxcoeff(2)*s(i1,i6)*s(i1,i5)/2._dp
!     write(6,*) boxcoeff(3)*s(i1,i6)*s(i2,i6)/2._dp

     
!=====Completed bubbles
      bubs=  -((L0(-s(i3,i4),-t(i1,i2,i6))*za(i1,i2)*
     -       za(i2,i4)*zb(i5,i3))/
     -     (t(i1,i2,i6)*za(i1,i5)*za(i1,i6)*za(i2,i6))
     -     ) - (L1(-s(i3,i4),-t(i1,i2,i6))*za(i1,i2)*
     -     za(i2,i5)*za(i3,i4)*zb(i5,i3)**2)/
     -   (2.*t(i1,i2,i6)**2*za(i1,i5)*za(i1,i6)*
     -     za(i2,i6)) - 
     -  (2*L0(-s(i2,i6),-t(i1,i2,i6))*za(i1,i3)*
     -     za(i1,i4)*za(i2,i4)*zb(i5,i3)*zb(i6,i1))/
     -   (t(i1,i2,i6)*za(i1,i5)*za(i1,i6)*za(i3,i4)*
     -     (za(i3,i5)*zb(i5,i3) + za(i4,i5)*zb(i5,i4))
     -     ) - (2*L0(-s(i2,i6),-t(i1,i2,i6))*
     -     za(i1,i4)**2*za(i2,i4)*zb(i5,i4)*zb(i6,i1))
     -    /(t(i1,i2,i6)*za(i1,i5)*za(i1,i6)*za(i3,i4)*
     -     (za(i3,i5)*zb(i5,i3) + za(i4,i5)*zb(i5,i4))
     -     ) + (L1(-s(i2,i6),-t(i1,i2,i6))*za(i1,i3)*
     -     za(i1,i4)**2*za(i2,i6)*zb(i5,i3)*
     -     zb(i6,i1)**2)/
     -   (2.*t(i1,i2,i6)**2*za(i1,i5)*za(i1,i6)*
     -     za(i3,i4)*(za(i3,i5)*zb(i5,i3) + 
     -       za(i4,i5)*zb(i5,i4))) + 
     -  (L1(-s(i2,i6),-t(i1,i2,i6))*za(i1,i4)**3*
     -     za(i2,i6)*zb(i5,i4)*zb(i6,i1)**2)/
     -   (2.*t(i1,i2,i6)**2*za(i1,i5)*za(i1,i6)*
     -     za(i3,i4)*(za(i3,i5)*zb(i5,i3) + 
     -       za(i4,i5)*zb(i5,i4)))
      
!      write(6,*) bubs
!      write(6,*) L0(-s(i3,i4),-t(i1,i2,i6))
!      write(6,*) L0(-s(i3,i4),-t(i1,i2,i6))*bubs
      
!      pause
      rat=czip
      
!-----total 
      Wgajprp34_Qu_mhv_lc=Vpole+boxes+bubs+rat

      write(6,*) 'Warning rational pieces for Wgaj_mhv_lc not there'
      
!======writeout
      write(6,*) 'Atree *im',Atree*im
      write(6,*) 'Vpole*im',Vpole/Atree*im
      write(6,*) 'Boxes*im',Boxes/Atree*im
      write(6,*) 'Bubs*im',Bubs/Atree*im
      write(6,*) 'CC *im',(Vpole+Bubs+Boxes)/Atree*im
      write(6,*) 'total*im',Wgajprp34_Qu_mhv_lc
      return
      end
      
