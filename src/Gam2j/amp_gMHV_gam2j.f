!-----collection of amplitudes for q(i1)^- q(i2)^+ g(i3)^- g(i4)^+ gamma(i5)
!=====C. Williams Spring 16

!========LEADING COLOR AMPLITUDE
      function amp_qqbggga_lc_gMHV(i1,i2,i3,i4,i5,za,zb)
      implicit none
      include 'types.f'
      complex(dp) :: amp_qqbggga_lc_gMHV
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

!------ boxcoeffs and ints
      ALO=-(za(i1,i3)**3/(za(i1,i4)*za(i1,i5)*za(i2,i5)*za(i3,i4)))


            
!-----pole expansion
      Vpole = (3*epinv*epinv2 + epinv*
     -   (1.5 + lnrat(musq,-s(i1,i4)) + 
     -     lnrat(musq,-s(i2,i3)) + lnrat(musq,-s(i3,i4))) + 
     -  (6 + lnrat(musq,-s(i1,i4))**2 + 
     -     lnrat(musq,-s(i2,i3))**2 + 
     -     3*lnrat(musq,-s(i3,i4)) + lnrat(musq,-s(i3,i4))**2)
     -     /2._dp)

      Vpole=Vpole*ALO

      boxcoeff(1)=ALO
      boxcoeff(2)= -((za(i1,i4)**2*za(i2,i3)**3)/
     -     (za(i1,i5)*za(i2,i4)**3*za(i2,i5)*za(i3,i4)))
      boxcoeff(3)= -((za(i1,i4)**2*za(i3,i5)**3)/
     -    (za(i1,i5)*za(i2,i5)*za(i3,i4)*za(i4,i5)**3))
      boxcoeff(4)=ALO
      
      Boxint(1)=Lsm1(-s(i2,i5),-s(i1,i4),-s(i2,i3),-s(i1,i4))
      Boxint(2)=Lsm1(-s(i3,i4),-s(i1,i5),-s(i2,i3),-s(i1,i5))
      Boxint(3)=Lsm1(-s(i1,i5),-s(i2,i3),-s(i1,i4),-s(i2,i3))
      Boxint(4)=Lsm1(-s(i3,i4),-s(i2,i5),-s(i1,i4),-s(i2,i5))
      
!------Boxes
      Boxes=czip
      do i=1,4
         Boxes=Boxes+boxcoeff(i)*Boxint(i)
      enddo

!-----(Completed) Bubs
      Bubs= -( lnrat(-s(i2,i3),-s(i3,i4))*
     -   ((za(i1,i3)*za(i3,i5))/(za(i2,i5)*za(i4,i5)**2) + 
     -     (3*za(i1,i3)**2)/(2.*za(i1,i4)*za(i2,i5)*za(i4,i5)))
     -    - (lnrat(-s(i1,i5),-s(i3,i4))*za(i1,i3)*
     -     (za(i1,i3)*za(i2,i4)*za(i4,i5) + 
     -       2*za(i1,i4)*(za(i2,i4)*za(i3,i5) + 
     -          za(i2,i3)*za(i4,i5))))/
     -   (2.*za(i1,i5)*za(i2,i4)**2*za(i4,i5)**2) - 
     -  (L0(-s(i2,i3),-s(i1,i5))*za(i1,i4)**2*za(i2,i3)*
     -     za(i3,i5)*zb(i4,i2))/
     -   (s(i1,i5)*za(i1,i5)*za(i2,i4)*za(i4,i5)**2) - 
     -  (L0(-s(i2,i3),-s(i1,i5))*za(i1,i4)**2*za(i2,i3)**2*
     -     zb(i4,i2))/
     -   (s(i1,i5)*za(i1,i5)*za(i2,i4)**2*za(i4,i5)) - 
     -  (L1(-s(i2,i3),-s(i1,i5))*za(i1,i4)**2*za(i2,i3)**2*
     -     zb(i4,i2)**2)/
     -   (2.*s(i1,i5)**2*za(i1,i5)*za(i2,i4)*za(i4,i5)) - 
     -  (L0(-s(i1,i5),-s(i3,i4))*za(i1,i4)*za(i2,i3)**3*
     -     zb(i5,i2))/
     -   (s(i3,i4)*za(i2,i4)**2*za(i2,i5)*za(i3,i4)) - 
     -  (L1(-s(i1,i5),-s(i3,i4))*za(i1,i5)*za(i2,i3)**3*
     -     zb(i5,i2)**2)/
     -   (2.*s(i3,i4)**2*za(i2,i4)*za(i2,i5)*za(i3,i4)) + 
     -  (L1(-s(i1,i4),-s(i2,i3))*za(i1,i4)*za(i3,i5)**2*
     -     zb(i5,i4)**2)/(2.*s(i2,i3)**2*za(i2,i5)*za(i4,i5))
     -   + (L0(-s(i1,i4),-s(i2,i3))*
     -     ((za(i1,i4)*za(i3,i5)**2*zb(i5,i4))/
     -        (za(i2,i5)*za(i4,i5)**2) + 
     -       (za(i1,i3)*za(i3,i5)*zb(i5,i4))/
     -     (za(i2,i5)*za(i4,i5))))/s(i2,i3))

      Rat=  za(i1,i3)**3/
     -   (4.*za(i1,i4)*za(i1,i5)*za(i2,i5)*za(i3,i4)) + 
     -  (za(i1,i4)**2*za(i2,i3)*zb(i4,i2))/
     -   (4.*za(i1,i5)*za(i2,i4)**2*za(i4,i5)*zb(i3,i2)) - 
     -  (za(i1,i2)*za(i1,i3)*za(i2,i3)*zb(i4,i2))/
     -   (4.*za(i1,i5)*za(i2,i4)*za(i2,i5)*za(i3,i4)*zb(i4,i3))
     -    + (za(i1,i2)**3*za(i3,i4)*zb(i4,i2))/
     -   (4.*za(i1,i4)*za(i1,i5)*za(i2,i4)**2*za(i2,i5)*
     -     zb(i4,i3)) + (za(i1,i2)*za(i1,i3)*za(i3,i4)*
     -     zb(i4,i3))/
     -   (4.*za(i1,i4)*za(i1,i5)*za(i2,i4)*za(i2,i5)*zb(i3,i1))
     -    - (za(i1,i3)**2*zb(i4,i2))/
     -   (4.*za(i1,i5)**2*za(i4,i5)*zb(i5,i1)) - 
     -  (za(i1,i3)*za(i1,i4)*za(i2,i3)*zb(i4,i2))/
     -   (4.*za(i1,i5)**2*za(i2,i4)*za(i4,i5)*zb(i5,i1)) - 
     -  (za(i1,i2)**2*za(i3,i4)*zb(i4,i2)**2)/
     -   (4.*za(i1,i5)**2*za(i2,i4)*za(i4,i5)*zb(i4,i3)*
     -     zb(i5,i1)) + (za(i1,i5)**2*za(i2,i3)*zb(i5,i2))/
     -   (4.*za(i1,i4)*za(i2,i5)**2*za(i4,i5)*zb(i3,i2)) - 
     -  (za(i1,i3)*za(i2,i3)**2*zb(i5,i2))/
     -   (2.*za(i2,i4)*za(i2,i5)*za(i3,i4)**2*zb(i4,i3)) - 
     -  (za(i1,i2)*za(i3,i5)**2*zb(i4,i1)*zb(i5,i2))/
     -   (4.*za(i2,i5)**2*za(i3,i4)*za(i4,i5)*zb(i2,i1)*
     -     zb(i4,i3)) - (za(i1,i2)**2*zb(i4,i2)*zb(i5,i2))/
     -   (4.*za(i1,i4)*za(i2,i4)*za(i2,i5)*zb(i3,i2)*zb(i4,i3))
     -    + (za(i1,i3)*zb(i4,i3)*zb(i5,i2))/
     -   (4.*za(i1,i4)*za(i2,i5)*zb(i3,i1)*zb(i3,i2)) + 
     -  (za(i2,i3)**3*zb(i5,i2))/
     -   (4.*za(i2,i4)*za(i2,i5)**2*za(i3,i4)*zb(i5,i1)) + 
     -  (za(i1,i2)*zb(i4,i1)**3*zb(i5,i2))/
     -   (4.*za(i2,i5)**2*zb(i2,i1)*zb(i3,i1)*zb(i4,i3)*
     -     zb(i5,i1)) + (za(i1,i2)**2*za(i1,i3)*zb(i3,i2)*
     -     zb(i5,i1))/
     -   (4.*za(i1,i4)*za(i1,i5)*za(i2,i4)*za(i2,i5)*zb(i3,i1)*
     -     zb(i5,i3)) + (za(i1,i2)*zb(i4,i1)*zb(i5,i2))/
     -   (4.*za(i2,i5)*za(i4,i5)*zb(i3,i1)*zb(i5,i3)) + 
     -  (za(i1,i2)**2*za(i3,i5)*zb(i4,i3)*zb(i5,i2))/
     -   (4.*za(i1,i4)*za(i1,i5)*za(i2,i5)**2*zb(i3,i1)*
     -     zb(i5,i3)) + (za(i2,i3)*zb(i4,i3)**2*zb(i5,i2))/
     -   (4.*za(i2,i5)**2*zb(i3,i1)*zb(i3,i2)*zb(i5,i3)) + 
     -  (za(i1,i2)*za(i3,i4)*zb(i4,i2)*zb(i5,i3)**2)/
     -   (4.*za(i1,i4)*za(i2,i4)**2*zb(i3,i1)*zb(i3,i2)*
     -     zb(i4,i3)) + (za(i1,i3)*za(i3,i5)*zb(i5,i4))/
     -   (4.*za(i2,i3)*za(i2,i5)*za(i4,i5)*zb(i3,i2)) - 
     -  (za(i3,i5)**2*zb(i5,i4))/
     -   (4.*za(i2,i5)*za(i4,i5)**2*zb(i4,i1)) - 
     -  (za(i1,i5)**2*za(i3,i4)*zb(i5,i4))/
     -   (4.*za(i1,i4)*za(i2,i5)*za(i4,i5)**2*zb(i4,i3)) - 
     -  (za(i1,i5)*zb(i2,i1)*zb(i4,i2)*zb(i5,i4))/
     -   (4.*za(i4,i5)**2*zb(i3,i2)*zb(i4,i1)*zb(i4,i3)) + 
     -  (za(i1,i3)*za(i3,i4)*zb(i5,i4))/
     -   (4.*za(i1,i5)*za(i2,i4)*za(i4,i5)*zb(i5,i1)) - 
     -  (za(i1,i4)*za(i2,i3)*zb(i5,i4))/
     -   (4.*za(i2,i4)*za(i2,i5)*za(i4,i5)*zb(i5,i3)) - 
     -  (za(i1,i4)*zb(i5,i2)*zb(i5,i4))/
     -   (4.*za(i2,i4)*za(i4,i5)*zb(i3,i2)*zb(i5,i3)) + 
     -  (za(i3,i4)*za(i3,i5)*zb(i5,i4)**2)/
     -   (4.*za(i2,i3)*za(i2,i5)*za(i4,i5)*zb(i3,i2)*zb(i5,i1))
     -    + (zb(i4,i2)*zb(i5,i4)**2)/
     -   (4.*za(i2,i5)*zb(i3,i2)*zb(i4,i3)*zb(i5,i1)) - 
     -  (zb(i2,i1)*zb(i4,i2)*zb(i5,i4)**2)/
     -   (4.*za(i4,i5)*zb(i3,i2)*zb(i4,i1)*zb(i4,i3)*zb(i5,i1))
      
      amp_qqbggga_lc_gMHV=Bubs+Boxes+Vpole+Rat
!      Write(6,*) 'CC  = ',amp_qqbggga_lc_gMHV*im
!
!      write(6,*) 'Boxes =',Boxes
!      write(6,*) 'Bubs = ',Bubs
!      write(6,*) 'Vpole = ',Vpole
!      write(6,*) 'RAT - CR = ',Rat
!      pause
      return
      end


      !========SUB LEADING COLOR AMPLITUDE
      function amp_qqbggga_slc_gMHV(i1,i2,i3,i4,i5,za,zb)
      implicit none
      include 'types.f'
      complex(dp) :: amp_qqbggga_slc_gMHV
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
      complex(dp) :: Boxint(6),boxcoeff(6)
      integer i 

!------ boxcoeffs and ints
      ALO=-(za(i1,i3)**3/(za(i1,i4)*za(i1,i5)*za(i2,i5)*za(i3,i4)))

       
!-----pole expansion
    
      Vpole =  (-3._dp - epinv*epinv2
     &     + epinv*(-1.5 - lnrat(musq,-s(i1,i2))) - 
     &     lnrat(musq,-s(i1,i2))**2/2.
     &     - (3*lnrat(musq,-s(i3,i4)))/2._dp)
      Vpole=Vpole*ALO


!      write(6,*) ALO

      
      Boxint(1)=Lsm1(-s(i3,i5),-s(i1,i2),-s(i3,i4),-s(i1,i2))
      Boxint(2)=Lsm1(-s(i3,i5),-s(i1,i4),-s(i2,i3),-s(i1,i4))
      Boxint(3)=Lsm1(-s(i3,i4),-s(i1,i5),-s(i2,i3),-s(i1,i5))
      Boxint(4)=Lsm1(-s(i1,i5),-s(i3,i4),-s(i1,i2),-s(i3,i4))
      Boxint(5)=Lsm1(-s(i1,i2),-s(i3,i4),-s(i2,i5),-s(i3,i4))
      Boxint(6)=Lsm1(-s(i1,i2),-s(i3,i5),-s(i1,i4),-s(i3,i5))

!-boxes
      boxcoeff(1)=-((za(i1,i5)**2*za(i3,i4)**2)
     &     /(za(i1,i4)*za(i2,i5)*za(i4,i5)**3))

      boxcoeff(2)=-((za(i1,i2)**2*za(i3,i5)**2)
     &     /(za(i1,i4)*za(i2,i5)**3*za(i4,i5)))

      boxcoeff(3)= (za(i1,i2)**2*za(i3,i4)**2)/
     -     (za(i1,i5)*za(i2,i4)**3*za(i4,i5))

      boxcoeff(4)= (za(i1,i2)**2*za(i3,i5)**3)/
     -     (za(i1,i5)*za(i2,i5)**3*za(i3,i4)*za(i4,i5))

      boxcoeff(5)= (za(i1,i3)**3/
     -     (za(i1,i4)*za(i1,i5)*za(i2,i5)*za(i3,i4)))

      boxcoeff(6)=  -((za(i1,i2)**2*za(i3,i4)**2)/
     -    (za(i1,i4)*za(i2,i4)**2*za(i2,i5)*
     -     za(i4,i5)))


      !------Boxes
      Boxes=czip
      do i=1,6
         Boxes=Boxes+boxcoeff(i)*Boxint(i)
      enddo
!-----sub leading completed bubbles

      
      Bubs=   lnrat(-s(i3,i5),-s(i3,i4))*
     -   ((za(i1,i2)*za(i1,i3)*za(i2,i3))/
     -      (za(i1,i4)*za(i2,i4)*za(i2,i5)**2) + 
     -     (za(i1,i3)*za(i3,i4))/
     -      (za(i2,i4)*za(i4,i5)**2)) + 
     -  lnrat(-s(i1,i4),-s(i3,i4))*
     -   (-za(i1,i3)**2/
     -      (2.*za(i1,i4)*za(i2,i5)*za(i4,i5)) - 
     -     (za(i1,i2)*za(i1,i3)*za(i3,i5))/
     -      (za(i1,i4)*za(i2,i5)**2*za(i4,i5))) + 
     -  (lnrat(-s(i1,i5),-s(i3,i4))*za(i1,i3)*
     -     (za(i1,i3)*za(i2,i4) + 
     -       2*za(i1,i2)*za(i3,i4)))/
     -   (2.*za(i1,i5)*za(i2,i4)**2*za(i4,i5)) + 
     -  (L0(-s(i1,i4),-s(i3,i5))*za(i1,i2)*
     -     za(i2,i3)**2*zb(i4,i2))/
     -   (s(i3,i5)*za(i2,i4)*za(i2,i5)**2) + 
     -  (L0(-s(i1,i2),-s(i3,i5))*za(i1,i2)*
     -     za(i3,i4)**2*zb(i4,i2))/
     -   (s(i3,i5)*za(i2,i4)*za(i4,i5)**2) - 
     -  (L0(-s(i2,i3),-s(i1,i5))*
     -     (-(za(i1,i3)**2*za(i2,i4)**2) + 
     -       za(i1,i2)**2*za(i3,i4)**2)*zb(i4,i2))/
     -   (s(i1,i5)*za(i1,i5)*za(i2,i4)**2*za(i4,i5))
     -    - (L1(-s(i2,i3),-s(i1,i5))*za(i1,i4)**2*
     -     za(i2,i3)**2*zb(i4,i2)**2)/
     -   (2.*s(i1,i5)**2*za(i1,i5)*za(i2,i4)*
     -     za(i4,i5)) - 
     -  (L0(-s(i1,i2),-s(i3,i4))*za(i1,i2)*
     -     za(i3,i5)**2*zb(i5,i2))/
     -   (s(i3,i4)*za(i2,i5)*za(i4,i5)**2) + 
     -  (2*L0(-s(i1,i2),-s(i3,i4))*za(i1,i2)*
     -     za(i2,i3)*za(i3,i5)**2*zb(i5,i2))/
     -   (s(i3,i4)*za(i2,i5)**2*za(i3,i4)*za(i4,i5))
     -    - (L0(-s(i2,i3),-s(i1,i4))*
     -     (za(i1,i3)**2*za(i2,i5)**2 - 
     -       za(i1,i2)**2*za(i3,i5)**2)*zb(i5,i2))/
     -   (s(i1,i4)*za(i1,i4)*za(i2,i5)**2*za(i4,i5))
     -    + (L1(-s(i1,i5),-s(i3,i4))*za(i1,i2)*
     -     za(i1,i5)*za(i2,i3)**3*zb(i2,i1)*
     -     zb(i5,i2))/
     -   (2.*s(i3,i4)**2*za(i2,i4)*za(i2,i5)**2*
     -     za(i3,i4)) + 
     -  (L1(-s(i1,i2),-s(i3,i4))*za(i1,i2)**2*
     -     za(i3,i5)**3*zb(i5,i1)*zb(i5,i2))/
     -   (s(i3,i4)**2*za(i2,i5)**2*za(i3,i4)*
     -     za(i4,i5)) + 
     -  (L1(-s(i2,i3),-s(i1,i4))*za(i1,i5)**2*
     -     za(i2,i3)**2*zb(i5,i2)**2)/
     -   (2.*s(i1,i4)**2*za(i1,i4)*za(i2,i5)*
     -     za(i4,i5)) + 
     -  (L0(-s(i1,i5),-s(i3,i4))*
     -     ((za(i1,i2)*za(i2,i3)**2*zb(i5,i2))/
     -        (za(i2,i4)**2*za(i2,i5)) + 
     -       (za(i1,i5)*za(i2,i3)**3*zb(i5,i2))/
     -        (2.*za(i2,i4)*za(i2,i5)**2*za(i3,i4))
     -        + (2*za(i1,i2)*za(i2,i3)**2*za(i3,i5)*
     -          zb(i5,i2))/
     -        (za(i2,i4)*za(i2,i5)**2*za(i3,i4))))/
     -     s(i3,i4)

!---- rational, needs simp but is completed
      Rat= za(i1,i3)**3/
     -   (4.*za(i1,i4)*za(i1,i5)*za(i2,i5)*za(i3,i4)) + 
     -  (za(i1,i4)**2*za(i2,i3)*zb(i4,i2))/
     -   (4.*za(i1,i5)*za(i2,i4)**2*za(i4,i5)*zb(i3,i2)) - 
     -  (za(i1,i2)*za(i1,i3)*za(i2,i3)*zb(i4,i2))/
     -   (4.*za(i1,i5)*za(i2,i4)*za(i2,i5)*za(i3,i4)*zb(i4,i3)) + 
     -  (za(i1,i2)**3*za(i3,i4)*zb(i4,i2))/
     -   (4.*za(i1,i4)*za(i1,i5)*za(i2,i4)**2*za(i2,i5)*zb(i4,i3))
     -   - (za(i1,i3)*za(i1,i4)*za(i2,i3)*zb(i4,i2))/
     -   (4.*za(i1,i5)**2*za(i2,i4)*za(i4,i5)*zb(i5,i1)) + 
     -  (za(i1,i3)*za(i1,i4)*zb(i4,i2)**2)/
     -   (4.*za(i1,i5)**2*za(i4,i5)*zb(i3,i2)*zb(i5,i1)) - 
     -  (za(i1,i2)**2*za(i3,i4)*zb(i4,i2)**2)/
     -   (4.*za(i1,i5)**2*za(i2,i4)*za(i4,i5)*zb(i4,i3)*zb(i5,i1))
     -   + (za(i1,i2)*za(i3,i5)**3*zb(i5,i2))/
     -   (4.*za(i1,i5)*za(i2,i5)**2*za(i3,i4)*za(i4,i5)*zb(i2,i1))
     -   + (za(i1,i3)*zb(i5,i2))/
     -   (4.*za(i1,i4)*za(i4,i5)*zb(i3,i1)) + 
     -  (za(i1,i3)*za(i1,i5)*za(i2,i3)*zb(i5,i2))/
     -   (4.*za(i1,i4)**2*za(i2,i5)*za(i4,i5)*zb(i4,i1)) - 
     -  (za(i1,i2)*za(i2,i3)**2*za(i3,i5)*zb(i5,i2))/
     -   (4.*za(i2,i4)*za(i2,i5)**2*za(i3,i4)**2*zb(i4,i3)) - 
     -  (za(i1,i2)*za(i2,i3)*za(i3,i5)**2*zb(i5,i2))/
     -   (za(i2,i5)**2*za(i3,i4)**2*za(i4,i5)*zb(i4,i3)) - 
     -  (za(i1,i2)*za(i1,i5)*za(i2,i3)**2*zb(i2,i1)*zb(i5,i2))/
     -   (4.*za(i2,i4)*za(i2,i5)**2*za(i3,i4)**2*zb(i3,i2)*
     -     zb(i4,i3)) + (za(i1,i2)*za(i1,i5)**2*zb(i4,i1)*
     -     zb(i5,i2))/
     -   (4.*za(i1,i4)*za(i2,i5)**2*za(i4,i5)*zb(i3,i2)*zb(i4,i3))
     -   + (za(i1,i2)*za(i2,i3)*za(i3,i5)*zb(i4,i2)*zb(i5,i2))/
     -   (4.*za(i2,i5)**2*za(i3,i4)**2*zb(i3,i2)*zb(i4,i3)) - 
     -  (za(i1,i2)*za(i3,i5)*zb(i4,i2)**2*zb(i5,i2))/
     -   (4.*za(i1,i5)*za(i2,i5)*za(i3,i4)*zb(i2,i1)*zb(i3,i2)*
     -     zb(i4,i3)) - (za(i1,i5)*za(i2,i3)*zb(i4,i3)*zb(i5,i2))/
     -   (4.*za(i1,i4)*za(i2,i5)**2*zb(i3,i1)*zb(i3,i2)) + 
     -  (za(i2,i3)**3*zb(i5,i2))/
     -   (4.*za(i2,i4)*za(i2,i5)**2*za(i3,i4)*zb(i5,i1)) + 
     -  (za(i1,i2)*zb(i4,i1)**3*zb(i5,i2))/
     -   (4.*za(i2,i5)**2*zb(i2,i1)*zb(i3,i1)*zb(i4,i3)*zb(i5,i1))
     -   - (za(i1,i2)*zb(i4,i2)**2*zb(i5,i2))/
     -   (4.*za(i1,i5)*za(i4,i5)*zb(i3,i2)*zb(i4,i3)*zb(i5,i1)) - 
     -  (za(i1,i5)*za(i2,i3)**2*zb(i5,i2)**2)/
     -   (4.*za(i1,i4)**2*za(i2,i5)*za(i4,i5)*zb(i4,i1)*zb(i5,i1))
     -   + (za(i1,i2)**2*za(i1,i3)*zb(i3,i2)*zb(i5,i1))/
     -   (4.*za(i1,i4)*za(i1,i5)*za(i2,i4)*za(i2,i5)*zb(i3,i1)*
     -     zb(i5,i3)) - (za(i1,i2)*za(i1,i3)*zb(i4,i3)*zb(i5,i1))/
     -   (4.*za(i1,i5)*za(i2,i4)*za(i2,i5)*zb(i3,i1)*zb(i5,i3)) + 
     -  (za(i1,i2)*zb(i4,i1)*zb(i5,i2))/
     -   (4.*za(i2,i5)*za(i4,i5)*zb(i3,i1)*zb(i5,i3)) + 
     -  (za(i1,i2)*za(i3,i4)*zb(i4,i2)*zb(i5,i3)**2)/
     -   (4.*za(i1,i4)*za(i2,i4)**2*zb(i3,i1)*zb(i3,i2)*zb(i4,i3))
     -   + (za(i1,i3)**2*zb(i5,i4))/
     -   (4.*za(i1,i4)*za(i2,i5)*za(i3,i4)*zb(i4,i3)) + 
     -  (za(i1,i3)*za(i3,i5)*zb(i5,i4))/
     -   (2.*za(i2,i5)*za(i3,i4)*za(i4,i5)*zb(i4,i3)) + 
     -  (za(i1,i3)*za(i3,i4)*zb(i5,i4))/
     -   (4.*za(i1,i5)*za(i2,i4)*za(i4,i5)*zb(i5,i1)) - 
     -  (za(i1,i5)*za(i3,i5)*zb(i5,i2)*zb(i5,i4))/
     -   (4.*za(i1,i4)*za(i2,i5)*za(i4,i5)*zb(i3,i2)*zb(i4,i1)) + 
     -  (za(i3,i4)*zb(i5,i2)*zb(i5,i4))/
     -   (4.*za(i2,i4)*za(i4,i5)*zb(i3,i2)*zb(i5,i1)) - 
     -  (za(i2,i3)*za(i3,i5)*zb(i5,i2)*zb(i5,i4))/
     -   (4.*za(i2,i5)*za(i3,i4)*za(i4,i5)*zb(i4,i3)*zb(i5,i1)) - 
     -  (zb(i4,i2)*zb(i5,i2)*zb(i5,i4))/
     -   (4.*za(i4,i5)*zb(i3,i2)*zb(i4,i3)*zb(i5,i1)) - 
     -  (za(i1,i2)*za(i1,i4)*zb(i4,i1)*zb(i5,i4))/
     -   (4.*za(i2,i4)*za(i2,i5)*za(i4,i5)*zb(i4,i3)*zb(i5,i3)) + 
     -  (za(i1,i2)*za(i2,i3)*zb(i5,i2)**2*zb(i5,i4))/
     -   (4.*za(i1,i4)**2*za(i2,i5)*zb(i4,i1)*zb(i5,i1)*zb(i5,i3))
     -   - (za(i1,i2)*zb(i5,i4)**2)/
     -   (4.*za(i2,i4)*za(i2,i5)*zb(i4,i3)*zb(i5,i3)) - 
     -  (zb(i5,i2)**2*zb(i5,i4)**2)/
     -   (4.*za(i1,i4)*zb(i3,i2)*zb(i4,i1)*zb(i5,i1)*zb(i5,i3))

      amp_qqbggga_slc_gMHV=-Bubs+Boxes+Vpole+Rat
!      write(6,*) 'CC  = ',amp_qqbggga_slc_gMHV*im

!      write(6,*) 'Boxes =',Boxes
!      write(6,*) 'Bubs = ',Bubs
!      write(6,*) 'Vpole = ',Vpole
!      write(6,*) 'RAT - CR = ',Rat
!      pause
      return
      end
