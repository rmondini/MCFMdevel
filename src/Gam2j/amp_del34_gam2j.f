!---  routine which provides the amplitude proportional to del(34)del(i1i2),
!==== otherwise known as A_{5;3}
!=====notation is q(i1)+qb(i2)+g(i3)+g(i4) +gamma(i5)
!==== cw april 2016

      function amp_qqbggga_del34_gaMHV(i1,i2,i3,i4,i5,za,zb)
      implicit none
      include 'types.f'
      complex(dp) :: amp_qqbggga_del34_gaMHV
      include 'constants.f'
      include 'mxpart.f' 
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'epinv.f'
      include 'scale.f'
      integer i1,i2,i3,i4,i5
      integer Nbox,Ntri,Nbub
      parameter(Nbox=3,Ntri=6,Nbub=2)
      complex(dp) :: ALO,Vpole,rat,ALOs
      complex(dp) :: L0,L1,Bubs,Boxes,Lsm1,lnrat
      complex(dp) :: Boxint(18),boxcoeff(18)
      integer i 

      ALO=za(i1,i5)**2/(za(i2,i3)*za(i3,i4)*za(i4,i1))
      ALOs=za(i1,i5)**2/(za(i2,i4)*za(i4,i3)*za(i3,i1))

      Vpole=(za(i1,i5)**2*((-lnrat(musq,-s(i1,i2))**2 + 
     -         lnrat(musq,-s(i1,i3))**2 + lnrat(musq,-s(i2,i4))**2 - 
     -         lnrat(musq,-s(i3,i4))**2)/(za(i1,i4)*za(i2,i3)) + 
     -      2*epinv*((-lnrat(musq,-s(i1,i2)) + lnrat(musq,-s(i1,i3)) + 
     -            lnrat(musq,-s(i2,i4)) - lnrat(musq,-s(i3,i4)))/
     -          (za(i1,i4)*za(i2,i3)) + 
     -         (lnrat(musq,-s(i1,i2)) - lnrat(musq,-s(i1,i4)) - 
     -            lnrat(musq,-s(i2,i3)) + lnrat(musq,-s(i3,i4)))/
     -          (za(i1,i3)*za(i2,i4))) + 
     -      (lnrat(musq,-s(i1,i2))**2 - lnrat(musq,-s(i1,i4))**2 - 
     -         lnrat(musq,-s(i2,i3))**2 + lnrat(musq,-s(i3,i4))**2)/
     -     (za(i1,i3)*za(i2,i4))))/(2.*za(i3,i4))

      Boxint(1)=Lsm1(-s(i2,i5),-s(i3,i4),-s(i1,i5),-s(i3,i4))
      Boxint(2)=Lsm1(-s(i3,i5),-s(i1,i2),-s(i4,i5),-s(i1,i2))
      Boxint(3)=Lsm1(-s(i2,i5),-s(i1,i3),-s(i4,i5),-s(i1,i3))
      Boxint(4)=Lsm1(-s(i2,i5),-s(i1,i4),-s(i3,i5),-s(i1,i4))
      Boxint(5)=Lsm1(-s(i2,i4),-s(i1,i5),-s(i2,i3),-s(i1,i5))
      Boxint(6)=Lsm1(-s(i3,i4),-s(i1,i5),-s(i2,i3),-s(i1,i5))
      Boxint(7)=Lsm1(-s(i3,i4),-s(i1,i5),-s(i2,i4),-s(i1,i5))
      Boxint(8)=Lsm1(-s(i1,i5),-s(i2,i3),-s(i4,i5),-s(i2,i3))
      Boxint(9)=Lsm1(-s(i1,i5),-s(i2,i4),-s(i3,i5),-s(i2,i4))
      Boxint(10)=Lsm1(-s(i1,i4),-s(i2,i5),-s(i1,i3),-s(i2,i5))
      Boxint(11)=Lsm1(-s(i3,i4),-s(i2,i5),-s(i1,i3),-s(i2,i5))
      Boxint(12)=Lsm1(-s(i3,i4),-s(i2,i5),-s(i1,i4),-s(i2,i5))
      Boxint(13)=Lsm1(-s(i1,i3),-s(i4,i5),-s(i1,i2),-s(i4,i5))
      Boxint(14)=Lsm1(-s(i1,i2),-s(i3,i5),-s(i2,i4),-s(i3,i5))
      Boxint(15)=Lsm1(-s(i1,i4),-s(i3,i5),-s(i1,i2),-s(i3,i5))
      Boxint(16)=Lsm1(-s(i1,i2),-s(i4,i5),-s(i2,i3),-s(i4,i5))
      Boxint(17)=Lsm1(-s(i1,i3),-s(i4,i5),-s(i2,i3),-s(i4,i5))
      Boxint(18)=Lsm1(-s(i2,i4),-s(i3,i5),-s(i1,i4),-s(i3,i5))

      boxcoeff(1)=-((za(i1,i2)*za(i1,i5)**2)
     &     /(za(i1,i3)*za(i1,i4)*za(i2,i3)*za(i2,i4)))

!      write(6,*) boxcoeff(1)*s(i2,i5)*s(i1,i5)/2._dp*im
      
      boxcoeff(2)= (za(i1,i4)**3*za(i2,i3)*za(i3,i5)**2 - 
     -    za(i1,i3)**3*za(i2,i4)*za(i4,i5)**2)/
     -     (za(i1,i3)*za(i1,i4)*za(i2,i3)*za(i2,i4)*za(i3,i4)**3)

!      write(6,*) boxcoeff(2)*s(i3,i5)*s(i4,i5)/2._dp*im
      
      boxcoeff(3)= (za(i1,i2)**2*za(i4,i5)**2)/
     -     (za(i1,i4)*za(i2,i3)*za(i2,i4)**2*za(i3,i4))
!      write(6,*) boxcoeff(3)*s(i2,i5)*s(i4,i5)/2._dp*im

      boxcoeff(4)= -((za(i1,i2)**2*za(i3,i5)**2)/
     -     (za(i1,i3)*za(i2,i3)**2*za(i2,i4)*za(i3,i4)))
!      write(6,*) boxcoeff(4)*s(i2,i5)*s(i3,i5)/2._dp*im

      boxcoeff(5)=(za(i1,i2)*za(i1,i5)**2)
     &     /(za(i1,i3)*za(i1,i4)*za(i2,i3)*za(i2,i4))
 !     write(6,*) boxcoeff(5)*s(i2,i4)*s(i2,i3)/2._dp*im

      boxcoeff(6)=-(za(i1,i5)**2/(za(i1,i4)*za(i2,i3)*za(i3,i4)))
 !     write(6,*) boxcoeff(6)*s(i3,i4)*s(i2,i3)/2._dp*im

      boxcoeff(7)=za(i1,i5)**2/(za(i1,i3)*za(i2,i4)*za(i3,i4))
!      write(6,*) boxcoeff(7)*s(i3,i4)*s(i2,i4)/2._dp*im

      boxcoeff(8)=-za(i1,i5)**2/(za(i1,i3)*za(i2,i4)*za(i3,i4))
!      write(6,*) boxcoeff(8)*s(i4,i5)*s(i1,i5)/2._dp*im

      boxcoeff(9)=za(i1,i5)**2/(za(i1,i4)*za(i2,i3)*za(i3,i4))
!      write(6,*) boxcoeff(9)*s(i3,i5)*s(i1,i5)/2._dp*im

      boxcoeff(10)= (-((za(i1,i4)**3*za(i3,i5)**2)/za(i2,i4)) + 
     -     (za(i1,i3)**3*za(i4,i5)**2)/za(i2,i3))
     &     /(za(i1,i3)*za(i1,i4)*za(i3,i4)**3)
!      write(6,*) boxcoeff(10)*s(i4,i1)*s(i1,i3)/2._dp*im

      boxcoeff(11)=za(i1,i5)**2/(za(i1,i3)*za(i2,i4)*za(i3,i4))
!      write(6,*) boxcoeff(11)*s(i3,i4)*s(i1,i3)/2._dp*im

      boxcoeff(12)=-(za(i1,i5)**2/(za(i1,i4)*za(i2,i3)*za(i3,i4)))
!      write(6,*) boxcoeff(12)*s(i3,i4)*s(i1,i4)/2._dp*im

      boxcoeff(13)=((za(i1,i2)**2*za(i3,i5)**2)
     &     /(za(i1,i3)*za(i2,i3)**2*za(i2,i4)*za(i3,i4)))
!      write(6,*) boxcoeff(13)*s(i1,i3)*s(i1,i2)/2._dp*im

      boxcoeff(14)=za(i1,i5)**2/(za(i1,i3)*za(i2,i4)*za(i3,i4))
      
      boxcoeff(15)=-((za(i1,i2)**2*za(i4,i5)**2)
     &     /(za(i1,i4)*za(i2,i3)*za(i2,i4)**2*za(i3,i4)))

      boxcoeff(16)=-(za(i1,i5)**2/(za(i1,i4)*za(i2,i3)*za(i3,i4)))

      boxcoeff(17)=(za(i1,i2)*za(i1,i5)**2)
     &     /(za(i1,i3)*za(i1,i4)*za(i2,i3)*za(i2,i4))
      

      boxcoeff(18)=(za(i1,i2)*za(i1,i5)**2)
     &     /(za(i1,i3)*za(i1,i4)*za(i2,i3)*za(i2,i4))


      Boxes=czip
      do i=1,18
         Boxes=Boxes+boxcoeff(i)*boxint(i)
      enddo

      amp_qqbggga_del34_gaMHV=Vpole+Boxes
!      write(6,*) 'del 34',amp_qqbggga_del34_gaMHV*im
      return
      end

!==== amplitude for q(i1)- qb(i2)+ g(i3)+ g(i4)- g(i5)
      
      function amp_qqbggga_del34_gMHV(i1,i2,i3,i4,i5,za,zb)
      implicit none
      include 'types.f'
      complex(dp) :: amp_qqbggga_del34_gMHV
      include 'constants.f'
      include 'mxpart.f' 
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'epinv.f'
      include 'scale.f'
      integer i1,i2,i3,i4,i5
      integer Nbox,Ntri,Nbub
      parameter(Nbox=3,Ntri=6,Nbub=2)
      complex(dp) :: ALO,Vpole,rat,ALOs
      complex(dp) :: L0,L1,Bubs,Boxes,Lsm1,lnrat
      complex(dp) :: Boxint(18),boxcoeff(18)
      integer i


      ALO=-((za(i1,i4)**2*za(i2,i4))
     &     /(za(i1,i5)*za(i2,i3)*za(i2,i5)*za(i3,i4)))
      ALOs=za(i1,i4)**3/(za(i1,i3)*za(i1,i5)*za(i2,i5)*za(i3,i4))

      Vpole= (za(i1,i4)**2*((lnrat(musq,-s(i1,i2))**2 - 
     -         lnrat(musq,-s(i1,i4))**2 - lnrat(musq,-s(i2,i3))**2 + 
     -         lnrat(musq,-s(i3,i4))**2)*za(i1,i4)*za(i2,i3) + 
     -      (-lnrat(musq,-s(i1,i2))**2 + lnrat(musq,-s(i1,i3))**2 + 
     -         lnrat(musq,-s(i2,i4))**2 - lnrat(musq,-s(i3,i4))**2)*
     -       za(i1,i3)*za(i2,i4) + 
     -      2*epinv*((lnrat(musq,-s(i1,i2)) - 
     -            lnrat(musq,-s(i1,i4)) - lnrat(musq,-s(i2,i3)) + 
     -            lnrat(musq,-s(i3,i4)))*za(i1,i4)*za(i2,i3) + 
     -         (-lnrat(musq,-s(i1,i2)) + lnrat(musq,-s(i1,i3)) + 
     -            lnrat(musq,-s(i2,i4)) - lnrat(musq,-s(i3,i4)))*
     -          za(i1,i3)*za(i2,i4))))/
     -  (2.*za(i1,i3)*za(i1,i5)*za(i2,i3)*za(i2,i5)*za(i3,i4))


      Boxint(1)=Lsm1(-s(i3,i4),-s(i1,i2),-s(i4,i5),-s(i1,i2))
      Boxint(2)=Lsm1(-s(i2,i5),-s(i1,i3),-s(i2,i4),-s(i1,i3))
      Boxint(3)=Lsm1(-s(i2,i4),-s(i1,i3),-s(i4,i5),-s(i1,i3))
      Boxint(4)=Lsm1(-s(i2,i3),-s(i1,i4),-s(i2,i5),-s(i1,i4))
      Boxint(5)=Lsm1(-s(i2,i4),-s(i1,i5),-s(i2,i3),-s(i1,i5))
      Boxint(6)=Lsm1(-s(i3,i4),-s(i1,i5),-s(i2,i3),-s(i1,i5))
      Boxint(7)=Lsm1(-s(i3,i4),-s(i1,i5),-s(i2,i4),-s(i1,i5))
      Boxint(8)=Lsm1(-s(i1,i5),-s(i2,i3),-s(i1,i4),-s(i2,i3))
      Boxint(9)=Lsm1(-s(i1,i4),-s(i2,i3),-s(i4,i5),-s(i2,i3))
      Boxint(10)=Lsm1(-s(i1,i5),-s(i2,i4),-s(i1,i3),-s(i2,i4))
      Boxint(11)=Lsm1(-s(i1,i4),-s(i2,i5),-s(i1,i3),-s(i2,i5))
      Boxint(12)=Lsm1(-s(i3,i4),-s(i2,i5),-s(i1,i3),-s(i2,i5))
      Boxint(13)=Lsm1(-s(i3,i4),-s(i2,i5),-s(i1,i4),-s(i2,i5))
      Boxint(14)=Lsm1(-s(i1,i2),-s(i3,i4),-s(i1,i5),-s(i3,i4))
      Boxint(15)=Lsm1(-s(i1,i2),-s(i3,i4),-s(i2,i5),-s(i3,i4))
      Boxint(16)=Lsm1(-s(i1,i2),-s(i4,i5),-s(i1,i3),-s(i4,i5))
      Boxint(17)=Lsm1(-s(i1,i2),-s(i4,i5),-s(i2,i3),-s(i4,i5))
      Boxint(18)=Lsm1(-s(i1,i3),-s(i4,i5),-s(i2,i3),-s(i4,i5))


      boxcoeff(1)=  (za(i1,i5)**3*za(i2,i3)*za(i3,i4)**2 - 
     -    za(i1,i3)**3*za(i2,i5)*za(i4,i5)**2)/
     -     (za(i1,i3)*za(i1,i5)*za(i2,i3)*za(i2,i5)*za(i3,i5)**3)

!      write(6,*) 'boxcoeff 1',boxcoeff(1)*s(i3,i4)*s(i4,i5)/2.*im

      boxcoeff(2)=(za(i1,i4)**2*za(i4,i5))
     &     /(za(i1,i5)*za(i2,i5)*za(i3,i4)*za(i3,i5))
      
 !     write(6,*) 'boxcoeff 2',boxcoeff(2)*s(i2,i4)*s(i2,i5)/2.*im
      
      boxcoeff(3)=(za(i1,i2)**2*za(i4,i5)**2)/
     -     (za(i1,i5)*za(i2,i3)*za(i2,i5)**2*za(i3,i5))
 !     write(6,*) 'boxcoeff 3',boxcoeff(3)*s(i2,i4)*s(i4,i5)/2.*im
      boxcoeff(4)=-(za(i1,i4)**3
     &     /(za(i1,i3)*za(i1,i5)*za(i2,i5)*za(i3,i4)))
 !     write(6,*) 'boxcoeff 4',boxcoeff(4)*s(i2,i3)*s(i2,i5)/2.*im
      boxcoeff(5)=za(i1,i4)**2/(za(i1,i5)*za(i2,i3)*za(i3,i5))
 !     write(6,*) 'boxcoeff 5',boxcoeff(5)*s(i2,i3)*s(i2,i4)/2.*im
     
      boxcoeff(6)= -((za(i1,i4)**2*za(i2,i4))/
     -     (za(i1,i5)*za(i2,i3)*za(i2,i5)*za(i3,i4)))
!      write(6,*) 'boxcoeff 6',boxcoeff(6)*s(i2,i3)*s(i3,i4)/2.*im
     
      boxcoeff(7)= (-(za(i1,i2)**2*za(i2,i5)*za(i3,i4)**3) + 
     -    za(i1,i3)**2*za(i2,i4)**3*za(i3,i5))/
     -     (za(i1,i5)*za(i2,i3)**3*za(i2,i5)*za(i3,i4)*za(i3,i5))
!      write(6,*) 'boxcoeff 7',boxcoeff(7)*s(i3,i4)*s(i2,i4)/2.*im
        boxcoeff(8)=-((za(i1,i4)**2*za(i4,i5))/
     -     (za(i1,i5)*za(i2,i5)*za(i3,i4)*za(i3,i5)))
!      write(6,*) 'boxcoeff 8',boxcoeff(8)*s(i1,i5)*s(i1,i4)/2.*im
      boxcoeff(9)=-(za(i1,i4)**2/(za(i1,i3)*za(i2,i5)*za(i3,i5)))
!      write(6,*) 'boxcoeff 9',boxcoeff(9)*s(i4,i5)*s(i1,i4)/2.*im
      boxcoeff(10)= ((za(i1,i3)**2*za(i2,i4)*za(i4,i5)**2)/
     -     (za(i1,i5)*za(i2,i3)*za(i2,i5)*za(i3,i4)*za(i3,i5)**2))
 !     write(6,*) 'boxcoeff 10',boxcoeff(10)*s(i1,i5)*s(i1,i3)/2.*im
      boxcoeff(11)=-(za(i1,i4)**2/(za(i1,i3)*za(i2,i5)*za(i3,i5)))
      boxcoeff(12)=za(i1,i4)**3/(za(i1,i3)*za(i1,i5)*za(i2,i5)*za(i3,i4))

      boxcoeff(13)= -((za(i1,i4)**2*za(i4,i5))/
     -     (za(i1,i5)*za(i2,i5)*za(i3,i4)*za(i3,i5)))

      boxcoeff(14)=-((za(i1,i2)**2*za(i4,i5)**2)/
     -     (za(i1,i5)*za(i2,i3)*za(i2,i5)**2*za(i3,i5)))

      boxcoeff(15)=za(i1,i4)**2/(za(i1,i3)*za(i2,i5)*za(i3,i5))

      boxcoeff(16)=(za(i1,i2)**2*za(i3,i4)**2)/
     -     (za(i1,i3)*za(i2,i3)**2*za(i2,i5)*za(i3,i5))

      boxcoeff(17)=-(za(i1,i4)**2/(za(i1,i5)*za(i2,i3)*za(i3,i5)))
      boxcoeff(18)=(za(i1,i2)*za(i1,i4)**2)
     &     /(za(i1,i3)*za(i1,i5)*za(i2,i3)*za(i2,i5))
      
      Boxes=czip
      do i=1,18
         Boxes=Boxes+boxcoeff(i)*boxint(i)
      enddo

      amp_qqbggga_del34_gMHV=Vpole+Boxes
!      write(6,*) 'del 34',amp_qqbggga_del34_gMHV*im
      return
      end
