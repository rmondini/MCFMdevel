
      subroutine nloxnlofeval(yy,zz,logmh2omu2,nloxnloeps0) 
      implicit none 
      include 'types.f'
      include 'constants.f'
      include 'nf.f'  
      include 'ewcouple.f' 
      include 'masses.f'
      include 'scet_const.f'
      integer i
      real(dp) xx

      real(dp) a1mipart(5),a2mipart(5),a1renorm(5),a2renorm(5),a1infrared(5),a2infrared(5)

!====
      real(dp) interfpart2oneoneeps(5),interfpart2oneonenoeps(5)
      real(dp) interfpart2twotwo(5),interfpart2twoone(5),interfpart2onetwo(5)

      real(dp) interfpart3oneoneeps(5),interfpart3oneonenoeps(5)
      real(dp) interfpart3twotwo(5),interfpart3onetwo(5)

      real(dp) interfpart1final(5),interfpart2final(5),interfpart3final(5),interffinal(5)

!====
      real(dp) infraredpart1oneoneeps(5),infraredpart1oneonenoeps(5)
      real(dp) infraredpart1twotwo(5),infraredpart1onetwo(5)

      real(dp) infraredpart2oneoneeps(5),infraredpart2oneonenoeps(5)
      real(dp) infraredpart2twotwo(5),infraredpart2twoone(5),infraredpart2onetwo(5)

      real(dp) infraredpart3oneoneeps(5),infraredpart3oneonenoeps(5)
      real(dp) infraredpart3twotwo(5),infraredpart3twoone(5),infraredpart3onetwo(5)

      real(dp) infraredpart1final(5),infraredpart2final(5),infraredpart3final(5)

!====
      real(dp),intent(in)::yy,zz,logmh2omu2
      real(dp),intent(out)::nloxnloeps0

c--- initialize to zero
      nloxnloeps0=0._dp
      interffinal(:)=0._dp

      interfpart2oneoneeps(:)=0._dp
      interfpart2oneonenoeps(:)=0._dp
      interfpart2twotwo(:)=0._dp
      interfpart2onetwo(:)=0._dp
      interfpart2twoone(:)=0._dp

      interfpart3oneoneeps(:)=0._dp
      interfpart3oneonenoeps(:)=0._dp
      interfpart3twotwo(:)=0._dp
      interfpart3onetwo(:)=0._dp

      interfpart1final(:)=0._dp
      interfpart2final(:)=0._dp
      interfpart3final(:)=0._dp

      infraredpart1oneoneeps(:)=0._dp
      infraredpart1oneonenoeps(:)=0._dp
      infraredpart1twotwo(:)=0._dp
      infraredpart1onetwo(:)=0._dp

      infraredpart2oneoneeps(:)=0._dp
      infraredpart2oneonenoeps(:)=0._dp
      infraredpart2twotwo(:)=0._dp
      infraredpart2onetwo(:)=0._dp
      infraredpart2twoone(:)=0._dp

      infraredpart3oneoneeps(:)=0._dp
      infraredpart3oneonenoeps(:)=0._dp
      infraredpart3twotwo(:)=0._dp
      infraredpart3onetwo(:)=0._dp
      infraredpart3twoone(:)=0._dp

      infraredpart1final(:)=0._dp
      infraredpart2final(:)=0._dp
      infraredpart3final(:)=0._dp

!==================================================================

c---- value for x
      xx = one-yy-zz

!==================================================================

!==== numerical evaluation of the MI part, infrared part, and renormalization part of A1,A2 after eps-expansion up to eps^2

      call acoeffnlomiparteval(yy,zz,logmh2omu2,a1mipart,a2mipart)
      call acoeffnloinfraredeval(yy,zz,logmh2omu2,a1infrared,a2infrared)
      call acoeffnlorenormeval(yy,zz,logmh2omu2,a1renorm,a2renorm)

!==== construction of all required pieces of the interference

!==== infrared part 1

      infraredpart1oneonenoeps(1)=four*yy*zz*(a1infrared(1)**2)
      infraredpart1oneonenoeps(2)=four*yy*zz*(two*a1infrared(1)*a1infrared(2))
      infraredpart1oneonenoeps(3)=four*yy*zz*(two*a1infrared(1)*a1infrared(3)+a1infrared(2)**2)
      infraredpart1oneonenoeps(4)=four*yy*zz*(two*a1infrared(1)*a1infrared(4)+two*a1infrared(2)*a1infrared(3))
      infraredpart1oneonenoeps(5)=four*yy*zz*(two*a1infrared(1)*a1infrared(5)+two*a1infrared(2)*a1infrared(4)+a1infrared(3)**2)

      do i=2,5
         infraredpart1oneoneeps(i)=-infraredpart1oneonenoeps(i-1)
      enddo

      infraredpart1twotwo(1)=two*xx**2*yy/zz*(a2infrared(1)**2)
      infraredpart1twotwo(2)=two*xx**2*yy/zz*(two*a2infrared(1)*a2infrared(2))
      infraredpart1twotwo(3)=two*xx**2*yy/zz*(two*a2infrared(1)*a2infrared(3)+a2infrared(2)**2)
      infraredpart1twotwo(4)=two*xx**2*yy/zz*(two*a2infrared(1)*a2infrared(4)+two*a2infrared(2)*a2infrared(3))
      infraredpart1twotwo(5)=two*xx**2*yy/zz*(two*a2infrared(1)*a2infrared(5)+two*a2infrared(2)*a2infrared(4)+a2infrared(3)**2)

      infraredpart1onetwo(1)=-two*xx*yy*(a1infrared(1)*a2infrared(1))
      infraredpart1onetwo(2)=-two*xx*yy*(a1infrared(1)*a2infrared(2)+a1infrared(2)*a2infrared(1))
      infraredpart1onetwo(3)=-two*xx*yy*(a1infrared(1)*a2infrared(3)+a1infrared(2)*a2infrared(2)+a1infrared(3)*a2infrared(1))
      infraredpart1onetwo(4)=-two*xx*yy*(a1infrared(1)*a2infrared(4)+a1infrared(2)*a2infrared(3)+a1infrared(3)*a2infrared(2)+a1infrared(4)*a2infrared(1))
      infraredpart1onetwo(5)=-two*xx*yy*(a1infrared(1)*a2infrared(5)+a1infrared(2)*a2infrared(4)+a1infrared(3)*a2infrared(3)+a1infrared(4)*a2infrared(2)+a1infrared(5)*a2infrared(1))

      do i=1,5
         infraredpart1final(i)=-(infraredpart1oneonenoeps(i)+infraredpart1oneoneeps(i)+infraredpart1twotwo(i)+two*infraredpart1onetwo(i))
      enddo

!==== infrared part 2

      infraredpart2oneonenoeps(1)=four*yy*zz*(a1mipart(1)*a1infrared(1))
      infraredpart2oneonenoeps(2)=four*yy*zz*(a1mipart(1)*a1infrared(2)+a1mipart(2)*a1infrared(1))
      infraredpart2oneonenoeps(3)=four*yy*zz*(a1mipart(1)*a1infrared(3)+a1mipart(2)*a1infrared(2)+a1mipart(3)*a1infrared(1))
      infraredpart2oneonenoeps(4)=four*yy*zz*(a1mipart(1)*a1infrared(4)+a1mipart(2)*a1infrared(3)+a1mipart(3)*a1infrared(2)+a1mipart(4)*a1infrared(1))
      infraredpart2oneonenoeps(5)=four*yy*zz*(a1mipart(1)*a1infrared(5)+a1mipart(2)*a1infrared(4)+a1mipart(3)*a1infrared(3)+a1mipart(4)*a1infrared(2)+a1mipart(5)*a1infrared(1))

      do i=2,5
         infraredpart2oneoneeps(i)=-infraredpart2oneonenoeps(i-1)
      enddo

      infraredpart2twotwo(1)=two*xx**2*yy/zz*(a2mipart(1)*a2infrared(1))
      infraredpart2twotwo(2)=two*xx**2*yy/zz*(a2mipart(1)*a2infrared(2)+a2mipart(2)*a2infrared(1))
      infraredpart2twotwo(3)=two*xx**2*yy/zz*(a2mipart(1)*a2infrared(3)+a2mipart(2)*a2infrared(2)+a2mipart(3)*a2infrared(1))
      infraredpart2twotwo(4)=two*xx**2*yy/zz*(a2mipart(1)*a2infrared(4)+a2mipart(2)*a2infrared(3)+a2mipart(3)*a2infrared(2)+a2mipart(4)*a2infrared(1))
      infraredpart2twotwo(5)=two*xx**2*yy/zz*(a2mipart(1)*a2infrared(5)+a2mipart(2)*a2infrared(4)+a2mipart(3)*a2infrared(3)+a2mipart(4)*a2infrared(2)+a2mipart(5)*a2infrared(1))

      infraredpart2onetwo(1)=-two*xx*yy*(a1mipart(1)*a2infrared(1))
      infraredpart2onetwo(2)=-two*xx*yy*(a1mipart(1)*a2infrared(2)+a1mipart(2)*a2infrared(1))
      infraredpart2onetwo(3)=-two*xx*yy*(a1mipart(1)*a2infrared(3)+a1mipart(2)*a2infrared(2)+a1mipart(3)*a2infrared(1))
      infraredpart2onetwo(4)=-two*xx*yy*(a1mipart(1)*a2infrared(4)+a1mipart(2)*a2infrared(3)+a1mipart(3)*a2infrared(2)+a1mipart(4)*a2infrared(1))
      infraredpart2onetwo(5)=-two*xx*yy*(a1mipart(1)*a2infrared(5)+a1mipart(2)*a2infrared(4)+a1mipart(3)*a2infrared(3)+a1mipart(4)*a2infrared(2)+a1mipart(5)*a2infrared(1))

      infraredpart2twoone(1)=-two*xx*yy*(a2mipart(1)*a1infrared(1))
      infraredpart2twoone(2)=-two*xx*yy*(a2mipart(1)*a1infrared(2)+a2mipart(2)*a1infrared(1))
      infraredpart2twoone(3)=-two*xx*yy*(a2mipart(1)*a1infrared(3)+a2mipart(2)*a1infrared(2)+a2mipart(3)*a1infrared(1))
      infraredpart2twoone(4)=-two*xx*yy*(a2mipart(1)*a1infrared(4)+a2mipart(2)*a1infrared(3)+a2mipart(3)*a1infrared(2)+a2mipart(4)*a1infrared(1))
      infraredpart2twoone(5)=-two*xx*yy*(a2mipart(1)*a1infrared(5)+a2mipart(2)*a1infrared(4)+a2mipart(3)*a1infrared(3)+a2mipart(4)*a1infrared(2)+a2mipart(5)*a1infrared(1))

      do i=1,5
         infraredpart2final(i)=two*(infraredpart2oneonenoeps(i)+infraredpart2oneoneeps(i)+infraredpart2twotwo(i)+infraredpart2twoone(i)+infraredpart2onetwo(i))
      enddo

!==== infrared part 3

      infraredpart3oneonenoeps(1)=four*yy*zz*(a1renorm(1)*a1infrared(1))
      infraredpart3oneonenoeps(2)=four*yy*zz*(a1renorm(1)*a1infrared(2)+a1renorm(2)*a1infrared(1))
      infraredpart3oneonenoeps(3)=four*yy*zz*(a1renorm(1)*a1infrared(3)+a1renorm(2)*a1infrared(2)+a1renorm(3)*a1infrared(1))
      infraredpart3oneonenoeps(4)=four*yy*zz*(a1renorm(1)*a1infrared(4)+a1renorm(2)*a1infrared(3)+a1renorm(3)*a1infrared(2)+a1renorm(4)*a1infrared(1))
      infraredpart3oneonenoeps(5)=four*yy*zz*(a1renorm(1)*a1infrared(5)+a1renorm(2)*a1infrared(4)+a1renorm(3)*a1infrared(3)+a1renorm(4)*a1infrared(2)+a1renorm(5)*a1infrared(1))

      do i=2,5
         infraredpart3oneoneeps(i)=-infraredpart3oneonenoeps(i-1)
      enddo

      infraredpart3twotwo(1)=two*xx**2*yy/zz*(a2renorm(1)*a2infrared(1))
      infraredpart3twotwo(2)=two*xx**2*yy/zz*(a2renorm(1)*a2infrared(2)+a2renorm(2)*a2infrared(1))
      infraredpart3twotwo(3)=two*xx**2*yy/zz*(a2renorm(1)*a2infrared(3)+a2renorm(2)*a2infrared(2)+a2renorm(3)*a2infrared(1))
      infraredpart3twotwo(4)=two*xx**2*yy/zz*(a2renorm(1)*a2infrared(4)+a2renorm(2)*a2infrared(3)+a2renorm(3)*a2infrared(2)+a2renorm(4)*a2infrared(1))
      infraredpart3twotwo(5)=two*xx**2*yy/zz*(a2renorm(1)*a2infrared(5)+a2renorm(2)*a2infrared(4)+a2renorm(3)*a2infrared(3)+a2renorm(4)*a2infrared(2)+a2renorm(5)*a2infrared(1))

      infraredpart3onetwo(1)=-two*xx*yy*(a1renorm(1)*a2infrared(1))
      infraredpart3onetwo(2)=-two*xx*yy*(a1renorm(1)*a2infrared(2)+a1renorm(2)*a2infrared(1))
      infraredpart3onetwo(3)=-two*xx*yy*(a1renorm(1)*a2infrared(3)+a1renorm(2)*a2infrared(2)+a1renorm(3)*a2infrared(1))
      infraredpart3onetwo(4)=-two*xx*yy*(a1renorm(1)*a2infrared(4)+a1renorm(2)*a2infrared(3)+a1renorm(3)*a2infrared(2)+a1renorm(4)*a2infrared(1))
      infraredpart3onetwo(5)=-two*xx*yy*(a1renorm(1)*a2infrared(5)+a1renorm(2)*a2infrared(4)+a1renorm(3)*a2infrared(3)+a1renorm(4)*a2infrared(2)+a1renorm(5)*a2infrared(1))

      infraredpart3twoone(1)=-two*xx*yy*(a2renorm(1)*a1infrared(1))
      infraredpart3twoone(2)=-two*xx*yy*(a2renorm(1)*a1infrared(2)+a2renorm(2)*a1infrared(1))
      infraredpart3twoone(3)=-two*xx*yy*(a2renorm(1)*a1infrared(3)+a2renorm(2)*a1infrared(2)+a2renorm(3)*a1infrared(1))
      infraredpart3twoone(4)=-two*xx*yy*(a2renorm(1)*a1infrared(4)+a2renorm(2)*a1infrared(3)+a2renorm(3)*a1infrared(2)+a2renorm(4)*a1infrared(1))
      infraredpart3twoone(5)=-two*xx*yy*(a2renorm(1)*a1infrared(5)+a2renorm(2)*a1infrared(4)+a2renorm(3)*a1infrared(3)+a2renorm(4)*a1infrared(2)+a2renorm(5)*a1infrared(1))

      do i=1,5
         infraredpart3final(i)=two*(infraredpart3oneonenoeps(i)+infraredpart3oneoneeps(i)+infraredpart3twotwo(i)+infraredpart3twoone(i)+infraredpart3onetwo(i))
      enddo

      call applyRe(infraredpart3final)

!==== interference part 1

      interfpart1final(3)=2*pi**2*(2*xx**2+2*xx*(yy+zz)+(yy+zz)**2)*(11*ca+18*cf-4*nf*tr)**2/(9*yy*zz)
      interfpart1final(4)=-2*pi**2*(yy+zz)**2*(11*ca+18*cf-4*nf*tr)**2/(9*yy*zz)

!==== interference part 2

      interfpart2oneonenoeps(1)=four*yy*zz*(a1mipart(1)*a1renorm(1))
      interfpart2oneonenoeps(2)=four*yy*zz*(a1mipart(1)*a1renorm(2)+a1mipart(2)*a1renorm(1))
      interfpart2oneonenoeps(3)=four*yy*zz*(a1mipart(1)*a1renorm(3)+a1mipart(2)*a1renorm(2)+a1mipart(3)*a1renorm(1))
      interfpart2oneonenoeps(4)=four*yy*zz*(a1mipart(1)*a1renorm(4)+a1mipart(2)*a1renorm(3)+a1mipart(3)*a1renorm(2)+a1mipart(4)*a1renorm(1))
      interfpart2oneonenoeps(5)=four*yy*zz*(a1mipart(1)*a1renorm(5)+a1mipart(2)*a1renorm(4)+a1mipart(3)*a1renorm(3)+a1mipart(4)*a1renorm(2)+a1mipart(5)*a1renorm(1))

      do i=2,5
         interfpart2oneoneeps(i)=-interfpart2oneonenoeps(i-1)
      enddo

      interfpart2twotwo(1)=two*xx**2*yy/zz*(a2mipart(1)*a2renorm(1))
      interfpart2twotwo(2)=two*xx**2*yy/zz*(a2mipart(1)*a2renorm(2)+a2mipart(2)*a2renorm(1))
      interfpart2twotwo(3)=two*xx**2*yy/zz*(a2mipart(1)*a2renorm(3)+a2mipart(2)*a2renorm(2)+a2mipart(3)*a2renorm(1))
      interfpart2twotwo(4)=two*xx**2*yy/zz*(a2mipart(1)*a2renorm(4)+a2mipart(2)*a2renorm(3)+a2mipart(3)*a2renorm(2)+a2mipart(4)*a2renorm(1))
      interfpart2twotwo(5)=two*xx**2*yy/zz*(a2mipart(1)*a2renorm(5)+a2mipart(2)*a2renorm(4)+a2mipart(3)*a2renorm(3)+a2mipart(4)*a2renorm(2)+a2mipart(5)*a2renorm(1))

      interfpart2onetwo(1)=-two*xx*yy*(a1mipart(1)*a2renorm(1))
      interfpart2onetwo(2)=-two*xx*yy*(a1mipart(1)*a2renorm(2)+a1mipart(2)*a2renorm(1))
      interfpart2onetwo(3)=-two*xx*yy*(a1mipart(1)*a2renorm(3)+a1mipart(2)*a2renorm(2)+a1mipart(3)*a2renorm(1))
      interfpart2onetwo(4)=-two*xx*yy*(a1mipart(1)*a2renorm(4)+a1mipart(2)*a2renorm(3)+a1mipart(3)*a2renorm(2)+a1mipart(4)*a2renorm(1))
      interfpart2onetwo(5)=-two*xx*yy*(a1mipart(1)*a2renorm(5)+a1mipart(2)*a2renorm(4)+a1mipart(3)*a2renorm(3)+a1mipart(4)*a2renorm(2)+a1mipart(5)*a2renorm(1))

      interfpart2twoone(1)=-two*xx*yy*(a2mipart(1)*a1renorm(1))
      interfpart2twoone(2)=-two*xx*yy*(a2mipart(1)*a1renorm(2)+a2mipart(2)*a1renorm(1))
      interfpart2twoone(3)=-two*xx*yy*(a2mipart(1)*a1renorm(3)+a2mipart(2)*a1renorm(2)+a2mipart(3)*a1renorm(1))
      interfpart2twoone(4)=-two*xx*yy*(a2mipart(1)*a1renorm(4)+a2mipart(2)*a1renorm(3)+a2mipart(3)*a1renorm(2)+a2mipart(4)*a1renorm(1))
      interfpart2twoone(5)=-two*xx*yy*(a2mipart(1)*a1renorm(5)+a2mipart(2)*a1renorm(4)+a2mipart(3)*a1renorm(3)+a2mipart(4)*a1renorm(2)+a2mipart(5)*a1renorm(1))

      do i=1,5
         interfpart2final(i)=interfpart2oneonenoeps(i)+interfpart2oneoneeps(i)+interfpart2twotwo(i)+interfpart2twoone(i)+interfpart2onetwo(i)
      enddo

      call applyRe(interfpart2final)

!==== interference part 3

      interfpart3oneonenoeps(1)=four*yy*zz*(a1mipart(1)**2)
      interfpart3oneonenoeps(2)=four*yy*zz*(two*a1mipart(1)*a1mipart(2))
      interfpart3oneonenoeps(3)=four*yy*zz*(two*a1mipart(1)*a1mipart(3)+a1mipart(2)**2)
      interfpart3oneonenoeps(4)=four*yy*zz*(two*a1mipart(1)*a1mipart(4)+two*a1mipart(2)*a1mipart(3))
      interfpart3oneonenoeps(5)=four*yy*zz*(two*a1mipart(1)*a1mipart(5)+two*a1mipart(2)*a1mipart(4)+a1mipart(3)**2)

      do i=2,5
         interfpart3oneoneeps(i)=-interfpart3oneonenoeps(i-1)
      enddo

      interfpart3twotwo(1)=two*xx**2*yy/zz*(a2mipart(1)**2)
      interfpart3twotwo(2)=two*xx**2*yy/zz*(two*a2mipart(1)*a2mipart(2))
      interfpart3twotwo(3)=two*xx**2*yy/zz*(two*a2mipart(1)*a2mipart(3)+a2mipart(2)**2)
      interfpart3twotwo(4)=two*xx**2*yy/zz*(two*a2mipart(1)*a2mipart(4)+two*a2mipart(2)*a2mipart(3))
      interfpart3twotwo(5)=two*xx**2*yy/zz*(two*a2mipart(1)*a2mipart(5)+two*a2mipart(2)*a2mipart(4)+a2mipart(3)**2)

      interfpart3onetwo(1)=-two*xx*yy*(a1mipart(1)*a2mipart(1))
      interfpart3onetwo(2)=-two*xx*yy*(a1mipart(1)*a2mipart(2)+a1mipart(2)*a2mipart(1))
      interfpart3onetwo(3)=-two*xx*yy*(a1mipart(1)*a2mipart(3)+a1mipart(2)*a2mipart(2)+a1mipart(3)*a2mipart(1))
      interfpart3onetwo(4)=-two*xx*yy*(a1mipart(1)*a2mipart(4)+a1mipart(2)*a2mipart(3)+a1mipart(3)*a2mipart(2)+a1mipart(4)*a2mipart(1))
      interfpart3onetwo(5)=-two*xx*yy*(a1mipart(1)*a2mipart(5)+a1mipart(2)*a2mipart(4)+a1mipart(3)*a2mipart(3)+a1mipart(4)*a2mipart(2)+a1mipart(5)*a2mipart(1))

      do i=1,5
         interfpart3final(i)=interfpart3oneonenoeps(i)+interfpart3oneoneeps(i)+interfpart3twotwo(i)+two*interfpart3onetwo(i)
      enddo

!==== combining all terms

!==== multiply by LO color: Nc*cf
!==== divide by (2*pi)**3
!==== overall factor is as**3*yb**2 [mu**(2*eps) exp(eugam*eps)/(4*pi)**eps]

!==== NLO^2 hard (e.g. rescaled by LO msq and (as/2/pi)**2)

      do i=1,5
         interffinal(i)=Nc*cf/(128._dp*pi**2)*yy*zz/(xx**2+one)*(interfpart1final(i)+two*interfpart2final(i)+interfpart3final(i)-infraredpart1final(i)-infraredpart2final(i)-infraredpart3final(i))
      enddo

      nloxnloeps0=interffinal(5)

!      write(6,99) ' eps -4: ',interffinal(1)
!      write(6,99) ' eps -3: ',interffinal(2)
!      write(6,99) ' eps -2: ',interffinal(3)
!      write(6,99) ' eps -1: ',interffinal(4)
!      write(6,99) ' eps  0: ',interffinal(5)

!      write(6,*)
!      write(6,*) '**********************************************'
!      write(6,*) '* NLO^2  *'
!      write(6,*) '**********************************************'
!      write(6,*)
!      write(6,99) ' y = ',yy
!      write(6,99) ' z = ',zz
!      write(6,99) ' eps 0 part: ',nloxnloeps0
    
      return 
      
   99 format(a35,e20.12)
      
      end

