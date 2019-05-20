
      subroutine nloheval(yy,zz,logmh2omu2,nloeps0total) 
      implicit none 
      include 'types.f'
      include 'constants.f' 
      include 'ewcouple.f' 
      include 'masses.f'
      include 'nf.f'
      include 'scet_const.f'
      real(dp),intent(in)::yy,zz,logmh2omu2
      real(dp),intent(out)::nloeps0total
      real(dp) xx,nloeps0(0:2)

      double precision Li2
      external Li2

c--- initialize to zero
      nloeps0(:)=0._dp
      nloeps0total=0._dp

c--- value for x
      xx = one-yy-zz

!==================================================================

!==== log(mh2omu2)^2 coefficient      
      nloeps0(2) = -17._dp/six

!==== log(mh2omu2)^1 coefficient      
      nloeps0(1) = -(17*yy**2 + 34*yy*zz + 17*zz**2 - 
     -     (2*xx**2 + 2*xx*(yy + zz) + (yy + zz)**2)*log(xx) + 
     -     9*(2*xx**2 + 2*xx*(yy + zz) + (yy + zz)**2)*log(yy) + 
     -     18*xx**2*log(zz) + 18*xx*yy*log(zz) + 9*yy**2*log(zz) + 
     -     18*xx*zz*log(zz) + 18*yy*zz*log(zz) + 9*zz**2*log(zz))/
     -  (3.*(xx**2+1))

!==== log(mh2omu2)^0 coefficient      
      nloeps0(0) = (-64*xx**2 + 70*pi**2*xx**2 - 44*xx*yy + 70*pi**2*xx*yy + 
     -    82*yy**2 + 35*pi**2*yy**2 - 44*xx*zz + 70*pi**2*xx*zz + 
     -    164*yy*zz + 70*pi**2*yy*zz + 82*zz**2 + 35*pi**2*zz**2 + 
     -    4*xx**2*log(xx)**2 + 4*xx*yy*log(xx)**2 + 
     -    2*yy**2*log(xx)**2 + 4*xx*zz*log(xx)**2 + 
     -    4*yy*zz*log(xx)**2 + 2*zz**2*log(xx)**2 - 
     -    16*xx**2*(Li2(xx) + log(1 - xx)*log(xx)) - 
     -    16*xx*yy*(Li2(xx) + log(1 - xx)*log(xx)) - 
     -    8*yy**2*(Li2(xx) + log(1 - xx)*log(xx)) - 
     -    16*xx*zz*(Li2(xx) + log(1 - xx)*log(xx)) - 
     -    16*yy*zz*(Li2(xx) + log(1 - xx)*log(xx)) - 
     -    8*zz**2*(Li2(xx) + log(1 - xx)*log(xx)) - 
     -    36*xx**2*log(yy)**2 - 36*xx*yy*log(yy)**2 - 
     -    18*yy**2*log(yy)**2 - 36*xx*zz*log(yy)**2 - 
     -    36*yy*zz*log(yy)**2 - 18*zz**2*log(yy)**2 + 
     -    64*xx**2*(Li2(yy) + log(1 - yy)*log(yy)) + 
     -    64*xx*yy*(Li2(yy) + log(1 - yy)*log(yy)) + 
     -    32*yy**2*(Li2(yy) + log(1 - yy)*log(yy)) + 
     -    64*xx*zz*(Li2(yy) + log(1 - yy)*log(yy)) + 
     -    64*yy*zz*(Li2(yy) + log(1 - yy)*log(yy)) + 
     -    32*zz**2*(Li2(yy) + log(1 - yy)*log(yy)) - 36*yy**2*log(zz) - 
     -    72*yy*zz*log(zz) - 36*zz**2*log(zz) - 36*xx**2*log(zz)**2 - 
     -    36*xx*yy*log(zz)**2 - 18*yy**2*log(zz)**2 - 
     -    36*xx*zz*log(zz)**2 - 36*yy*zz*log(zz)**2 - 
     -    18*zz**2*log(zz)**2 - 
     -    36*log(yy)*((yy + zz)**2 + 
     -       (2*xx**2 + 2*xx*(yy + zz) + (yy + zz)**2)*log(zz)) + 
     -    4*log(xx)*((yy + zz)**2 + 
     -       (2*xx**2 + 2*xx*(yy + zz) + (yy + zz)**2)*log(yy) + 
     -       (2*xx**2 + 2*xx*(yy + zz) + (yy + zz)**2)*log(zz)) + 
     -    64*xx**2*(Li2(zz) + log(1 - zz)*log(zz)) + 
     -    64*xx*yy*(Li2(zz) + log(1 - zz)*log(zz)) + 
     -    32*yy**2*(Li2(zz) + log(1 - zz)*log(zz)) + 
     -    64*xx*zz*(Li2(zz) + log(1 - zz)*log(zz)) + 
     -    64*yy*zz*(Li2(zz) + log(1 - zz)*log(zz)) + 
     -    32*zz**2*(Li2(zz) + log(1 - zz)*log(zz)))/
     -  (12.*(xx**2+1))

      nloeps0total=logmh2omu2**2*nloeps0(2)+logmh2omu2*nloeps0(1)+nloeps0(0)

!===================================================================

!      write(6,*)
!      write(6,*) '**********************************************'
!      write(6,*) '* NLOh  *'
!      write(6,*) '**********************************************'
!      write(6,*)
!      write(6,99) ' y = ',yy
!      write(6,99) ' z = ',zz
!      write(6,99) ' eps 0 part, log(mh2omu2)^2 coeff: ',nloeps0(2)
!      write(6,99) ' eps 0 part, log(mh2omu2)^1 coeff: ',nloeps0(1)
!      write(6,99) ' eps 0 part, log(mh2omu2)^0 coeff: ',nloeps0(0)
!      write(6,*)
!      write(6,99) ' eps 0 part, total: ',nloeps0total

      return 
      
   99 format(a35,e15.8)
      
      end

