
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
      real(dp) xx,nloeps0(3)

      double precision Li2
      external Li2

c--- initialize to zero
      nloeps0(:)=0._dp
      nloeps0total=0._dp

c--- value for x
      xx = one-yy-zz

!==================================================================

!==== log(mh2omu2)^2 coefficient      
      nloeps0(3) = 0

!==== log(mh2omu2)^1 coefficient      
      nloeps0(2) = -47._dp/six

!==== log(mh2omu2)^0 coefficient      
      nloeps0(1) = -(96*xx**2 + 14*pi**2*xx**2 + 66*xx*yy + 14*pi**2*xx*yy + 
     -     18*yy**2 + 7*pi**2*yy**2 + 66*xx*zz + 14*pi**2*xx*zz + 
     -     36*yy*zz + 14*pi**2*yy*zz + 18*zz**2 + 7*pi**2*zz**2 + 
     -     24*xx**2*(Li2(xx) + log(1 - xx)*log(xx)) + 
     -     24*xx*yy*(Li2(xx) + log(1 - xx)*log(xx)) + 
     -     12*yy**2*(Li2(xx) + log(1 - xx)*log(xx)) + 
     -     24*xx*zz*(Li2(xx) + log(1 - xx)*log(xx)) + 
     -     24*yy*zz*(Li2(xx) + log(1 - xx)*log(xx)) + 
     -     12*zz**2*(Li2(xx) + log(1 - xx)*log(xx)) - 
     -     96*xx**2*(Li2(yy) + log(1 - yy)*log(yy)) - 
     -     96*xx*yy*(Li2(yy) + log(1 - yy)*log(yy)) - 
     -     48*yy**2*(Li2(yy) + log(1 - yy)*log(yy)) - 
     -     96*xx*zz*(Li2(yy) + log(1 - yy)*log(yy)) - 
     -     96*yy*zz*(Li2(yy) + log(1 - yy)*log(yy)) - 
     -     48*zz**2*(Li2(yy) + log(1 - yy)*log(yy)) + 
     -     150*xx**2*log(zz) + 150*xx*yy*log(zz) + 75*yy**2*log(zz) + 
     -     150*xx*zz*log(zz) + 150*yy*zz*log(zz) + 75*zz**2*log(zz) - 
     -     3*(2*xx**2 + 2*xx*(yy + zz) + (yy + zz)**2)*log(xx)*
     -      (3 + 2*log(yy) + 2*log(zz)) + 
     -     3*(2*xx**2 + 2*xx*(yy + zz) + (yy + zz)**2)*log(yy)*
     -      (25 + 18*log(zz)) - 
     -     96*xx**2*(Li2(zz) + log(1 - zz)*log(zz)) - 
     -     96*xx*yy*(Li2(zz) + log(1 - zz)*log(zz)) - 
     -     48*yy**2*(Li2(zz) + log(1 - zz)*log(zz)) - 
     -     96*xx*zz*(Li2(zz) + log(1 - zz)*log(zz)) - 
     -     96*yy*zz*(Li2(zz) + log(1 - zz)*log(zz)) - 
     -     48*zz**2*(Li2(zz) + log(1 - zz)*log(zz)))/
     -  (18.*(xx**2+1))

      nloeps0total=logmh2omu2**2*nloeps0(3)+logmh2omu2*nloeps0(2)+nloeps0(1)

!===================================================================

!      write(6,*)
!      write(6,*) '**********************************************'
!      write(6,*) '* NLOh  *'
!      write(6,*) '**********************************************'
!      write(6,*)
!      write(6,99) ' y = ',yy
!      write(6,99) ' z = ',zz
!      write(6,99) ' eps 0 part, log(mh2omu2)^2 coeff: ',nloeps0(3)
!      write(6,99) ' eps 0 part, log(mh2omu2)^1 coeff: ',nloeps0(2)
!      write(6,99) ' eps 0 part, log(mh2omu2)^0 coeff: ',nloeps0(1)
!      write(6,*)
!      write(6,99) ' eps 0 part, total: ',nloeps0total

      return 
      
   99 format(a35,e15.8)
      
      end

