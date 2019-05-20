
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
      nloeps0(2) = -2*cf - 1/(2.*Nc)

!==== log(mh2omu2)^1 coefficient      
      nloeps0(1) = -((yy**2 + 4*cf*Nc*yy**2 + 2*yy*zz + 8*cf*Nc*yy*zz + zz**2 + 
     -      4*cf*Nc*zz**2 - (2*xx**2 + 2*xx*(yy + zz) + (yy + zz)**2)*
     -       log(xx) + (1 + 2*cf*Nc)*
     -       (2*xx**2 + 2*xx*(yy + zz) + (yy + zz)**2)*log(yy) + 
     -      2*xx**2*log(zz) + 4*cf*Nc*xx**2*log(zz) + 2*xx*yy*log(zz) + 
     -      4*cf*Nc*xx*yy*log(zz) + yy**2*log(zz) + 
     -      2*cf*Nc*yy**2*log(zz) + 2*xx*zz*log(zz) + 
     -      4*cf*Nc*xx*zz*log(zz) + 2*yy*zz*log(zz) + 
     -      4*cf*Nc*yy*zz*log(zz) + zz**2*log(zz) + 
     -      2*cf*Nc*zz**2*log(zz))/
     -    (Nc*(xx**2+1)))

!==== log(mh2omu2)^0 coefficient      
      nloeps0(0) = (-48*cf*Nc*xx**2 + 18*pi**2*xx**2 + 48*cf*Nc*pi**2*xx**2 + 
     -    12*xx*yy - 36*cf*Nc*xx*yy + 18*pi**2*xx*yy + 
     -    48*cf*Nc*pi**2*xx*yy + 12*yy**2 + 22*ca*Nc*yy**2 + 
     -    24*cf*Nc*yy**2 + 9*pi**2*yy**2 + 24*cf*Nc*pi**2*yy**2 - 
     -    8*Nc*nf*tr*yy**2 + 12*xx*zz - 36*cf*Nc*xx*zz + 
     -    18*pi**2*xx*zz + 48*cf*Nc*pi**2*xx*zz + 24*yy*zz + 
     -    44*ca*Nc*yy*zz + 48*cf*Nc*yy*zz + 18*pi**2*yy*zz + 
     -    48*cf*Nc*pi**2*yy*zz - 16*Nc*nf*tr*yy*zz + 12*zz**2 + 
     -    22*ca*Nc*zz**2 + 24*cf*Nc*zz**2 + 9*pi**2*zz**2 + 
     -    24*cf*Nc*pi**2*zz**2 - 8*Nc*nf*tr*zz**2 + 
     -    12*xx**2*log(xx)**2 + 12*xx*yy*log(xx)**2 + 
     -    6*yy**2*log(xx)**2 + 12*xx*zz*log(xx)**2 + 
     -    12*yy*zz*log(xx)**2 + 6*zz**2*log(xx)**2 - 
     -    48*xx**2*(Li2(xx) + log(1 - xx)*log(xx)) - 
     -    48*xx*yy*(Li2(xx) + log(1 - xx)*log(xx)) - 
     -    24*yy**2*(Li2(xx) + log(1 - xx)*log(xx)) - 
     -    48*xx*zz*(Li2(xx) + log(1 - xx)*log(xx)) - 
     -    48*yy*zz*(Li2(xx) + log(1 - xx)*log(xx)) - 
     -    24*zz**2*(Li2(xx) + log(1 - xx)*log(xx)) - 
     -    12*xx**2*log(yy)**2 - 24*cf*Nc*xx**2*log(yy)**2 - 
     -    12*xx*yy*log(yy)**2 - 24*cf*Nc*xx*yy*log(yy)**2 - 
     -    6*yy**2*log(yy)**2 - 12*cf*Nc*yy**2*log(yy)**2 - 
     -    12*xx*zz*log(yy)**2 - 24*cf*Nc*xx*zz*log(yy)**2 - 
     -    12*yy*zz*log(yy)**2 - 24*cf*Nc*yy*zz*log(yy)**2 - 
     -    6*zz**2*log(yy)**2 - 12*cf*Nc*zz**2*log(yy)**2 + 
     -    48*cf*Nc*xx**2*(Li2(yy) + log(1 - yy)*log(yy)) + 
     -    48*cf*Nc*xx*yy*(Li2(yy) + log(1 - yy)*log(yy)) + 
     -    24*cf*Nc*yy**2*(Li2(yy) + log(1 - yy)*log(yy)) + 
     -    48*cf*Nc*xx*zz*(Li2(yy) + log(1 - yy)*log(yy)) + 
     -    48*cf*Nc*yy*zz*(Li2(yy) + log(1 - yy)*log(yy)) + 
     -    24*cf*Nc*zz**2*(Li2(yy) + log(1 - yy)*log(yy)) - 
     -    12*yy**2*log(zz) - 24*cf*Nc*yy**2*log(zz) - 
     -    24*yy*zz*log(zz) - 48*cf*Nc*yy*zz*log(zz) - 
     -    12*zz**2*log(zz) - 24*cf*Nc*zz**2*log(zz) - 
     -    12*xx**2*log(zz)**2 - 24*cf*Nc*xx**2*log(zz)**2 - 
     -    12*xx*yy*log(zz)**2 - 24*cf*Nc*xx*yy*log(zz)**2 - 
     -    6*yy**2*log(zz)**2 - 12*cf*Nc*yy**2*log(zz)**2 - 
     -    12*xx*zz*log(zz)**2 - 24*cf*Nc*xx*zz*log(zz)**2 - 
     -    12*yy*zz*log(zz)**2 - 24*cf*Nc*yy*zz*log(zz)**2 - 
     -    6*zz**2*log(zz)**2 - 12*cf*Nc*zz**2*log(zz)**2 - 
     -    12*(1 + 2*cf*Nc)*log(yy)*
     -     ((yy + zz)**2 + (2*xx**2 + 2*xx*(yy + zz) + (yy + zz)**2)*
     -        log(zz)) + 12*log(xx)*
     -     ((yy + zz)**2 + (2*xx**2 + 2*xx*(yy + zz) + (yy + zz)**2)*
     -        log(yy) + (2*xx**2 + 2*xx*(yy + zz) + (yy + zz)**2)*
     -        log(zz)) + 48*cf*Nc*xx**2*
     -     (Li2(zz) + log(1 - zz)*log(zz)) + 
     -    48*cf*Nc*xx*yy*(Li2(zz) + log(1 - zz)*log(zz)) + 
     -    24*cf*Nc*yy**2*(Li2(zz) + log(1 - zz)*log(zz)) + 
     -    48*cf*Nc*xx*zz*(Li2(zz) + log(1 - zz)*log(zz)) + 
     -    48*cf*Nc*yy*zz*(Li2(zz) + log(1 - zz)*log(zz)) + 
     -    24*cf*Nc*zz**2*(Li2(zz) + log(1 - zz)*log(zz)))/
     -  (12.*Nc*(xx**2+1))

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

