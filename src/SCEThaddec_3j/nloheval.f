
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
      nloeps0(1) = (log(xx) - (1 + 2*cf*Nc)*(log(yy) + log(zz)))/Nc

!==== log(mh2omu2)^0 coefficient      
      nloeps0(0) = (-16*cf*Nc*xx**2 + 6*pi**2*xx**2 + 16*cf*Nc*pi**2*xx**2 + 
     -    4*xx*yy - 12*cf*Nc*xx*yy + 6*pi**2*xx*yy + 
     -    16*cf*Nc*pi**2*xx*yy + 4*yy**2 - 4*cf*Nc*yy**2 + 
     -    3*pi**2*yy**2 + 8*cf*Nc*pi**2*yy**2 + 4*xx*zz - 
     -    12*cf*Nc*xx*zz + 6*pi**2*xx*zz + 16*cf*Nc*pi**2*xx*zz + 
     -    8*yy*zz - 8*cf*Nc*yy*zz + 6*pi**2*yy*zz + 
     -    16*cf*Nc*pi**2*yy*zz + 4*zz**2 - 4*cf*Nc*zz**2 + 
     -    3*pi**2*zz**2 + 8*cf*Nc*pi**2*zz**2 - 
     -    8*(2*xx**2 + 2*xx*(yy + zz) + (yy + zz)**2)*Li2(xx) + 
     -    8*cf*Nc*(2*xx**2 + 2*xx*(yy + zz) + (yy + zz)**2)*Li2(yy) + 
     -    16*cf*Nc*xx**2*Li2(zz) + 16*cf*Nc*xx*yy*Li2(zz) + 
     -    8*cf*Nc*yy**2*Li2(zz) + 16*cf*Nc*xx*zz*Li2(zz) + 
     -    16*cf*Nc*yy*zz*Li2(zz) + 8*cf*Nc*zz**2*Li2(zz) - 
     -    16*xx**2*log(1 - xx)*log(xx) - 16*xx*yy*log(1 - xx)*log(xx) - 
     -    8*yy**2*log(1 - xx)*log(xx) - 16*xx*zz*log(1 - xx)*log(xx) - 
     -    16*yy*zz*log(1 - xx)*log(xx) - 8*zz**2*log(1 - xx)*log(xx) + 
     -    4*xx**2*log(xx)**2 + 4*xx*yy*log(xx)**2 + 
     -    2*yy**2*log(xx)**2 + 4*xx*zz*log(xx)**2 + 
     -    4*yy*zz*log(xx)**2 + 2*zz**2*log(xx)**2 + 
     -    8*xx**2*log(xx)*log(yy) + 8*xx*yy*log(xx)*log(yy) + 
     -    4*yy**2*log(xx)*log(yy) + 8*xx*zz*log(xx)*log(yy) + 
     -    8*yy*zz*log(xx)*log(yy) + 4*zz**2*log(xx)*log(yy) + 
     -    16*cf*Nc*xx**2*log(1 - yy)*log(yy) + 
     -    16*cf*Nc*xx*yy*log(1 - yy)*log(yy) + 
     -    8*cf*Nc*yy**2*log(1 - yy)*log(yy) + 
     -    16*cf*Nc*xx*zz*log(1 - yy)*log(yy) + 
     -    16*cf*Nc*yy*zz*log(1 - yy)*log(yy) + 
     -    8*cf*Nc*zz**2*log(1 - yy)*log(yy) - 4*xx**2*log(yy)**2 - 
     -    8*cf*Nc*xx**2*log(yy)**2 - 4*xx*yy*log(yy)**2 - 
     -    8*cf*Nc*xx*yy*log(yy)**2 - 2*yy**2*log(yy)**2 - 
     -    4*cf*Nc*yy**2*log(yy)**2 - 4*xx*zz*log(yy)**2 - 
     -    8*cf*Nc*xx*zz*log(yy)**2 - 4*yy*zz*log(yy)**2 - 
     -    8*cf*Nc*yy*zz*log(yy)**2 - 2*zz**2*log(yy)**2 - 
     -    4*cf*Nc*zz**2*log(yy)**2 + 8*xx**2*log(xx)*log(zz) + 
     -    8*xx*yy*log(xx)*log(zz) + 4*yy**2*log(xx)*log(zz) + 
     -    8*xx*zz*log(xx)*log(zz) + 8*yy*zz*log(xx)*log(zz) + 
     -    4*zz**2*log(xx)*log(zz) - 8*xx**2*log(yy)*log(zz) - 
     -    16*cf*Nc*xx**2*log(yy)*log(zz) - 8*xx*yy*log(yy)*log(zz) - 
     -    16*cf*Nc*xx*yy*log(yy)*log(zz) - 4*yy**2*log(yy)*log(zz) - 
     -    8*cf*Nc*yy**2*log(yy)*log(zz) - 8*xx*zz*log(yy)*log(zz) - 
     -    16*cf*Nc*xx*zz*log(yy)*log(zz) - 8*yy*zz*log(yy)*log(zz) - 
     -    16*cf*Nc*yy*zz*log(yy)*log(zz) - 4*zz**2*log(yy)*log(zz) - 
     -    8*cf*Nc*zz**2*log(yy)*log(zz) + 
     -    16*cf*Nc*xx**2*log(1 - zz)*log(zz) + 
     -    16*cf*Nc*xx*yy*log(1 - zz)*log(zz) + 
     -    8*cf*Nc*yy**2*log(1 - zz)*log(zz) + 
     -    16*cf*Nc*xx*zz*log(1 - zz)*log(zz) + 
     -    16*cf*Nc*yy*zz*log(1 - zz)*log(zz) + 
     -    8*cf*Nc*zz**2*log(1 - zz)*log(zz) - 4*xx**2*log(zz)**2 - 
     -    8*cf*Nc*xx**2*log(zz)**2 - 4*xx*yy*log(zz)**2 - 
     -    8*cf*Nc*xx*yy*log(zz)**2 - 2*yy**2*log(zz)**2 - 
     -    4*cf*Nc*yy**2*log(zz)**2 - 4*xx*zz*log(zz)**2 - 
     -    8*cf*Nc*xx*zz*log(zz)**2 - 4*yy*zz*log(zz)**2 - 
     -    8*cf*Nc*yy*zz*log(zz)**2 - 2*zz**2*log(zz)**2 - 
     -    4*cf*Nc*zz**2*log(zz)**2)/
     -  (4.*Nc*(xx**2+1))

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

