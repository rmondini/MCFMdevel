
      subroutine acoeffnloinfraredeval(yy,zz,logmh2omu2,a1infrared,a2infrared) 
      implicit none 
      include 'types.f'
      include 'constants.f' 
      include 'ewcouple.f' 
      include 'masses.f' 
      include 'nf.f'
      include 'scet_const.f'
      real(dp) xx
      real(dp),intent(in)::yy,zz,logmh2omu2
      real(dp),intent(out)::a1infrared(5),a2infrared(5)

      double precision Li2,Li3,Li4
      external Li2,Li3,Li4

c--- initialize to zero
      a1infrared(:)=0._dp
      a2infrared(:)=0._dp

      xx=one-yy-zz

!==================================================================

!==== eps^(-2)
      a1infrared(1) = -(((ca + 2*cf)*pi*rt2*(yy + zz))/(yy*zz)) 

!==== eps^(-1)
      a1infrared(2) = (pi*(yy + zz)*(-11*ca - 18*cf + 6*(ca + 2*cf)*logmh2omu2 + 
     -      2*nf - 6*(ca - 2*cf)*log(xx) + 6*ca*log(yy) + 6*ca*log(zz)))
     -   /(3.*rt2*yy*zz) 

!==== eps^( 0)
      a1infrared(3) = (pi*(yy + zz)*(-6*(ca + 2*cf)*logmh2omu2**2 + ca*pi**2 + 
     -      2*cf*pi**2 - 18*(ca - 2*cf)*log(xx) + 
     -      6*(ca - 2*cf)*log(xx)**2 + 20*ca*log(yy) - 2*nf*log(yy) - 
     -      6*ca*log(yy)**2 + 20*ca*log(zz) - 2*nf*log(zz) - 
     -      6*ca*log(zz)**2 + 
     -      2*logmh2omu2*(11*ca + 18*cf - 2*nf + 
     -         6*(ca - 2*cf)*log(xx) - 6*ca*log(yy) - 6*ca*log(zz))))/
     -  (6.*rt2*yy*zz) 

!==== eps^( 1)
      a1infrared(4) = -(pi*(yy + zz)*(-12*(ca + 2*cf)*logmh2omu2**3 - 11*ca*pi**2 - 
     -       18*cf*pi**2 + 2*nf*pi**2 - 24*ca*zeta3 - 48*cf*zeta3 - 
     -       6*(ca - 2*cf)*pi**2*log(xx) - 54*(ca - 2*cf)*log(xx)**2 + 
     -       12*(ca - 2*cf)*log(xx)**3 + 6*ca*pi**2*log(yy) + 
     -       60*ca*log(yy)**2 - 6*nf*log(yy)**2 - 12*ca*log(yy)**3 + 
     -       6*ca*pi**2*log(zz) + 60*ca*log(zz)**2 - 6*nf*log(zz)**2 - 
     -       12*ca*log(zz)**3 + 
     -       6*logmh2omu2**2*
     -        (11*ca + 18*cf - 2*nf + 6*(ca - 2*cf)*log(xx) - 
     -          6*ca*log(yy) - 6*ca*log(zz)) + 
     -       6*logmh2omu2*(ca*pi**2 + 2*cf*pi**2 - 
     -          18*(ca - 2*cf)*log(xx) + 6*(ca - 2*cf)*log(xx)**2 + 
     -          20*ca*log(yy) - 2*nf*log(yy) - 6*ca*log(yy)**2 + 
     -          20*ca*log(zz) - 2*nf*log(zz) - 6*ca*log(zz)**2)))/
     -  (36.*rt2*yy*zz) 

!==== eps^( 2)
      a1infrared(5) = -(pi*(yy + zz)*(60*(ca + 2*cf)*logmh2omu2**4 + ca*pi**4 + 
     -       2*cf*pi**4 - 880*ca*zeta3 - 1440*cf*zeta3 + 160*nf*zeta3 - 
     -       60*(ca - 2*cf)*(3*pi**2 + 8*zeta3)*log(xx) + 
     -       60*(ca - 2*cf)*pi**2*log(xx)**2 + 
     -       360*(ca - 2*cf)*log(xx)**3 - 60*(ca - 2*cf)*log(xx)**4 + 
     -       200*ca*pi**2*log(yy) - 20*nf*pi**2*log(yy) + 
     -       480*ca*zeta3*log(yy) - 60*ca*pi**2*log(yy)**2 - 
     -       400*ca*log(yy)**3 + 40*nf*log(yy)**3 + 60*ca*log(yy)**4 + 
     -       200*ca*pi**2*log(zz) - 20*nf*pi**2*log(zz) + 
     -       480*ca*zeta3*log(zz) - 60*ca*pi**2*log(zz)**2 - 
     -       400*ca*log(zz)**3 + 40*nf*log(zz)**3 + 60*ca*log(zz)**4 - 
     -       40*logmh2omu2**3*
     -        (11*ca + 18*cf - 2*nf + 6*(ca - 2*cf)*log(xx) - 
     -          6*ca*log(yy) - 6*ca*log(zz)) - 
     -       60*logmh2omu2**2*
     -        (ca*pi**2 + 2*cf*pi**2 - 18*(ca - 2*cf)*log(xx) + 
     -          6*(ca - 2*cf)*log(xx)**2 + (20*ca - 2*nf)*log(yy) - 
     -          6*ca*log(yy)**2 + 20*ca*log(zz) - 2*nf*log(zz) - 
     -          6*ca*log(zz)**2) + 
     -       20*logmh2omu2*(11*ca*pi**2 + 18*cf*pi**2 - 2*nf*pi**2 + 
     -          24*ca*zeta3 + 48*cf*zeta3 + 
     -          6*(ca - 2*cf)*pi**2*log(xx) + 
     -          54*(ca - 2*cf)*log(xx)**2 - 12*(ca - 2*cf)*log(xx)**3 - 
     -          6*ca*pi**2*log(yy) + (-60*ca + 6*nf)*log(yy)**2 + 
     -          12*ca*log(yy)**3 - 6*ca*pi**2*log(zz) - 
     -          60*ca*log(zz)**2 + 6*nf*log(zz)**2 + 12*ca*log(zz)**3)))
     -   /(720.*rt2*yy*zz) 

!==== eps^(-2)
      a2infrared(1) = (2*(ca + 2*cf)*pi*rt2)/yy

!==== eps^(-1)
      a2infrared(2) = -(pi*rt2*(-11*ca - 18*cf + 6*(ca + 2*cf)*logmh2omu2 + 2*nf - 
     -       6*(ca - 2*cf)*log(xx) + 6*ca*log(yy) + 6*ca*log(zz)))/
     -  (3.*yy)

!==== eps^( 0)
      a2infrared(3) = -(pi*(-6*(ca + 2*cf)*logmh2omu2**2 + ca*pi**2 + 2*cf*pi**2 - 
     -       18*(ca - 2*cf)*log(xx) + 6*(ca - 2*cf)*log(xx)**2 + 
     -       20*ca*log(yy) - 2*nf*log(yy) - 6*ca*log(yy)**2 + 
     -       20*ca*log(zz) - 2*nf*log(zz) - 6*ca*log(zz)**2 + 
     -       2*logmh2omu2*(11*ca + 18*cf - 2*nf + 
     -          6*(ca - 2*cf)*log(xx) - 6*ca*log(yy) - 6*ca*log(zz))))/
     -  (3.*rt2*yy)

!==== eps^( 1)
      a2infrared(4) = (pi*(-12*(ca + 2*cf)*logmh2omu2**3 - 11*ca*pi**2 - 
     -      18*cf*pi**2 + 2*nf*pi**2 - 24*ca*zeta3 - 48*cf*zeta3 - 
     -      6*(ca - 2*cf)*pi**2*log(xx) - 54*(ca - 2*cf)*log(xx)**2 + 
     -      12*(ca - 2*cf)*log(xx)**3 + 6*ca*pi**2*log(yy) + 
     -      60*ca*log(yy)**2 - 6*nf*log(yy)**2 - 12*ca*log(yy)**3 + 
     -      6*ca*pi**2*log(zz) + 60*ca*log(zz)**2 - 6*nf*log(zz)**2 - 
     -      12*ca*log(zz)**3 + 
     -      6*logmh2omu2**2*(11*ca + 18*cf - 2*nf + 
     -         6*(ca - 2*cf)*log(xx) - 6*ca*log(yy) - 6*ca*log(zz)) + 
     -      6*logmh2omu2*(ca*pi**2 + 2*cf*pi**2 - 
     -         18*(ca - 2*cf)*log(xx) + 6*(ca - 2*cf)*log(xx)**2 + 
     -         20*ca*log(yy) - 2*nf*log(yy) - 6*ca*log(yy)**2 + 
     -         20*ca*log(zz) - 2*nf*log(zz) - 6*ca*log(zz)**2)))/
     -  (18.*rt2*yy)

!==== eps^( 2)
      a2infrared(5) = (pi*(60*(ca + 2*cf)*logmh2omu2**4 + ca*pi**4 + 2*cf*pi**4 - 
     -      880*ca*zeta3 - 1440*cf*zeta3 + 160*nf*zeta3 - 
     -      60*(ca - 2*cf)*(3*pi**2 + 8*zeta3)*log(xx) + 
     -      60*(ca - 2*cf)*pi**2*log(xx)**2 + 
     -      360*(ca - 2*cf)*log(xx)**3 - 60*(ca - 2*cf)*log(xx)**4 + 
     -      200*ca*pi**2*log(yy) - 20*nf*pi**2*log(yy) + 
     -      480*ca*zeta3*log(yy) - 60*ca*pi**2*log(yy)**2 - 
     -      400*ca*log(yy)**3 + 40*nf*log(yy)**3 + 60*ca*log(yy)**4 + 
     -      200*ca*pi**2*log(zz) - 20*nf*pi**2*log(zz) + 
     -      480*ca*zeta3*log(zz) - 60*ca*pi**2*log(zz)**2 - 
     -      400*ca*log(zz)**3 + 40*nf*log(zz)**3 + 60*ca*log(zz)**4 - 
     -      40*logmh2omu2**3*
     -       (11*ca + 18*cf - 2*nf + 6*(ca - 2*cf)*log(xx) - 
     -         6*ca*log(yy) - 6*ca*log(zz)) - 
     -      60*logmh2omu2**2*
     -       (ca*pi**2 + 2*cf*pi**2 - 18*(ca - 2*cf)*log(xx) + 
     -         6*(ca - 2*cf)*log(xx)**2 + (20*ca - 2*nf)*log(yy) - 
     -         6*ca*log(yy)**2 + 20*ca*log(zz) - 2*nf*log(zz) - 
     -         6*ca*log(zz)**2) + 
     -      20*logmh2omu2*(11*ca*pi**2 + 18*cf*pi**2 - 2*nf*pi**2 + 
     -         24*ca*zeta3 + 48*cf*zeta3 + 
     -         6*(ca - 2*cf)*pi**2*log(xx) + 
     -         54*(ca - 2*cf)*log(xx)**2 - 12*(ca - 2*cf)*log(xx)**3 - 
     -         6*ca*pi**2*log(yy) + (-60*ca + 6*nf)*log(yy)**2 + 
     -         12*ca*log(yy)**3 - 6*ca*pi**2*log(zz) - 
     -         60*ca*log(zz)**2 + 6*nf*log(zz)**2 + 12*ca*log(zz)**3)))/
     -  (360.*rt2*yy) 

!===================================================================

      return 
      
   99 format(a35,e15.8)
      
      end

