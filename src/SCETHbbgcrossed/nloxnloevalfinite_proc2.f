
      subroutine nloxnloevalfinite_proc2(uu,vv,musq,nloxnloeps0) 
      implicit none 
      include 'types.f'
      include 'constants.f'
      include 'nf.f'  
      include 'ewcouple.f' 
      include 'masses.f'
      include 'scet_const.f'
      integer i
      real(dp) xx,yy,zz
      real(dp) a1born,a2born

      complex(dp) a1finitelogcoeff(3),a2finitelogcoeff(3)
      complex(dp) a1finite,a2finite
      complex(dp) lnrat
      real(dp) H1(0:1),H2(0:1,0:1),H3(0:1,0:1,0:1),H4(0:1,0:1,0:1,0:1)
      real(dp) G1(0:3),G2(0:3,0:3),G3(0:3,0:3,0:3),G4(0:3,0:3,0:3,0:3)

      real(dp) musq
      real(dp),intent(in)::uu,vv
      real(dp),intent(out)::nloxnloeps0

c---- initialize to zero
      nloxnloeps0=0._dp

c---- value for x,y,z
      yy = -uu/vv
      zz = one/vv
      xx = one-yy-zz

c---- form factors at LO with overall factor sqrt(as/(2*pi))
      a1born=two*rt2*pi*(1/yy+1/zz)
      a2born=two*rt2*pi*(-2/yy)

!==================================================================

      call tdhpl(uu,vv,4,G1,G2,G3,G4,H1,H2,H3,H4)

!==== MSbar-renormalized coefficients for hard function

      a1finitelogcoeff(1) = -cf - 1/(4.*Nc)
      a1finitelogcoeff(2) = -((0,2)*cf*Nc*pi + (1 + 2*cf*Nc)*G1(0) - G1(2) - 
     -     (1 + 4*cf*Nc)*H1(0) + H1(1))/(2.*Nc)
      a1finitelogcoeff(3) = (12 - 12*cf*Nc - 5*pi**2 + 20*cf*Nc*pi**2 - (0,24)*pi*G1(1) + 
     -    (0,24)*pi*G1(2) + (0,24)*cf*Nc*pi*G1(3) - 12*G2(0,0) - 
     -    24*cf*Nc*G2(0,0) + 12*G2(0,2) - 24*G2(1,2) + 12*G2(2,0) + 
     -    12*G2(2,2) + 24*cf*Nc*G2(3,0) + 
     -    12*((0,4)*cf*Nc*pi + 2*G1(1) - 3*G1(2) - 2*cf*Nc*G1(3))*
     -     H1(0) + 12*G1(0)*((0,-2)*cf*Nc*pi + (1 + 4*cf*Nc)*H1(0) - 
     -       H1(1)) - (0,24)*pi*H1(1) + 24*G1(1)*H1(1) - 
     -    12*G1(2)*H1(1) - 12*H2(0,0) - 72*cf*Nc*H2(0,0) + 12*H2(0,1) + 
     -    36*H2(1,0) + 24*cf*Nc*H2(1,0) + 12*H2(1,1))/(24.*Nc) 

      a2finitelogcoeff(1) = -cf - 1/(4.*Nc) 
      a2finitelogcoeff(2) = -((0,2)*cf*Nc*pi + (1 + 2*cf*Nc)*G1(0) - G1(2) - 
     -     (1 + 4*cf*Nc)*H1(0) + H1(1))/(2.*Nc)  
      a2finitelogcoeff(3) = (-24*cf*Nc - 5*pi**2 + 20*cf*Nc*pi**2 - (0,24)*pi*G1(1) + 
     -    (0,24)*pi*G1(2) + (0,24)*cf*Nc*pi*G1(3) - 12*G2(0,0) - 
     -    24*cf*Nc*G2(0,0) + 12*G2(0,2) - 24*G2(1,2) + 12*G2(2,0) + 
     -    12*G2(2,2) + 24*cf*Nc*G2(3,0) + 
     -    12*((0,4)*cf*Nc*pi + 2*G1(1) - 3*G1(2) - 2*cf*Nc*G1(3))*
     -     H1(0) + 12*G1(0)*((0,-2)*cf*Nc*pi + (1 + 4*cf*Nc)*H1(0) - 
     -       H1(1)) - (0,24)*pi*H1(1) + 24*G1(1)*H1(1) - 
     -    12*G1(2)*H1(1) - 12*H2(0,0) - 72*cf*Nc*H2(0,0) + 12*H2(0,1) + 
     -    36*H2(1,0) + 24*cf*Nc*H2(1,0) + 12*H2(1,1))/(24.*Nc)  

      a1finite = a1finitelogcoeff(1)*lnrat(-hmass**2,musq)**2
     & + a1finitelogcoeff(2)*lnrat(-hmass**2,musq)
     & + a1finitelogcoeff(3)
      a1finite = a1born*a1finite

      a2finite = a2finitelogcoeff(1)*lnrat(-hmass**2,musq)**2
     & + a2finitelogcoeff(2)*lnrat(-hmass**2,musq)
     & + a2finitelogcoeff(3)
      a2finite = a2born*a2finite

!      write(*,*) 'a1finite nloxnlo proc2 = ',a1finite
!      write(*,*) 'a2finite nloxnlo proc2 = ',a2finite

!==== construction of all required pieces of the interference

      nloxnloeps0=
     & + four*yy*zz*(conjg(a1finite)*a1finite)
     & + two*xx**2*yy/zz*(conjg(a2finite)*a2finite)
     & - two*xx*yy*(conjg(a1finite)*a2finite)
     & - two*xx*yy*(conjg(a2finite)*a1finite)

!==== multiply by LO color: Nc*cf
!==== divide by (2*pi)**3
!==== overall factor is as**3*yb**2
!==== NLO^2 hard (e.g. rescaled by LO msq and (as/2/pi)**2)

!==== no overall -1:
!==== for proc1 because we are crossing two fermions from final to initial
!==== for proc2 we are crossing one fermion from final to initial so
!==== there would be a -1 in finitepartsum but then it would be canceled when rescaling by LO msq

      nloxnloeps0=1/(32._dp*pi**2)*yy*zz/(xx**2+one)*nloxnloeps0

      return 
      
   99 format(a35,e20.12)
      
      end

