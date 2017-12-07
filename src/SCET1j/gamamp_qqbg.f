! CW May 16
!     routine for
! ----- q(i1)+qb(i2) => gamma(i3) + g(i4) 
!
! Cross-checked with JC implementation, July 19, 2016
      subroutine gamamp_qqb_gag(order,p,i1,i2,i3,i4,amps0,hard1,hard2)
      implicit none
      include 'types.f'
      include 'constants.f' 
      include 'nf.f' 
      include 'mxpart.f' 
      include 'sprods_com.f' 
      include 'scale.f' 
      include 'ewcouple.f' 
      include 'qcdcouple.f'
      include 'scet_const.f'
      real(dp):: p(mxpart,4)
      integer:: order,i1,i2,i3,i4
      real(dp):: ss,tt,uu
      real(dp):: BigS,BigT,BigU,BigX,BigY,xx,yy
      real(dp) amps0,hard1,hard2,jchard1,jchard2
      complex(dp):: coeff1fin,coeff1C0,Eq13curlybr
! set this to .true. to compare CW and JC expressions
      logical, parameter:: crosscheck=.false.
      amps0=zip
      hard1=zip
      hard2=zip
      call dotem(4,p,s)

!=====Mandelstam invariants 
      ss=s(i1,i2)
      tt=s(i2,i3)
      uu=s(i1,i3)
      
!----- Big X and Big Y defined for physical invariants, so real(dp) is fine 
      BigX=log(-tt/ss)
      BigY=log(-uu/ss)
      BigS=log(ss/musq) 
      BigU=log(-uu/musq)
      BigT=log(-tt/musq)
      xx=-tt/ss
      yy=-uu/ss 

!======LO (matches Nigels)
      if (order >= 0) then
      amps0=8._dp*xn*cf*(tt/uu+uu/tt)
      endif

      if(order >= 1) then 
!-----this is the hard function for 1-loop (2xRe(0x1)+C0)
!-----  Nigel only bit
       hard1=   (2*CF*(6*(1 - 2*BigS + 
     -         (3 + 4*BigS + BigT + BigU - BigX**2 - BigY**2)*
     -          CA**2 + 2*(-2 - 4*BigS + BigX**2 + BigY**2)*CA*CF)*
     -       ss**2 + (6*(3 + BigS - BigS**2) + 7*Pi**2 - 
     -         CA**2*(54 - 6*BigS**2 + 14*BigT + 8*BigU + 
     -            12*BigX*BigY - 6*BigY**2 + Pi**2) + 
     -         2*CA*((BigT + BigU)*NF + 
     -            CF*(12 - 6*BigS**2 + 6*BigX**2 + 7*Pi**2)))*tt**2
     -       + 6*(-12 + BigT + BigU)*CA*(CA - 2*CF)*tt*uu + 
     -      (6*(3 + BigS - BigS**2) + 7*Pi**2 - 
     -         CA**2*(54 - 6*BigS**2 + 8*BigT + 14*BigU - 
     -            6*BigX**2 + 12*BigX*BigY + Pi**2) + 
     -         2*CA*((BigT + BigU)*NF + 
     -            CF*(12 - 6*BigS**2 + 6*BigY**2 + 7*Pi**2)))*uu**2
     -       + 6*CA*ss*(6*CF*(tt - BigU*tt + uu - BigT*uu) + 
     -         CA*((-2 + 2*BigT + 3*BigU)*tt + 
     -        (-2 + 3*BigT + 2*BigU)*uu))))/(3.*tt*uu)
!------ additional piece
       hard1=hard1+
     &      (-2*CA*CF*(CA*(18*BigS - 6*BigS**2 - 20*BigT + 6*BigT**2 - 
     -         20*BigU + 6*BigU**2 + 5*Pi**2) + 
     -      2*((BigT + BigU)*NF + CF*(6*(-3 + BigS)*BigS - 7*Pi**2))
     -      )*(tt**2 + uu**2))/(3.*tt*uu)

      endif

      if(order ==2) then
         hard2=(CF*(2*(9*CF**2*(-6*(-3 + BigS)*BigS + 7*Pi**2)*
     -          ((2*(-6 - 5*BigS - 3*BigS**2 + 6*BigU + 4*BigX**2 + 
     -                  2*BigY**2) + 7*Pi**2)*tt**2 + 
     -            8*(1 - 4*BigS + BigT + BigU + BigX**2 + BigY**2)*
     -             tt*uu + (2*
     -                (-6 - 5*BigS - 3*BigS**2 + 6*BigT + 
     -                  2*BigX**2 + 4*BigY**2) + 7*Pi**2)*uu**2)*xn
     -          - NF*(3*(BigT + BigU)*
     -             ((-6*(-4 + BigS + BigS**2) + 7*Pi**2)*tt**2 + 
     -               12*(1 - 2*BigS)*tt*uu + 
     -               (-6*(-4 + BigS + BigS**2) + 7*Pi**2)*uu**2) + 
     -            NF*(9*BigT**2 + BigT*(-20 + 6*BigU) + 
     -               BigU*(-20 + 9*BigU) + Pi**2)*(tt**2 + uu**2)*xn
     -              + 3*(BigT + BigU)*
     -             ((3*(-4 + 5*BigT + 2*BigT**2 + 5*BigU + 
     -                     2*BigU**2) - 2*Pi**2)*tt**2 + 
     -               12*(1 + 2*BigT + 2*BigU)*tt*uu + 
     -               (3*(-4 + 5*BigT + 2*BigT**2 + 5*BigU + 
     -                     2*BigU**2) - 2*Pi**2)*uu**2)*xn**2) + 
     -         CF*(-21*Pi**4*(tt**2 + uu**2)*(-7 + 2*xn**2) + 
     -            Pi**2*(252*tt*uu*
     -                (1 - 2*BigS + (1 + 2*BigT + 2*BigU)*xn**2) + 
     -               tt**2*(-252*(-2 - BigS + BigS**2) + 
     -                  2*(33*BigS - 7*(16 + 3*BigT + 3*BigU))*NF*
     -                   xn + 
     -                  9*(-12*BigS + 4*BigS**2 + 
     -                     7*
     -                      (-4 + 5*BigT + 2*BigT**2 + 5*BigU + 
     -                       2*BigU**2))*xn**2) + 
     -               uu**2*(-252*(-2 - BigS + BigS**2) + 
     -                  2*(33*BigS - 7*(16 + 3*BigT + 3*BigU))*NF*
     -                   xn + 
     -                  9*(-12*BigS + 4*BigS**2 + 
     -                     7*
     -                      (-4 + 5*BigT + 2*BigT**2 + 5*BigU + 
     -                       2*BigU**2))*xn**2)) + 
     -            6*(18*BigS**4*(tt**2 + uu**2) - 
     -               4*BigS**3*
     -                (-18*tt*uu + tt**2*(9 + NF*xn) + 
     -                  uu**2*(9 + NF*xn)) + 
     -               3*BigS*
     -                (4*tt*uu*
     -                   (9 + 4*(BigT + BigU)*NF*xn + 
     -                     9*(1 + 2*BigT + 2*BigU)*xn**2) + 
     -                  tt**2*
     -                   (72 + 4*(-5 + 2*BigT + 2*BigU)*NF*xn + 
     -                     9*
     -                      (-4 + 5*BigT + 2*BigT**2 + 5*BigU + 
     -                       2*BigU**2)*xn**2) + 
     -                  uu**2*
     -                   (72 + 4*(-5 + 2*BigT + 2*BigU)*NF*xn + 
     -                     9*
     -                      (-4 + 5*BigT + 2*BigT**2 + 5*BigU + 
     -                       2*BigU**2)*xn**2)) - 
     -               BigS**2*
     -                (36*tt*uu*(7 + (1 + 2*BigT + 2*BigU)*xn**2) + 
     -                  tt**2*
     -                   (126 - 2*(19 + 3*BigT + 3*BigU)*NF*xn + 
     -                     9*
     -                      (-4 + 5*BigT + 2*BigT**2 + 5*BigU + 
     -                       2*BigU**2)*xn**2) + 
     -                  uu**2*
     -                   (126 - 2*(19 + 3*BigT + 3*BigU)*NF*xn + 
     -                     9*
     -                      (-4 + 5*BigT + 2*BigT**2 + 5*BigU + 
     -                       2*BigU**2)*xn**2)) - 
     -               6*NF*xn*
     -                (2*(BigT**2 + 
     -                     BigU*(1 + BigU + BigX**2 + BigY**2) + 
     -                     BigT*(1 + 2*BigU + BigX**2 + BigY**2))*
     -                   tt*uu + 
     -                  tt**2*
     -                   (3*BigU**2 + 
     -                     BigU*(-3 + 2*BigX**2 + BigY**2) + 
     -                     BigT*
     -                      (-3 + 3*BigU + 2*BigX**2 + BigY**2) - 
     -                     Zeta3) + 
     -                  uu**2*
     -                   (3*BigT**2 + 
     -                     BigU*(-3 + BigX**2 + 2*BigY**2) + 
     -                     BigT*
     -                      (-3 + 3*BigU + BigX**2 + 2*BigY**2) - 
     -                     Zeta3))))) + 
     -      CA**2*xn*(57*Pi**4*(tt**2 + uu**2) - 
     -         3*Pi**2*((415 - 109*BigS + 66*BigS**2 - 66*BigT - 
     -               66*BigT**2 - 66*BigU - 66*BigU**2 - 
     -               30*BigX**2 - 60*BigX*BigY)*tt**2 + 
     -            60*(-2 + 4*BigS - 3*BigT - 3*BigU - BigX**2 - 
     -               BigY**2)*tt*uu + 
     -            (415 - 109*BigS + 66*BigS**2 - 66*BigT - 
     -               66*BigT**2 - 66*BigU - 66*BigU**2 - 
     -               60*BigX*BigY - 30*BigY**2)*uu**2) + 
     -         2*(36*(-10*BigT + 3*BigT**2 + BigU*(-10 + 3*BigU))*
     -             (2 + 3*BigT + 3*BigU + BigX**2 + BigY**2)*tt*uu
     -             + 81*BigS**4*(tt**2 + uu**2) - 
     -            6*BigS**3*(29*tt**2 - 72*tt*uu + 29*uu**2) + 
     -            9*BigS*((-98 + 125*BigT + 12*BigT**2 + 125*BigU + 
     -                  12*BigU**2 + 18*BigX**2 + 36*BigX*BigY)*
     -                tt**2 + 
     -               4*(67*BigT - 12*BigT**2 + 67*BigU - 
     -                  12*BigU**2 + 9*(2 + BigX**2 + BigY**2))*tt*
     -                uu + (-98 + 125*BigT + 12*BigT**2 + 
     -                  125*BigU + 12*BigU**2 + 36*BigX*BigY + 
     -                  18*BigY**2)*uu**2) - 
     -            3*BigS**2*
     -             ((-62 - 15*BigT + 54*BigT**2 - 15*BigU + 
     -                  54*BigU**2 + 18*BigX**2 + 36*BigX*BigY)*
     -                tt**2 + 
     -               36*(14 + 3*BigT + 3*BigU + BigX**2 + BigY**2)*
     -                tt*uu + 
     -               (-62 - 15*BigT + 54*BigT**2 - 15*BigU + 
     -                  54*BigU**2 + 36*BigX*BigY + 18*BigY**2)*
     -                uu**2) + 
     -            tt**2*(21*BigT**3 + 81*BigT**4 + 21*BigU**3 + 
     -               81*BigU**4 + 
     -               6*BigU**2*(-229 + 9*BigX**2 + 18*BigX*BigY) - 
     -               20*BigU*(-49 + 9*BigX**2 + 18*BigX*BigY) + 
     -               3*BigT**2*
     -                (-458 - 15*BigU + 54*BigU**2 + 18*BigX**2 + 
     -                  36*BigX*BigY) - 
     -               5*BigT*
     -                (-196 + 300*BigU + 9*BigU**2 + 36*BigX**2 + 
     -                  72*BigX*BigY) - 99*Zeta3) + 
     -            uu**2*(21*BigT**3 + 81*BigT**4 + 21*BigU**3 + 
     -               81*BigU**4 + 
     -               6*BigU**2*(-229 + 18*BigX*BigY + 9*BigY**2) - 
     -               20*BigU*(-49 + 18*BigX*BigY + 9*BigY**2) + 
     -               3*BigT**2*
     -                (-458 - 15*BigU + 54*BigU**2 + 36*BigX*BigY + 
     -                  18*BigY**2) - 
     -               5*BigT*
     -                (-196 + 300*BigU + 9*BigU**2 + 72*BigX*BigY + 
     -                  36*BigY**2) - 99*Zeta3))) - 
     -      CA*(3*Pi**4*(tt**2 + uu**2)*
     -          (35 + 151*CF*xn - 10*xn**2) + 
     -         Pi**2*(tt**2*
     -             (-6*(-60 - 48*BigS + 36*BigS**2 + 70*BigT - 
     -                  21*BigT**2 + 70*BigU - 21*BigU**2) + 
     -               2*(-1400 + 579*BigS - 576*BigS**2 - 105*BigT + 
     -                  378*BigT**2 + 165*BigU + 378*BigU**2 + 
     -                  306*BigX**2 + 252*BigX*BigY + 90*BigY**2)*
     -                CF*xn + (-245 + 66*BigS)*NF*xn + 
     -               3*(-60 - 36*BigS + 12*BigS**2 + 115*BigT + 
     -                  18*BigT**2 + 115*BigU + 18*BigU**2)*xn**2)
     -             + uu**2*(-6*
     -                (-60 - 48*BigS + 36*BigS**2 + 70*BigT - 
     -                  21*BigT**2 + 70*BigU - 21*BigU**2) + 
     -               2*(-1400 + 579*BigS - 576*BigS**2 + 165*BigT + 
     -                  378*BigT**2 - 105*BigU + 378*BigU**2 + 
     -                  90*BigX**2 + 252*BigX*BigY + 306*BigY**2)*
     -                CF*xn + (-245 + 66*BigS)*NF*xn + 
     -               3*(-60 - 36*BigS + 12*BigS**2 + 115*BigT + 
     -                  18*BigT**2 + 115*BigU + 18*BigU**2)*xn**2)
     -             - 36*tt*uu*
     -             (2*(-19 + 48*BigS - 26*BigT - 26*BigU - 
     -                  12*BigX**2 - 12*BigY**2)*CF*xn - 
     -               5*(1 - 2*BigS + (1 + 2*BigT + 2*BigU)*xn**2)))
     -          + 2*(54*BigS**4*(tt**2 + uu**2)*(1 + 6*CF*xn) - 
     -            12*BigS**3*
     -             (tt**2*(9 + 47*CF*xn + NF*xn) + 
     -               uu**2*(9 + 47*CF*xn + NF*xn) - 
     -               18*tt*(uu + 8*CF*uu*xn)) + 
     -            9*BigS*(4*tt*uu*
     -                (9 + 20*BigT - 6*BigT**2 + 20*BigU - 
     -                  6*BigU**2 + 
     -                  2*(27 + 76*BigT - 12*BigT**2 + 76*BigU - 
     -                     12*BigU**2 + 18*BigX**2 + 18*BigY**2)*CF*
     -                   xn + 4*(BigT + BigU)*NF*xn + 
     -                  9*(1 + 2*BigT + 2*BigU)*xn**2) + 
     -               tt**2*(72 + 20*BigT - 6*BigT**2 + 20*BigU - 
     -                  6*BigU**2 + 
     -                  2*(-152 + 125*BigT + 12*BigT**2 + 
     -                     179*BigU + 12*BigU**2 + 54*BigX**2 + 
     -                     36*BigX*BigY + 18*BigY**2)*CF*xn + 
     -                  4*(-5 + 2*BigT + 2*BigU)*NF*xn + 
     -                  9*(-4 + 5*BigT + 2*BigT**2 + 5*BigU + 
     -                     2*BigU**2)*xn**2) + 
     -               uu**2*(72 + 20*BigT - 6*BigT**2 + 20*BigU - 
     -                  6*BigU**2 + 
     -                  2*(-152 + 179*BigT + 12*BigT**2 + 
     -                     125*BigU + 12*BigU**2 + 18*BigX**2 + 
     -                     36*BigX*BigY + 54*BigY**2)*CF*xn + 
     -                  4*(-5 + 2*BigT + 2*BigU)*NF*xn + 
     -                  9*(-4 + 5*BigT + 2*BigT**2 + 5*BigU + 
     -                     2*BigU**2)*xn**2)) - 
     -            3*BigS**2*
     -             (36*tt*uu*
     -                (7 + 2*
     -                   (27 + 4*BigT + 4*BigU + 2*BigX**2 + 
     -                     2*BigY**2)*CF*xn + 
     -                  (1 + 2*BigT + 2*BigU)*xn**2) + 
     -               tt**2*(6*
     -                   (21 - 10*BigT + 3*BigT**2 - 10*BigU + 
     -                     3*BigU**2) + 
     -                  2*(19 - 15*BigT + 54*BigT**2 + 39*BigU + 
     -                     54*BigU**2 + 54*BigX**2 + 36*BigX*BigY + 
     -                     18*BigY**2)*CF*xn - 
     -                  2*(19 + 3*BigT + 3*BigU)*NF*xn + 
     -                  9*(-4 + 5*BigT + 2*BigT**2 + 5*BigU + 
     -                     2*BigU**2)*xn**2) + 
     -               uu**2*(6*
     -                   (21 - 10*BigT + 3*BigT**2 - 10*BigU + 
     -                     3*BigU**2) + 
     -                  2*(19 + 39*BigT + 54*BigT**2 - 15*BigU + 
     -                     54*BigU**2 + 18*BigX**2 + 36*BigX*BigY + 
     -                     54*BigY**2)*CF*xn - 
     -                  2*(19 + 3*BigT + 3*BigU)*NF*xn + 
     -                  9*(-4 + 5*BigT + 2*BigT**2 + 5*BigU + 
     -                     2*BigU**2)*xn**2)) + 
     -            36*tt*uu*(6*BigT**3*xn*(CF + xn) + 
     -               BigT**2*
     -                (3 + 2*(-7 + 3*BigU + 3*BigX**2 + 3*BigY**2)*
     -                   CF*xn - 3*NF*xn + (-17 + 6*BigU)*xn**2) + 
     -               BigT*(-10 + 
     -                  2*(-20*BigU + 3*BigU**2 - 
     -                     10*(1 + BigX**2 + BigY**2))*CF*xn - 
     -                  (2 + 6*BigU + BigX**2 + BigY**2)*NF*xn + 
     -                  2*(-5 - 20*BigU + 3*BigU**2)*xn**2) + 
     -               BigU*(-10 + 3*BigU + 
     -                  2*(3*BigU**2 - 10*(1 + BigX**2 + BigY**2) + 
     -                     BigU*(-7 + 3*BigX**2 + 3*BigY**2))*CF*xn
     -                   - (2 + 3*BigU + BigX**2 + BigY**2)*NF*xn + 
     -                  (-10 - 17*BigU + 6*BigU**2)*xn**2)) + 
     -            tt**2*(54*BigT**4*xn**2 + 54*BigU**4*xn**2 - 
     -               3*BigT**3*xn*(2*NF + 15*xn) - 
     -               3*BigU**3*xn*(-108*CF + 2*NF + 15*xn) - 
     -               2*BigU*
     -                (360 + 180*(-3 + 2*BigX**2 + BigY**2)*CF*xn + 
     -                  (-149 + 9*BigX**2 + 18*BigX*BigY)*NF*xn - 
     -                  180*xn**2) + 
     -               6*BigU**2*
     -                (36 + 18*(-13 + 2*BigX**2 + BigY**2)*CF*xn - 
     -                  43*NF*xn - 93*xn**2) - 
     -               BigT*(720 + 
     -                  360*(-3 + 3*BigU + 2*BigX**2 + BigY**2)*CF*
     -                   xn + 
     -                  2*(-149 + 105*BigU + 9*BigU**2 + 
     -                     9*BigX**2 + 18*BigX*BigY)*NF*xn + 
     -                  45*(-8 + 20*BigU + BigU**2)*xn**2) + 
     -               3*BigT**2*
     -                (72 + 
     -                  36*(-3 + 3*BigU + 2*BigX**2 + BigY**2)*CF*
     -                   xn - 2*(43 + 3*BigU)*NF*xn + 
     -                  3*(-62 - 5*BigU + 12*BigU**2)*xn**2) + 
     -               18*(11*CF - NF)*xn*Zeta3) + 
     -            uu**2*(54*BigT**4*xn**2 + 54*BigU**4*xn**2 - 
     -               3*BigU**3*xn*(2*NF + 15*xn) - 
     -               3*BigT**3*xn*(-108*CF + 2*NF + 15*xn) - 
     -               2*BigU*
     -                (360 + 180*(-3 + BigX**2 + 2*BigY**2)*CF*xn + 
     -                  (-149 + 18*BigX*BigY + 9*BigY**2)*NF*xn - 
     -                  180*xn**2) + 
     -               6*BigU**2*
     -                (36 + 18*(-3 + BigX**2 + 2*BigY**2)*CF*xn - 
     -                  43*NF*xn - 93*xn**2) + 
     -               BigT*(-720 + 
     -                  36*(-30*BigU + 9*BigU**2 - 
     -                     10*(-3 + BigX**2 + 2*BigY**2))*CF*xn - 
     -                  2*(-149 + 105*BigU + 9*BigU**2 + 
     -                     18*BigX*BigY + 9*BigY**2)*NF*xn - 
     -                  45*(-8 + 20*BigU + BigU**2)*xn**2) + 
     -               3*BigT**2*
     -                (72 + 36*(-13 + BigX**2 + 2*BigY**2)*CF*xn - 
     -                  2*(43 + 3*BigU)*NF*xn + 
     -                  3*(-62 - 5*BigU + 12*BigU**2)*xn**2) + 
     -               18*(11*CF - NF)*xn*Zeta3)))))/(54.*tt*uu)
      endif

      if (crosscheck) then
c--- JC implementation
      if (order >= 1) then
      coeff1fin =  + im*pi * (
     &     - 6.D0*CF*xx**(-1)*yy
     &     - 8.D0*CF
     &     - 6.D0*CF*xx*yy**(-1)
     &     - 8.D0*CF*BigY*xx**(-1)*yy
     &     - 8.D0*CF*BigY
     &     - 4.D0*CF*BigY*xx*yy**(-1)
     &     - 4.D0*CF*BigX*xx**(-1)*yy
     &     - 8.D0*CF*BigX
     &     - 8.D0*CF*BigX*xx*yy**(-1)
     &     )
      coeff1fin = coeff1fin - 2.D0/3.D0*TR*nf*BigU*xx**(-1)*yy
     &     - 2.D0/3.D0*TR*nf*BigU*xx*yy**(-1)
     &     - 2.D0/3.D0*TR*nf*BigT*xx**(-1)*yy
     &     - 2.D0/3.D0*TR*nf*BigT*xx*yy**(-1)
     &     + 7.D0*CA*xx**(-1)*yy
     &     + 7.D0*CA*xx*yy**(-1)
     &     + 10.D0/3.D0*CA*BigU*xx**(-1)*yy
     &     + 1.D0/3.D0*CA*BigU*xx*yy**(-1)
     &     + 1.D0/3.D0*CA*BigT*xx**(-1)*yy
     &     + 10.D0/3.D0*CA*BigT*xx*yy**(-1)
     &     - CA*BigY**2*xx**(-1)*yy
     &     - CA*BigY**2*xx*yy**(-1)
     &     + 2.D0*CA*BigX*BigY*xx**(-1)*yy
     &     + 2.D0*CA*BigX*BigY*xx*yy**(-1)
     &     - CA*BigX**2*xx**(-1)*yy
     &
      coeff1fin = coeff1fin - CA*BigX**2*xx*yy**(-1)
     &     - CA*pisq*xx**(-1)*yy
     &     - CA*pisq*xx*yy**(-1)
     &     + 2.D0/3.D0*CA**2*TR*nf*BigU*xx**(-1)*yy
     &     + 2.D0/3.D0*CA**2*TR*nf*BigU*xx*yy**(-1)
     &     + 2.D0/3.D0*CA**2*TR*nf*BigT*xx**(-1)*yy
     &     + 2.D0/3.D0*CA**2*TR*nf*BigT*xx*yy**(-1)
     &     - 7.D0*CA**3*xx**(-1)*yy
     &     - 7.D0*CA**3*xx*yy**(-1)
     &     - 10.D0/3.D0*CA**3*BigU*xx**(-1)*yy
     &     - 1.D0/3.D0*CA**3*BigU*xx*yy**(-1)
     &     - 1.D0/3.D0*CA**3*BigT*xx**(-1)*yy
     &     - 10.D0/3.D0*CA**3*BigT*xx*yy**(-1)
     &     + CA**3*BigY**2*xx**(-1)*yy
     &     + CA**3*BigY**2*xx*yy**(-1)
     &
      coeff1fin = coeff1fin - 2.D0*CA**3*BigX*BigY*xx**(-1)*yy
     &     - 2.D0*CA**3*BigX*BigY*xx*yy**(-1)
     &     + CA**3*BigX**2*xx**(-1)*yy
     &     + CA**3*BigX**2*xx*yy**(-1)
     &     + CA**3*pisq*xx**(-1)*yy
     &     + CA**3*pisq*xx*yy**(-1)
     &     + 14.D0*CF*xx**(-1)*yy
     &     + 14.D0*CF*xx*yy**(-1)
     &     - 4.D0*CF*BigU
     &     - 6.D0*CF*BigU*xx*yy**(-1)
     &     - 6.D0*CF*BigT*xx**(-1)*yy
     &     - 4.D0*CF*BigT
     &     + 6.D0*CF*BigS*xx**(-1)*yy
     &     + 8.D0*CF*BigS
     &     + 6.D0*CF*BigS*xx*yy**(-1)
     &
      coeff1fin = coeff1fin - 4.D0*CF*BigY**2*xx**(-1)*yy
     &     - 4.D0*CF*BigY**2
     &     - 2.D0*CF*BigY**2*xx*yy**(-1)
     &     - 2.D0*CF*BigX**2*xx**(-1)*yy
     &     - 4.D0*CF*BigX**2
     &     - 4.D0*CF*BigX**2*xx*yy**(-1)
     &

      coeff1C0 =  + im*pi * (
     &     + 3.D0/4.D0*CA
     &     - 1.D0/2.D0*CA*BigS
     &     - 3.D0/2.D0*CF
     &     + CF*BigS
     &     )
      coeff1C0 = coeff1C0 - 1.D0/6.D0*TR*nf*BigU
     &     - 1.D0/6.D0*TR*nf*BigT
     &     + 5.D0/6.D0*CA*BigU
     &     - 1.D0/4.D0*CA*BigU**2
     &     + 5.D0/6.D0*CA*BigT
     &     - 1.D0/4.D0*CA*BigT**2
     &     - 3.D0/4.D0*CA*BigS
     &     + 1.D0/4.D0*CA*BigS**2
     &     - 5.D0/24.D0*CA*pisq
     &     + 3.D0/2.D0*CF*BigS
     &     - 1.D0/2.D0*CF*BigS**2
     &     + 7.D0/12.D0*CF*pisq

! factor of two for interference
      jchard1=two*(coeff1fin+coeff1C0*amps0)
      
      write(6,*) 'qqb: hard1,jchard1,hard1/jchard1',hard1,jchard1,hard1/jchard1
!      pause
      
      endif
      
      if (order == 2) then
      
      Eq13curlybr =  + im*pi * (
     &     + 67.D0/24.D0
     &     + 5.D0/8.D0*BigU
     &     - 3.D0/16.D0*BigU**2
     &     + 5.D0/8.D0*BigT
     &     - 3.D0/16.D0*BigT**2
     &     - 233.D0/72.D0*BigS
     &     - 5.D0/12.D0*BigS*BigU
     &     + 1.D0/8.D0*BigS*BigU**2
     &     - 5.D0/12.D0*BigS*BigT
     &     + 1.D0/8.D0*BigS*BigT**2
     &     + 11.D0/24.D0*BigS**2
     &     - 17.D0/96.D0*pisq
     &     + 1.D0/24.D0*pisq*BigS
     &     )
      Eq13curlybr = Eq13curlybr + CA * (
     &     - 167.D0/108.D0*TR*nf*BigU
     &     + 41.D0/72.D0*TR*nf*BigU**2
     &     - 1.D0/72.D0*TR*nf*BigU**3
     &     - 167.D0/108.D0*TR*nf*BigT
     &     - 5.D0/18.D0*TR*nf*BigT*BigU
     &     + 1.D0/24.D0*TR*nf*BigT*BigU**2
     &     + 41.D0/72.D0*TR*nf*BigT**2
     &     + 1.D0/24.D0*TR*nf*BigT**2*BigU
     &     - 1.D0/72.D0*TR*nf*BigT**3
     &     + 5.D0/6.D0*TR*nf*BigS
     &     + 1.D0/8.D0*TR*nf*BigS*BigU
     &     + 1.D0/8.D0*TR*nf*BigS*BigT
     &     - 19.D0/36.D0*TR*nf*BigS**2
     &     - 1.D0/24.D0*TR*nf*BigS**2*BigU
     &     )
      Eq13curlybr = Eq13curlybr + CA * (
     &     - 1.D0/24.D0*TR*nf*BigS**2*BigT
     &     + 1.D0/18.D0*TR*nf*BigS**3
     &     + 1.D0/12.D0*TR*nf*zeta3
     &     + 245.D0/432.D0*TR*nf*pisq
     &     + 7.D0/144.D0*TR*nf*pisq*BigU
     &     + 7.D0/144.D0*TR*nf*pisq*BigT
     &     - 11.D0/72.D0*TR*nf*pisq*BigS
     &     )
      Eq13curlybr = Eq13curlybr + CA*im*pi * (
     &     - 5.D0/6.D0*TR*nf
     &     - 1.D0/8.D0*TR*nf*BigU
     &     - 1.D0/8.D0*TR*nf*BigT
     &     + 19.D0/18.D0*TR*nf*BigS
     &     + 1.D0/12.D0*TR*nf*BigS*BigU
     &     + 1.D0/12.D0*TR*nf*BigS*BigT
     &     - 1.D0/6.D0*TR*nf*BigS**2
     &     + 1.D0/24.D0*TR*nf*pisq
     &     )
      Eq13curlybr = Eq13curlybr + CA**2 * (
     &     + 335.D0/108.D0*BigU
     &     - 97.D0/72.D0*BigU**2
     &     - 1.D0/18.D0*BigU**3
     &     + 1.D0/32.D0*BigU**4
     &     + 335.D0/108.D0*BigT
     &     + 25.D0/36.D0*BigT*BigU
     &     - 5.D0/24.D0*BigT*BigU**2
     &     - 97.D0/72.D0*BigT**2
     &     - 5.D0/24.D0*BigT**2*BigU
     &     + 1.D0/16.D0*BigT**2*BigU**2
     &     - 1.D0/18.D0*BigT**3
     &     + 1.D0/32.D0*BigT**4
     &     + 9.D0/32.D0*BigS**2
     &     - 3.D0/16.D0*BigS**3
     &     )
      Eq13curlybr = Eq13curlybr + CA**2 * (
     &     + 1.D0/32.D0*BigS**4
     &     - 11.D0/24.D0*zeta3
     &     - 143.D0/432.D0*pisq
     &     - 1.D0/32.D0*pisq*BigU
     &     + 1.D0/48.D0*pisq*BigU**2
     &     - 1.D0/32.D0*pisq*BigT
     &     + 1.D0/48.D0*pisq*BigT**2
     &     + 19.D0/32.D0*pisq*BigS
     &     - 19.D0/96.D0*pisq*BigS**2
     &     + 49.D0/1152.D0*pisq**2
     &     )
      Eq13curlybr = Eq13curlybr + CA**2*im*pi * (
     &     - 9.D0/16.D0*BigS
     &     + 9.D0/16.D0*BigS**2
     &     - 1.D0/8.D0*BigS**3
     &     - 7.D0/32.D0*pisq
     &     + 7.D0/48.D0*pisq*BigS
     &     )
      Eq13curlybr = Eq13curlybr + CF * (
     &     - 5.D0/3.D0*TR*nf*BigS
     &     - 1.D0/4.D0*TR*nf*BigS*BigU
     &     - 1.D0/4.D0*TR*nf*BigS*BigT
     &     + 19.D0/18.D0*TR*nf*BigS**2
     &     + 1.D0/12.D0*TR*nf*BigS**2*BigU
     &     + 1.D0/12.D0*TR*nf*BigS**2*BigT
     &     - 1.D0/9.D0*TR*nf*BigS**3
     &     + 1.D0/6.D0*TR*nf*zeta3
     &     - 28.D0/27.D0*TR*nf*pisq
     &     - 7.D0/72.D0*TR*nf*pisq*BigU
     &     - 7.D0/72.D0*TR*nf*pisq*BigT
     &     + 11.D0/36.D0*TR*nf*pisq*BigS
     &     )
      Eq13curlybr = Eq13curlybr + CF*im*pi * (
     &     + 5.D0/3.D0*TR*nf
     &     + 1.D0/4.D0*TR*nf*BigU
     &     + 1.D0/4.D0*TR*nf*BigT
     &     - 19.D0/9.D0*TR*nf*BigS
     &     - 1.D0/6.D0*TR*nf*BigS*BigU
     &     - 1.D0/6.D0*TR*nf*BigS*BigT
     &     + 1.D0/3.D0*TR*nf*BigS**2
     &     - 1.D0/12.D0*TR*nf*pisq
     &     )
      Eq13curlybr = Eq13curlybr + CF*CA * (
     &     - 9.D0/8.D0*BigS**2
     &     + 3.D0/4.D0*BigS**3
     &     - 1.D0/8.D0*BigS**4
     &     + 9.D0/8.D0*pisq
     &     - 19.D0/8.D0*pisq*BigS
     &     + 19.D0/24.D0*pisq*BigS**2
     &     - 49.D0/288.D0*pisq**2
     &     )
      Eq13curlybr = Eq13curlybr + CF*CA*im*pi * (
     &     + 9.D0/4.D0*BigS
     &     - 9.D0/4.D0*BigS**2
     &     + 1.D0/2.D0*BigS**3
     &     + 7.D0/8.D0*pisq
     &     - 7.D0/12.D0*pisq*BigS
     &     )
      Eq13curlybr = Eq13curlybr + CF**2 * (
     &     + 9.D0/8.D0*BigS**2
     &     - 3.D0/4.D0*BigS**3
     &     + 1.D0/8.D0*BigS**4
     &     - 9.D0/8.D0*pisq
     &     + 19.D0/8.D0*pisq*BigS
     &     - 19.D0/24.D0*pisq*BigS**2
     &     + 49.D0/288.D0*pisq**2
     &     )
      Eq13curlybr = Eq13curlybr + CF**2*im*pi * (
     &     - 9.D0/4.D0*BigS
     &     + 9.D0/4.D0*BigS**2
     &     - 1.D0/2.D0*BigS**3
     &     - 7.D0/8.D0*pisq
     &     + 7.D0/12.D0*pisq*BigS
     &     )
      Eq13curlybr = Eq13curlybr - 67.D0/24.D0*BigS
     &     - 5.D0/8.D0*BigS*BigU
     &     + 3.D0/16.D0*BigS*BigU**2
     &     - 5.D0/8.D0*BigS*BigT
     &     + 3.D0/16.D0*BigS*BigT**2
     &     + 233.D0/144.D0*BigS**2
     &     + 5.D0/24.D0*BigS**2*BigU
     &     - 1.D0/16.D0*BigS**2*BigU**2
     &     + 5.D0/24.D0*BigS**2*BigT
     &     - 1.D0/16.D0*BigS**2*BigT**2
     &     - 11.D0/72.D0*BigS**3
     &     + 11.D0/48.D0*zeta3
     &     - 691.D0/432.D0*pisq
     &     - 35.D0/144.D0*pisq*BigU
     &     + 7.D0/96.D0*pisq*BigU**2
     &
      Eq13curlybr = Eq13curlybr - 35.D0/144.D0*pisq*BigT
     &     + 7.D0/96.D0*pisq*BigT**2
     &     + 139.D0/288.D0*pisq*BigS
     &     - 1.D0/48.D0*pisq*BigS**2
     &     + 11.D0/576.D0*pisq**2
     &     + 5.D0/27.D0*TR**2*nf**2*BigU
     &     - 1.D0/24.D0*TR**2*nf**2*BigU**2
     &     + 5.D0/27.D0*TR**2*nf**2*BigT
     &     + 1.D0/36.D0*TR**2*nf**2*BigT*BigU
     &     - 1.D0/24.D0*TR**2*nf**2*BigT**2
     &     - 1.D0/108.D0*TR**2*nf**2*pisq

! factor of two for interferences
      jchard2=two*(coeff1C0*coeff1fin+Eq13curlybr*amps0)

! extra pieces, c.f. gagaAGTYassemble1loop
      jchard2=jchard2+two*coeff1C0*conjg(coeff1fin)
     & +coeff1C0*conjg(coeff1C0)*amps0

!      write(6,*) '2-loop bit',two*coeff1C0*coeff1fin+two*Eq13curlybr*amps0
!      write(6,*) '1-loop bit',two*coeff1C0*conjg(coeff1fin)+coeff1C0*conjg(coeff1C0)*amps0
      write(6,*) 'qqb: hard2,jchard2,hard2/jchard2',hard2,jchard2,hard2/jchard2
!      pause
      
      endif
      endif
      
      return
      end

