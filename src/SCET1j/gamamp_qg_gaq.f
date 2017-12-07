!     routine for
! ----- q(i1)+g(i2) => gamma(i3) + q(i4) 
!
! Cross-checked with JC implementation, July 20, 2016
      subroutine gamamp_qg_gaq(order,p,i1,i2,i3,i4,amps0,hard1,hard2)
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
      real(dp) amps0,hard1,hard2,xtag,jchard1,jchard2
      complex(dp):: coeff1fin,coeff1C0,Eq13curlybr
! set this to .true. to compare CW and JC expressions
      logical, parameter:: crosscheck=.false.
      amps0=zip
      hard1=zip
      hard2=zip
!==== this controls sign of im in box crossings,
      xtag=1._dp
!=====xtag = -1 swaps s<-> u in box
!      xtag=-1._dp
      
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
      amps0=8._dp*xn*cf*(ss/uu+uu/ss)
      endif

      if(order >= 1) then 
!-----this is the hard function for 1-loop (2xRe(0x1)+C0)
       hard1= (2*CF*((6*(3 + BigT - BigT**2) + Pi**2 + 
     -         4*CA*CF*(6 + 9*BigT - 6*BigT**2 + 3*BigX**2 + 
     -            Pi**2) + 6*CA**2*
     -          (-9 + BigS - BigS**2 - 3*BigT + 2*BigT**2 + 
     -            2*BigU - BigU**2 - BigX**2 + BigY**2 + Pi**2))*
     -       ss**2 - 6*(-1 + 2*BigT - 
     -         2*CA*CF*(-2 - 4*BigT + 2*BigX**2 - 2*BigX*BigY + 
     -            BigY**2 + Pi**2) + 
     -         CA**2*(-3 - BigS - 4*BigT - BigU + 2*BigX**2 - 
     -            2*BigX*BigY + BigY**2 + Pi**2))*tt**2 + 
     -      6*CA*((-2 + 3*BigS + 2*BigU)*CA - 6*(-1 + BigS)*CF)*tt*
     -       uu + (6*(3 + BigT - BigT**2) - 
     -         6*(9 - 2*BigS + BigS**2 + 3*BigT - 2*BigT**2 - 
     -            BigU + BigU**2 + BigX**2 - 2*BigX*BigY)*CA**2 + 
     -         Pi**2 + 4*CA*CF*
     -          (3*(2 + 3*BigT - 2*BigT**2 + BigX**2 - 
     -               2*BigX*BigY + BigY**2) + 4*Pi**2))*uu**2 + 
     -      6*CA*ss*(-2*CF*(3*(-1 + BigU)*tt + 
     -            (-12 + BigS + BigU)*uu) + 
     -         CA*((-2 + 2*BigS + 3*BigU)*tt + 
     -            (-12 + BigS + BigU)*uu))))/(3.*ss*uu)

      endif

      if(order==2) then
         hard2=(CF*(CA**2*xn*(3*Pi**4*(19*ss**2 - 84*ss*uu - 23*uu**2) - 
     -         3*Pi**2*((-395 + 246*BigS + 78*BigS**2 - 49*BigT - 
     -               78*BigT**2 + 114*BigU + 78*BigU**2 + 
     -               126*BigX**2 - 84*BigX*BigY)*ss**2 - 
     -            12*(-14 - 41*BigS + 6*BigS**2 + 46*BigT - 
     -               6*BigT**2 - 41*BigU + 6*BigU**2 - 14*BigX**2 + 
     -               14*BigX*BigY - 7*BigY**2)*ss*uu + 
     -            (-395 + 366*BigS + 42*BigS**2 - 157*BigT - 
     -               42*BigT**2 + 234*BigU + 42*BigU**2 + 
     -               126*BigX**2 - 168*BigX*BigY + 42*BigY**2)*uu**2
     -            ) + 2*(36*
     -             (-10*BigS + 3*BigS**2 + 9*BigT - 3*BigT**2 + 
     -               BigU*(-10 + 3*BigU))*
     -             (2 + 3*BigS - 4*BigT + 3*BigU + 2*BigX**2 - 
     -               2*BigX*BigY + BigY**2)*ss*uu + 
     -            ss**2*(21*BigS**3 + 81*BigS**4 - 174*BigT**3 + 
     -               81*BigT**4 + 980*BigU - 1374*BigU**2 + 
     -               21*BigU**3 + 81*BigU**4 - 540*BigU*BigX**2 + 
     -               162*BigU**2*BigX**2 + 360*BigU*BigX*BigY - 
     -               108*BigU**2*BigX*BigY + 
     -               9*BigT*
     -                (-98 + 125*BigU + 12*BigU**2 + 54*BigX**2 - 
     -                  36*BigX*BigY) - 
     -               3*BigT**2*
     -                (-62 - 15*BigU + 54*BigU**2 + 54*BigX**2 - 
     -                  36*BigX*BigY) - 
     -               3*BigS**2*
     -                (458 - 36*BigT + 54*BigT**2 + 15*BigU - 
     -                  54*BigU**2 - 54*BigX**2 + 36*BigX*BigY) + 
     -               5*BigS*
     -                (196 + 225*BigT + 9*BigT**2 - 300*BigU - 
     -                  9*BigU**2 - 108*BigX**2 + 72*BigX*BigY) - 
     -               99*Zeta3) + 
     -            uu**2*(21*BigS**3 + 81*BigS**4 - 174*BigT**3 + 
     -               81*BigT**4 + 980*BigU - 1374*BigU**2 + 
     -               21*BigU**3 + 81*BigU**4 - 540*BigU*BigX**2 + 
     -               162*BigU**2*BigX**2 + 720*BigU*BigX*BigY - 
     -               216*BigU**2*BigX*BigY - 180*BigU*BigY**2 + 
     -               54*BigU**2*BigY**2 + 
     -               5*BigS*
     -                (196 + 225*BigT + 9*BigT**2 - 300*BigU - 
     -                  9*BigU**2 - 108*BigX**2 + 144*BigX*BigY - 
     -                  36*BigY**2) - 
     -               3*BigS**2*
     -                (458 - 36*BigT + 54*BigT**2 + 15*BigU - 
     -                  54*BigU**2 - 54*BigX**2 + 72*BigX*BigY - 
     -                  18*BigY**2) + 
     -               9*BigT*
     -                (-98 + 125*BigU + 12*BigU**2 + 54*BigX**2 - 
     -                  72*BigX*BigY + 18*BigY**2) - 
     -               3*BigT**2*
     -                (-62 - 15*BigU + 54*BigU**2 + 54*BigX**2 - 
     -                  72*BigX*BigY + 18*BigY**2) - 99*Zeta3))) + 
     -      CA*(3*Pi**4*(144*CF*ss*uu*xn + 
     -            ss**2*(7 + 113*CF*xn - 56*xn**2) + 
     -            uu**2*(7 + 185*CF*xn - 56*xn**2)) + 
     -         Pi**2*(-36*ss*uu*
     -             (-7 + 14*BigT + 
     -               2*(-5 - 24*BigS + 6*BigS**2 + 60*BigT - 
     -                  12*BigT**2 - 24*BigU + 6*BigU**2 - 
     -                  12*BigX**2 + 12*BigX*BigY - 6*BigY**2)*CF*xn
     -                 - 2*(BigS + BigU)*NF*xn - 
     -               7*(1 + 2*BigS + 2*BigU)*xn**2) - 
     -            ss**2*(6*(-84 - 10*BigS + 3*BigS**2 + 30*BigT + 
     -                  18*BigT**2 - 10*BigU + 3*BigU**2) + 
     -               2*(430 - 375*BigS + 162*BigS**2 + 507*BigT + 
     -                  72*BigT**2 - 753*BigU + 162*BigU**2 - 
     -                  324*BigX**2 + 216*BigX*BigY - 126*BigY**2)*
     -                CF*xn + (289 - 72*BigS - 6*BigT)*NF*xn - 
     -               3*(-84 - 55*BigS + 90*BigS**2 + 144*BigT - 
     -                  48*BigT**2 - 55*BigU + 90*BigU**2)*xn**2) + 
     -            uu**2*(-6*
     -                (-84 - 10*BigS + 3*BigS**2 + 30*BigT + 
     -                  18*BigT**2 - 10*BigU + 3*BigU**2) + 
     -               2*(-430 + 1113*BigS - 270*BigS**2 - 
     -                  1155*BigT + 144*BigT**2 + 735*BigU - 
     -                  270*BigU**2 + 324*BigX**2 - 432*BigX*BigY + 
     -                  234*BigY**2)*CF*xn + 
     -               (-289 + 108*BigS + 6*BigT + 36*BigU)*NF*xn + 
     -               3*(-84 - 55*BigS + 90*BigS**2 + 144*BigT - 
     -                  48*BigT**2 - 55*BigU + 90*BigU**2)*xn**2))
     -          - 2*(36*ss*uu*
     -             (6*BigS**3*xn*(CF + xn) + 
     -               6*BigT**3*(1 + 8*CF*xn) - 
     -               3*BigT**2*
     -                (7 + 2*
     -                   (27 + 4*BigU + 4*BigX**2 - 4*BigX*BigY + 
     -                     2*BigY**2)*CF*xn + (1 + 2*BigU)*xn**2) + 
     -               BigT*(9 + 20*BigU - 6*BigU**2 + 
     -                  2*(76*BigU - 12*BigU**2 + 
     -                     9*
     -                      (3 + 4*BigX**2 - 4*BigX*BigY + 
     -                       2*BigY**2))*CF*xn + 4*BigU*NF*xn + 
     -                  9*(1 + 2*BigU)*xn**2) + 
     -               BigS**2*
     -                (3 - 6*BigT - 
     -                  2*(7 + 12*BigT - 3*BigU - 6*BigX**2 + 
     -                     6*BigX*BigY - 3*BigY**2)*CF*xn - 
     -                  3*NF*xn + (-17 + 6*BigU)*xn**2) - 
     -               BigS*(10 - 20*BigT + 
     -                  2*(-76*BigT + 12*BigT**2 + 20*BigU - 
     -                     3*BigU**2 + 
     -                     10*
     -                      (1 + 2*BigX**2 - 2*BigX*BigY + BigY**2))
     -                    *CF*xn + 
     -                  (2 - 4*BigT + 6*BigU + 2*BigX**2 - 
     -                     2*BigX*BigY + BigY**2)*NF*xn + 
     -                  2*(5 - 9*BigT + 3*BigT**2 + 20*BigU - 
     -                     3*BigU**2)*xn**2) + 
     -               BigU*(-10 + 3*BigU + 
     -                  2*(3*BigU**2 - 
     -                     10*
     -                      (1 + 2*BigX**2 - 2*BigX*BigY + BigY**2)
     -                      + BigU*
     -                      (-7 + 6*BigX**2 - 6*BigX*BigY + 
     -                       3*BigY**2))*CF*xn - 
     -                  (2 + 3*BigU + 2*BigX**2 - 2*BigX*BigY + 
     -                     BigY**2)*NF*xn + 
     -                  (-10 - 17*BigU + 6*BigU**2)*xn**2)) + 
     -            ss**2*(-720*BigU + 216*BigU**2 + 
     -               1080*BigU*CF*xn - 1404*BigU**2*CF*xn + 
     -               324*BigU**3*CF*xn - 1080*BigU*BigX**2*CF*xn + 
     -               324*BigU**2*BigX**2*CF*xn + 
     -               720*BigU*BigX*BigY*CF*xn - 
     -               216*BigU**2*BigX*BigY*CF*xn - 
     -               360*BigU*BigY**2*CF*xn + 
     -               108*BigU**2*BigY**2*CF*xn + 298*BigU*NF*xn - 
     -               258*BigU**2*NF*xn - 6*BigU**3*NF*xn - 
     -               54*BigU*BigX**2*NF*xn + 
     -               36*BigU*BigX*BigY*NF*xn + 54*BigS**4*xn**2 + 
     -               360*BigU*xn**2 - 558*BigU**2*xn**2 - 
     -               45*BigU**3*xn**2 + 54*BigU**4*xn**2 - 
     -               3*BigS**3*xn*(2*NF + 15*xn) + 
     -               54*BigT**4*(1 + 6*CF*xn) - 
     -               12*BigT**3*(9 + 47*CF*xn + NF*xn) - 
     -               3*BigS**2*
     -                (18*(-4 + BigT + BigT**2) + 
     -                  36*(3 - 2*BigT + 3*BigT**2 - 3*BigU - 
     -                     3*BigX**2 + 2*BigX*BigY - BigY**2)*CF*xn
     -                   + 2*(43 + 3*BigU)*NF*xn + 
     -                  3*(62 - 18*BigT + 6*BigT**2 + 5*BigU - 
     -                     12*BigU**2)*xn**2) + 
     -               BigS*(180*(-4 + BigT + BigT**2) + 
     -                  90*(25*BigT + BigT**2 - 
     -                     4*
     -                      (-3 + 3*BigU + 3*BigX**2 - 
     -                       2*BigX*BigY + BigY**2))*CF*xn + 
     -                  2*(149 + 36*BigT + 9*BigT**2 - 105*BigU - 
     -                     9*BigU**2 - 27*BigX**2 + 18*BigX*BigY)*
     -                   NF*xn - 
     -                  45*(-8 - 9*BigT + 3*BigT**2 + 20*BigU + 
     -                     BigU**2)*xn**2) + 
     -               9*BigT*
     -                (72 + 20*BigU - 6*BigU**2 + 
     -                  2*(179*BigU + 12*BigU**2 + 
     -                     2*
     -                      (-76 + 54*BigX**2 - 36*BigX*BigY + 
     -                       9*BigY**2))*CF*xn + 
     -                  4*(-5 + 2*BigU)*NF*xn + 
     -                  9*(-4 + 5*BigU + 2*BigU**2)*xn**2) - 
     -               3*BigT**2*
     -                (6*(21 - 10*BigU + 3*BigU**2) + 
     -                  2*(19 + 39*BigU + 54*BigU**2 + 
     -                     108*BigX**2 - 72*BigX*BigY + 18*BigY**2)*
     -                   CF*xn - 2*(19 + 3*BigU)*NF*xn + 
     -                  9*(-4 + 5*BigU + 2*BigU**2)*xn**2) + 
     -               198*CF*xn*Zeta3 - 18*NF*xn*Zeta3) + 
     -            uu**2*(-720*BigU + 216*BigU**2 + 
     -               1080*BigU*CF*xn - 324*BigU**2*CF*xn - 
     -               1080*BigU*BigX**2*CF*xn + 
     -               324*BigU**2*BigX**2*CF*xn + 
     -               1440*BigU*BigX*BigY*CF*xn - 
     -               432*BigU**2*BigX*BigY*CF*xn - 
     -               720*BigU*BigY**2*CF*xn + 
     -               216*BigU**2*BigY**2*CF*xn + 298*BigU*NF*xn - 
     -               258*BigU**2*NF*xn - 6*BigU**3*NF*xn - 
     -               54*BigU*BigX**2*NF*xn + 
     -               72*BigU*BigX*BigY*NF*xn - 
     -               18*BigU*BigY**2*NF*xn + 54*BigS**4*xn**2 + 
     -               360*BigU*xn**2 - 558*BigU**2*xn**2 - 
     -               45*BigU**3*xn**2 + 54*BigU**4*xn**2 - 
     -               3*BigS**3*xn*(-108*CF + 2*NF + 15*xn) + 
     -               54*BigT**4*(1 + 6*CF*xn) - 
     -               12*BigT**3*(9 + 47*CF*xn + NF*xn) - 
     -               3*BigS**2*
     -                (18*(-4 + BigT + BigT**2) + 
     -                  36*(13 - 2*BigT + 3*BigT**2 - 3*BigX**2 + 
     -                     4*BigX*BigY - 2*BigY**2)*CF*xn + 
     -                  2*(43 + 3*BigU)*NF*xn + 
     -                  3*(62 - 18*BigT + 6*BigT**2 + 5*BigU - 
     -                     12*BigU**2)*xn**2) + 
     -               BigS*(180*(-4 + BigT + BigT**2) - 
     -                  18*(-60 - 179*BigT + 13*BigT**2 + 60*BigU - 
     -                     18*BigU**2 + 60*BigX**2 - 80*BigX*BigY + 
     -                     40*BigY**2)*CF*xn + 
     -                  2*(149 + 36*BigT + 9*BigT**2 - 105*BigU - 
     -                     9*BigU**2 - 27*BigX**2 + 36*BigX*BigY - 
     -                     9*BigY**2)*NF*xn - 
     -                  45*(-8 - 9*BigT + 3*BigT**2 + 20*BigU + 
     -                     BigU**2)*xn**2) + 
     -               9*BigT*
     -                (72 + 20*BigU - 6*BigU**2 + 
     -                  2*(125*BigU + 12*BigU**2 + 
     -                     2*
     -                      (-76 + 54*BigX**2 - 72*BigX*BigY + 
     -                       27*BigY**2))*CF*xn + 
     -                  4*(-5 + 2*BigU)*NF*xn + 
     -                  9*(-4 + 5*BigU + 2*BigU**2)*xn**2) - 
     -               3*BigT**2*
     -                (6*(21 - 10*BigU + 3*BigU**2) + 
     -                  2*(19 - 15*BigU + 54*BigU**2 + 
     -                     108*BigX**2 - 144*BigX*BigY + 54*BigY**2)
     -                    *CF*xn - 2*(19 + 3*BigU)*NF*xn + 
     -                  9*(-4 + 5*BigU + 2*BigU**2)*xn**2) + 
     -               198*CF*xn*Zeta3 - 18*NF*xn*Zeta3))) + 
     -      2*(9*CF**2*(-6*(-3 + BigT)*BigT + Pi**2)*
     -          ((-2*(6 + 5*BigT + 3*BigT**2 - 6*BigU - 6*BigX**2 + 
     -                  4*BigX*BigY - 2*BigY**2) + 5*Pi**2)*ss**2 + 
     -            8*(1 + BigS - 4*BigT + BigU + 2*BigX**2 - 
     -               2*BigX*BigY + BigY**2 + Pi**2)*ss*uu + 
     -            (2*(-6 + 6*BigS - 5*BigT - 3*BigT**2 + 
     -                  6*BigX**2 - 8*BigX*BigY + 4*BigY**2) + 
     -               9*Pi**2)*uu**2)*xn + 
     -         NF*(-3*(BigS + BigU)*
     -             ((-6*(-4 + BigT + BigT**2) + Pi**2)*ss**2 + 
     -               12*(1 - 2*BigT)*ss*uu + 
     -               (-6*(-4 + BigT + BigT**2) + Pi**2)*uu**2) + 
     -            NF*(-9*BigS**2 + BigS*(20 - 6*BigU) + 
     -               (20 - 9*BigU)*BigU + 5*Pi**2)*(ss**2 + uu**2)*
     -             xn - 3*(BigS + BigU)*
     -             ((3*(-4 + 5*BigS + 2*BigS**2 + 5*BigU + 
     -                     2*BigU**2) - 8*Pi**2)*ss**2 + 
     -               12*(1 + 2*BigS + 2*BigU)*ss*uu + 
     -               (3*(-4 + 5*BigS + 2*BigS**2 + 5*BigU + 
     -                     2*BigU**2) - 8*Pi**2)*uu**2)*xn**2) + 
     -         CF*(-3*Pi**4*(ss**2 + uu**2)*(-1 + 8*xn**2) + 
     -            Pi**2*(36*ss*uu*
     -                (1 - 2*BigT - 2*(BigS + BigU)*NF*xn + 
     -                  (1 + 2*BigS + 2*BigU)*xn**2) + 
     -               ss**2*(36*(2 + BigT - BigT**2) - 
     -                  2*(-2 + 21*BigS + 3*BigT + 21*BigU)*NF*xn + 
     -                  9*(-4 + 5*BigS + 2*BigS**2 - 48*BigT + 
     -                     16*BigT**2 + 5*BigU + 2*BigU**2)*xn**2)
     -                + uu**2*
     -                (36*(2 + BigT - BigT**2) - 
     -                  2*(-2 + 39*BigS + 3*BigT + 39*BigU)*NF*xn + 
     -                  9*(-4 + 5*BigS + 2*BigS**2 - 48*BigT + 
     -                     16*BigT**2 + 5*BigU + 2*BigU**2)*xn**2))
     -             - 6*(-18*(-3 + BigT)*BigT*
     -                ((-4 + BigT + BigT**2)*ss**2 + 
     -                  2*(-1 + 2*BigT)*ss*uu + 
     -                  (-4 + BigT + BigT**2)*uu**2) + 
     -               9*(-3 + BigT)*BigT*
     -                ((-4 + 5*BigS + 2*BigS**2 + 5*BigU + 
     -                     2*BigU**2)*ss**2 + 
     -                  4*(1 + 2*BigS + 2*BigU)*ss*uu + 
     -                  (-4 + 5*BigS + 2*BigS**2 + 5*BigU + 
     -                     2*BigU**2)*uu**2)*xn**2 - 
     -               2*NF*xn*
     -                (-6*(BigS + BigU)*
     -                   (1 + BigS - 4*BigT + BigU + 2*BigX**2 - 
     -                     2*BigX*BigY + BigY**2)*ss*uu + 
     -                  ss**2*
     -                   (-2*BigT**3 + 9*BigU - 9*BigU**2 + 
     -                     6*BigT*(-5 + 2*BigU) + 
     -                     BigT**2*(19 + 3*BigU) - 9*BigU*BigX**2 + 
     -                     6*BigU*BigX*BigY - 3*BigU*BigY**2 + 
     -                     3*BigS*
     -                      (3 + 4*BigT + BigT**2 - 3*BigU - 
     -                       3*BigX**2 + 2*BigX*BigY - BigY**2) + 
     -                     3*Zeta3) + 
     -                  uu**2*
     -                   (-9*BigS**2 - 2*BigT**3 + 
     -                     6*BigT*(-5 + 2*BigU) + 
     -                     BigT**2*(19 + 3*BigU) + 
     -                     3*BigS*
     -                      (3 + 4*BigT + BigT**2 - 3*BigU - 
     -                       3*BigX**2 + 4*BigX*BigY - 2*BigY**2) + 
     -                     3*
     -                      (BigU*
     -                       (3 - 3*BigX**2 + 4*BigX*BigY - 
     -                       2*BigY**2) + Zeta3))))))))/(54.*ss*uu)
      endif

      if (crosscheck) then
c--- JC implementation
      if (order >= 1) then
      coeff1fin =  + im*pi * (
     &     - 2.D0/3.D0*TR*nf
     &     - 2.D0/3.D0*TR*nf*yy
     &     - 2.D0/3.D0*TR*nf*xx*yy**(-1)
     &     )
      coeff1fin = coeff1fin + CA * (
     &     - 7.D0
     &     - 7.D0*yy
     &     - 7.D0*xx*yy**(-1)
     &     - 1.D0/3.D0*BigU
     &     - 10.D0/3.D0*BigU*yy
     &     - 1.D0/3.D0*BigU*xx*yy**(-1)
     &     - 10.D0/3.D0*BigS
     &     - 1.D0/3.D0*BigS*yy
     &     - 10.D0/3.D0*BigS*xx*yy**(-1)
     &     + BigY**2
     &     + BigY**2*yy
     &     + BigY**2*xx*yy**(-1)
     &     )
      coeff1fin = coeff1fin + CA*im*pi * (
     &     + 10.D0/3.D0
     &     + 1.D0/3.D0*yy
     &     + 10.D0/3.D0*xx*yy**(-1)
     &     + 2.D0*BigY
     &     + 2.D0*BigY*yy
     &     + 2.D0*BigY*xx*yy**(-1)
     &     )
      coeff1fin = coeff1fin + CA**2 * (
     &     - 2.D0/3.D0*TR*nf*BigU
     &     - 2.D0/3.D0*TR*nf*BigU*yy
     &     - 2.D0/3.D0*TR*nf*BigU*xx*yy**(-1)
     &     - 2.D0/3.D0*TR*nf*BigS
     &     - 2.D0/3.D0*TR*nf*BigS*yy
     &     - 2.D0/3.D0*TR*nf*BigS*xx*yy**(-1)
     &     )
      coeff1fin = coeff1fin + CA**2*im*pi * (
     &     + 2.D0/3.D0*TR*nf
     &     + 2.D0/3.D0*TR*nf*yy
     &     + 2.D0/3.D0*TR*nf*xx*yy**(-1)
     &     )
      coeff1fin = coeff1fin + CA**3 * (
     &     + 7.D0
     &     + 7.D0*yy
     &     + 7.D0*xx*yy**(-1)
     &     + 1.D0/3.D0*BigU
     &     + 10.D0/3.D0*BigU*yy
     &     + 1.D0/3.D0*BigU*xx*yy**(-1)
     &     + 10.D0/3.D0*BigS
     &     + 1.D0/3.D0*BigS*yy
     &     + 10.D0/3.D0*BigS*xx*yy**(-1)
     &     - BigY**2
     &     - BigY**2*yy
     &     - BigY**2*xx*yy**(-1)
     &     )
      coeff1fin = coeff1fin + CA**3*im*pi * (
     &     - 10.D0/3.D0
     &     - 1.D0/3.D0*yy
     &     - 10.D0/3.D0*xx*yy**(-1)
     &     - 2.D0*BigY
     &     - 2.D0*BigY*yy
     &     - 2.D0*BigY*xx*yy**(-1)
     &     )
      coeff1fin = coeff1fin + CF * (
     &     - 14.D0
     &     - 14.D0*yy
     &     - 14.D0*xx*yy**(-1)
     &     + 2.D0*BigU
     &     + 6.D0*BigU*xx*yy**(-1)
     &     + 2.D0*BigT
     &     - 6.D0*BigT*yy
     &     - 6.D0*BigT*xx*yy**(-1)
     &     - 4.D0*BigS
     &     + 6.D0*BigS*yy
     &     - 2.D0*BigY**2
     &     + 4.D0*BigY**2*yy
     &     + 2.D0*BigY**2*xx*yy**(-1)
     &     + 4.D0*BigX*BigY
     &     )
      coeff1fin = coeff1fin + CF * (
     &     - 8.D0*BigX*BigY*yy
     &     - 4.D0*BigX*BigY*xx*yy**(-1)
     &     - 2.D0*BigX**2
     &     + 6.D0*BigX**2*yy
     &     + 6.D0*BigX**2*xx*yy**(-1)
     &     - 2.D0*pisq
     &     + 4.D0*pisq*yy
     &     + 2.D0*pisq*xx*yy**(-1)
     &     )
      coeff1fin = coeff1fin + CF*im*pi * (
     &     + 4.D0
     &     - 6.D0*yy
     &     + 4.D0*BigX*yy
     &     + 8.D0*BigX*xx*yy**(-1)
     &     )
      coeff1fin = coeff1fin + 2.D0/3.D0*TR*nf*BigU
     &     + 2.D0/3.D0*TR*nf*BigU*yy
     &     + 2.D0/3.D0*TR*nf*BigU*xx*yy**(-1)
     &     + 2.D0/3.D0*TR*nf*BigS
     &     + 2.D0/3.D0*TR*nf*BigS*yy
     &     + 2.D0/3.D0*TR*nf*BigS*xx*yy**(-1)
     &

      coeff1C0 =  + im*pi * (
     &     + 1.D0/6.D0*TR*nf
     &     )
      coeff1C0 = coeff1C0 + CA * (
     &     + 5.D0/6.D0*BigU
     &     - 1.D0/4.D0*BigU**2
     &     - 3.D0/4.D0*BigT
     &     + 1.D0/4.D0*BigT**2
     &     + 5.D0/6.D0*BigS
     &     - 1.D0/4.D0*BigS**2
     &     + 7.D0/24.D0*pisq
     &     )
      coeff1C0 = coeff1C0 + CA*im*pi * (
     &     - 5.D0/6.D0
     &     + 1.D0/2.D0*BigS
     &     )
      coeff1C0 = coeff1C0 + CF * (
     &     + 3.D0/2.D0*BigT
     &     - 1.D0/2.D0*BigT**2
     &     + 1.D0/12.D0*pisq
     &     )
      coeff1C0 = coeff1C0 - 1.D0/6.D0*TR*nf*BigU
     &     - 1.D0/6.D0*TR*nf*BigS

! factor of two for interference
      jchard1=two*(coeff1fin+coeff1C0*amps0)
      
      write(6,*) 'qg: hard1,jchard1,hard1/jchard1',hard1,jchard1,hard1/jchard1
!      pause
      
      endif
      
      if (order == 2) then

      Eq13curlybr =  + im*pi * (
     &     + 5.D0/8.D0*BigT
     &     - 5.D0/24.D0*BigT**2
     &     - 3.D0/8.D0*BigS*BigT
     &     + 1.D0/8.D0*BigS*BigT**2
     &     + 5.D0/144.D0*pisq
     &     - 1.D0/48.D0*pisq*BigS
     &     - 5.D0/27.D0*TR**2*nf**2
     &     - 1.D0/36.D0*TR**2*nf**2*BigU
     &     + 1.D0/12.D0*TR**2*nf**2*BigS
     &     )
      Eq13curlybr = Eq13curlybr + CA * (
     &     - 167.D0/108.D0*TR*nf*BigU
     &     + 41.D0/72.D0*TR*nf*BigU**2
     &     - 1.D0/72.D0*TR*nf*BigU**3
     &     + 5.D0/6.D0*TR*nf*BigT
     &     + 1.D0/8.D0*TR*nf*BigT*BigU
     &     - 19.D0/36.D0*TR*nf*BigT**2
     &     - 1.D0/24.D0*TR*nf*BigT**2*BigU
     &     + 1.D0/18.D0*TR*nf*BigT**3
     &     - 167.D0/108.D0*TR*nf*BigS
     &     - 5.D0/18.D0*TR*nf*BigS*BigU
     &     + 1.D0/24.D0*TR*nf*BigS*BigU**2
     &     + 1.D0/8.D0*TR*nf*BigS*BigT
     &     - 1.D0/24.D0*TR*nf*BigS*BigT**2
     &     + 41.D0/72.D0*TR*nf*BigS**2
     &     )
      Eq13curlybr = Eq13curlybr + CA * (
     &     + 1.D0/24.D0*TR*nf*BigS**2*BigU
     &     - 1.D0/72.D0*TR*nf*BigS**3
     &     + 1.D0/12.D0*TR*nf*zeta3
     &     - 229.D0/432.D0*TR*nf*pisq
     &     - 5.D0/144.D0*TR*nf*pisq*BigU
     &     + 1.D0/72.D0*TR*nf*pisq*BigT
     &     + 7.D0/144.D0*TR*nf*pisq*BigS
     &     )
      Eq13curlybr = Eq13curlybr + CA*im*pi * (
     &     + 167.D0/108.D0*TR*nf
     &     + 5.D0/18.D0*TR*nf*BigU
     &     - 1.D0/24.D0*TR*nf*BigU**2
     &     - 1.D0/8.D0*TR*nf*BigT
     &     + 1.D0/24.D0*TR*nf*BigT**2
     &     - 41.D0/36.D0*TR*nf*BigS
     &     - 1.D0/12.D0*TR*nf*BigS*BigU
     &     + 1.D0/24.D0*TR*nf*BigS**2
     &     - 1.D0/48.D0*TR*nf*pisq
     &     )
      Eq13curlybr = Eq13curlybr + CA**2 * (
     &     + 335.D0/108.D0*BigU
     &     - 97.D0/72.D0*BigU**2
     &     - 1.D0/18.D0*BigU**3
     &     + 1.D0/32.D0*BigU**4
     &     + 9.D0/32.D0*BigT**2
     &     - 3.D0/16.D0*BigT**3
     &     + 1.D0/32.D0*BigT**4
     &     + 335.D0/108.D0*BigS
     &     + 25.D0/36.D0*BigS*BigU
     &     - 5.D0/24.D0*BigS*BigU**2
     &     - 97.D0/72.D0*BigS**2
     &     - 5.D0/24.D0*BigS**2*BigU
     &     + 1.D0/16.D0*BigS**2*BigU**2
     &     - 1.D0/18.D0*BigS**3
     &     )
      Eq13curlybr = Eq13curlybr + CA**2 * (
     &     + 1.D0/32.D0*BigS**4
     &     - 11.D0/24.D0*zeta3
     &     + 1121.D0/864.D0*pisq
     &     + 17.D0/96.D0*pisq*BigU
     &     - 1.D0/24.D0*pisq*BigU**2
     &     + 1.D0/32.D0*pisq*BigT
     &     - 1.D0/96.D0*pisq*BigT**2
     &     + 13.D0/96.D0*pisq*BigS
     &     - 1.D0/6.D0*pisq*BigS**2
     &     + 13.D0/1152.D0*pisq**2
     &     )
      Eq13curlybr = Eq13curlybr + CA**2*im*pi * (
     &     - 335.D0/108.D0
     &     - 25.D0/36.D0*BigU
     &     + 5.D0/24.D0*BigU**2
     &     + 97.D0/36.D0*BigS
     &     + 5.D0/12.D0*BigS*BigU
     &     - 1.D0/8.D0*BigS*BigU**2
     &     + 1.D0/6.D0*BigS**2
     &     - 1.D0/8.D0*BigS**3
     &     - 7.D0/288.D0*pisq
     &     + 1.D0/12.D0*pisq*BigS
     &     )
      Eq13curlybr = Eq13curlybr + CF * (
     &     - 5.D0/3.D0*TR*nf*BigT
     &     - 1.D0/4.D0*TR*nf*BigT*BigU
     &     + 19.D0/18.D0*TR*nf*BigT**2
     &     + 1.D0/12.D0*TR*nf*BigT**2*BigU
     &     - 1.D0/9.D0*TR*nf*BigT**3
     &     - 1.D0/4.D0*TR*nf*BigS*BigT
     &     + 1.D0/12.D0*TR*nf*BigS*BigT**2
     &     + 1.D0/6.D0*TR*nf*zeta3
     &     + 1.D0/54.D0*TR*nf*pisq
     &     - 1.D0/72.D0*TR*nf*pisq*BigU
     &     - 1.D0/36.D0*TR*nf*pisq*BigT
     &     - 1.D0/72.D0*TR*nf*pisq*BigS
     &     )
      Eq13curlybr = Eq13curlybr + CF*im*pi * (
     &     + 1.D0/4.D0*TR*nf*BigT
     &     - 1.D0/12.D0*TR*nf*BigT**2
     &     + 1.D0/72.D0*TR*nf*pisq
     &     )
      Eq13curlybr = Eq13curlybr + CF*CA * (
     &     - 9.D0/8.D0*BigT**2
     &     + 3.D0/4.D0*BigT**3
     &     - 1.D0/8.D0*BigT**4
     &     - 1.D0/8.D0*pisq*BigT
     &     + 1.D0/24.D0*pisq*BigT**2
     &     - 1.D0/288.D0*pisq**2
     &     )
      Eq13curlybr = Eq13curlybr + CF**2 * (
     &     + 9.D0/8.D0*BigT**2
     &     - 3.D0/4.D0*BigT**3
     &     + 1.D0/8.D0*BigT**4
     &     + 1.D0/8.D0*pisq*BigT
     &     - 1.D0/24.D0*pisq*BigT**2
     &     + 1.D0/288.D0*pisq**2
     &     )
      Eq13curlybr = Eq13curlybr - 67.D0/24.D0*BigT
     &     - 5.D0/8.D0*BigT*BigU
     &     + 3.D0/16.D0*BigT*BigU**2
     &     + 233.D0/144.D0*BigT**2
     &     + 5.D0/24.D0*BigT**2*BigU
     &     - 1.D0/16.D0*BigT**2*BigU**2
     &     - 11.D0/72.D0*BigT**3
     &     - 5.D0/8.D0*BigS*BigT
     &     + 5.D0/24.D0*BigS*BigT**2
     &     + 3.D0/16.D0*BigS**2*BigT
     &     - 1.D0/16.D0*BigS**2*BigT**2
     &     + 11.D0/48.D0*zeta3
     &     + 1.D0/54.D0*pisq
     &     - 5.D0/144.D0*pisq*BigU
     &     + 1.D0/96.D0*pisq*BigU**2
     &
      Eq13curlybr = Eq13curlybr - 47.D0/288.D0*pisq*BigT
     &     + 1.D0/24.D0*pisq*BigT**2
     &     - 5.D0/144.D0*pisq*BigS
     &     + 1.D0/96.D0*pisq*BigS**2
     &     - 7.D0/576.D0*pisq**2
     &     + 5.D0/27.D0*TR**2*nf**2*BigU
     &     - 1.D0/24.D0*TR**2*nf**2*BigU**2
     &     + 5.D0/27.D0*TR**2*nf**2*BigS
     &     + 1.D0/36.D0*TR**2*nf**2*BigS*BigU
     &     - 1.D0/24.D0*TR**2*nf**2*BigS**2
     &     + 7.D0/216.D0*TR**2*nf**2*pisq

! factor of two for interferences
      jchard2=two*(coeff1C0*coeff1fin+Eq13curlybr*amps0)

! extra pieces, c.f. gagaAGTYassemble1loop
      jchard2=jchard2+two*coeff1C0*conjg(coeff1fin)
     & +coeff1C0*conjg(coeff1C0)*amps0

!      write(6,*) '2-loop bit',two*coeff1C0*coeff1fin+two*Eq13curlybr*amps0
!      write(6,*) '1-loop bit',two*coeff1C0*conjg(coeff1fin)+coeff1C0*conjg(coeff1C0)*amps0
      write(6,*) 'qg: hard2,jchard2,hard2/jchard2',hard2,jchard2,hard2/jchard2
      pause
      
      endif
      endif

      return
      end

