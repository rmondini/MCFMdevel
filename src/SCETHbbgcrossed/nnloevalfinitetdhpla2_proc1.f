
      subroutine nnloevalfinitetdhpla2_proc1(uu,vv,musq,a2finite) 
      implicit none 
      include 'types.f'
      include 'constants.f'
      include 'nf.f' 
      include 'ewcouple.f' 
      include 'masses.f'
      include 'scet_const.f' 
      integer i,i1,i2,i3,i4
      real(dp) xx,yy,zz
      real(dp) a2born
      complex(dp) a2finitelogcoeff(5)

      real(dp),intent(in)::uu,vv,musq
      complex(dp),intent(out)::a2finite

      complex(dp) lnrat

      real(dp) H1(0:1),H2(0:1,0:1),H3(0:1,0:1,0:1),H4(0:1,0:1,0:1,0:1)
      real(dp) G1(0:3),G2(0:3,0:3),G3(0:3,0:3,0:3),G4(0:3,0:3,0:3,0:3)

c---- value for x,y,z
      xx = one/vv
      yy = -uu/vv
      zz = one-xx-yy

c---- form factor at LO with overall factor sqrt(as/(2*pi))
      a2born=two*rt2*pi*(-2/yy)

      call tdhpl(uu,vv,4,G1,G2,G3,G4,H1,H2,H3,H4)

      a2finitelogcoeff(1) = 289._dp/288._dp

      a2finitelogcoeff(2) = (391 + (0,918)*pi + 459*G1(0) + 459*G1(2) - 867*H1(0) - 459*H1(1))/216. 

      a2finitelogcoeff(3) = (-4318 + (0,2484)*pi - 19773*zeta2 + 1242*G1(0) + 
     -    (0,3780)*pi*G1(0) - (0,816)*pi*G1(1) + 1242*G1(2) + 
     -    (0,3780)*pi*G1(2) - (0,816)*pi*G1(3) + 1890*G2(0,0) + 
     -    1890*G2(0,2) - 816*G2(1,2) + 1890*G2(2,0) + 1890*G2(2,2) - 
     -    816*G2(3,0) - 2346*H1(0) - (0,6324)*pi*H1(0) - 
     -    3570*G1(0)*H1(0) + 816*G1(1)*H1(0) - 3570*G1(2)*H1(0) + 
     -    816*G1(3)*H1(0) - 1242*H1(1) - (0,3780)*pi*H1(1) - 
     -    1890*G1(0)*H1(1) + 816*G1(1)*H1(1) - 1890*G1(2)*H1(1) + 
     -    6018*H2(0,0) + 2754*H2(0,1) + 3774*H2(1,0) + 1890*H2(1,1))/
     -  432.  

      a2finitelogcoeff(4) = (12280 - (0,13716)*pi - 40716*zeta2 - (0,25758)*pi*zeta2 - 
     -    11394*zeta3 - 6858*G1(0) + (0,7452)*pi*G1(0) - 
     -    47871*zeta2*G1(0) - (0,3312)*pi*G1(1) + 15552*zeta2*G1(1) - 
     -    6858*G1(2) + (0,7452)*pi*G1(2) - 47871*zeta2*G1(2) - 
     -    (0,3312)*pi*G1(3) + 15552*zeta2*G1(3) + 3726*G2(0,0) + 
     -    (0,8748)*pi*G2(0,0) - (0,1296)*pi*G2(0,1) + 3726*G2(0,2) + 
     -    (0,8748)*pi*G2(0,2) - (0,1296)*pi*G2(0,3) - 
     -    (0,1296)*pi*G2(1,0) - 3312*G2(1,2) - (0,3888)*pi*G2(1,2) + 
     -    3726*G2(2,0) + (0,8748)*pi*G2(2,0) - (0,1296)*pi*G2(2,1) + 
     -    3726*G2(2,2) + (0,8748)*pi*G2(2,2) - (0,1296)*pi*G2(2,3) - 
     -    3312*G2(3,0) - (0,3888)*pi*G2(3,0) - (0,1296)*pi*G2(3,2) + 
     -    4374*G3(0,0,0) + 4374*G3(0,0,2) - 1296*G3(0,1,2) + 
     -    4374*G3(0,2,0) + 4374*G3(0,2,2) - 1296*G3(0,3,0) - 
     -    1296*G3(1,0,2) - 1296*G3(1,2,0) - 2592*G3(1,2,2) + 
     -    4374*G3(2,0,0) + 4374*G3(2,0,2) - 1296*G3(2,1,2) + 
     -    4374*G3(2,2,0) + 4374*G3(2,2,2) - 1296*G3(2,3,0) - 
     -    2592*G3(3,0,0) - 1296*G3(3,0,2) - 1296*G3(3,2,0) + 
     -    12954*H1(0) - (0,10764)*pi*H1(0) + 74871*zeta2*H1(0) - 
     -    7038*G1(0)*H1(0) - (0,15228)*pi*G1(0)*H1(0) + 
     -    3312*G1(1)*H1(0) + (0,5040)*pi*G1(1)*H1(0) - 
     -    7038*G1(2)*H1(0) - (0,15228)*pi*G1(2)*H1(0) + 
     -    3312*G1(3)*H1(0) + (0,5040)*pi*G1(3)*H1(0) - 
     -    8262*G2(0,0)*H1(0) + 1296*G2(0,1)*H1(0) - 
     -    8262*G2(0,2)*H1(0) + 1296*G2(0,3)*H1(0) + 
     -    1296*G2(1,0)*H1(0) + 3744*G2(1,2)*H1(0) - 
     -    8262*G2(2,0)*H1(0) + 1296*G2(2,1)*H1(0) - 
     -    8262*G2(2,2)*H1(0) + 1296*G2(2,3)*H1(0) + 
     -    3744*G2(3,0)*H1(0) + 1296*G2(3,2)*H1(0) + 6858*H1(1) - 
     -    (0,7452)*pi*H1(1) + 47871*zeta2*H1(1) - 3726*G1(0)*H1(1) - 
     -    (0,8748)*pi*G1(0)*H1(1) + 3312*G1(1)*H1(1) + 
     -    (0,3888)*pi*G1(1)*H1(1) - 3726*G1(2)*H1(1) - 
     -    (0,8748)*pi*G1(2)*H1(1) + (0,1296)*pi*G1(3)*H1(1) - 
     -    4374*G2(0,0)*H1(1) + 1296*G2(0,1)*H1(1) - 
     -    4374*G2(0,2)*H1(1) + 1296*G2(1,0)*H1(1) + 
     -    2592*G2(1,2)*H1(1) - 4374*G2(2,0)*H1(1) + 
     -    1296*G2(2,1)*H1(1) - 4374*G2(2,2)*H1(1) + 
     -    1296*G2(3,0)*H1(1) + 10350*H2(0,0) + (0,24012)*pi*H2(0,0) + 
     -    14454*G1(0)*H2(0,0) - 4896*G1(1)*H2(0,0) + 
     -    14454*G1(2)*H2(0,0) - 4896*G1(3)*H2(0,0) + 3726*H2(0,1) + 
     -    (0,12636)*pi*H2(0,1) + 6966*G1(0)*H2(0,1) - 
     -    3744*G1(1)*H2(0,1) + 6966*G1(2)*H2(0,1) - 
     -    1296*G1(3)*H2(0,1) + 7866*H2(1,0) + (0,15876)*pi*H2(1,0) + 
     -    8586*G1(0)*H2(1,0) - 3744*G1(1)*H2(1,0) + 
     -    8586*G1(2)*H2(1,0) - 1296*G1(3)*H2(1,0) + 3726*H2(1,1) + 
     -    (0,8748)*pi*H2(1,1) + 4374*G1(0)*H2(1,1) - 
     -    2592*G1(1)*H2(1,1) + 4374*G1(2)*H2(1,1) - 22950*H3(0,0,0) - 
     -    9558*H3(0,0,1) - 12618*H3(0,1,0) - 5670*H3(0,1,1) - 
     -    15678*H3(1,0,0) - 7290*H3(1,0,1) - 8910*H3(1,1,0) - 
     -    4374*H3(1,1,1))/648.  

      a2finitelogcoeff(5) = (137573*pi**4)/51840. + 
     -  ((0,6)*pi*vv + (0,6)*pi*vv**2 + zeta2/3. - (vv*zeta2)/3. + 
     -     (0,0.3333333333333333)*pi*H1(0) - 6*vv*H1(0) - 
     -     (0,0.3333333333333333)*pi*vv*H1(0) - 6*vv**2*H1(0))/
     -   ((1 - uu - vv)*(uu + vv)) + 
     -  ((0,-6)*pi*vv**2 - (0,6)*pi*vv**3 + 6*vv**2*H1(0) + 
     -     6*vv**3*H1(0))/((1 - uu)*(1 - uu - vv)*(uu + vv)) + 
     -  (-zeta2/3. - (17*vv*zeta2)/9. - 
     -     (0,4.444444444444445)*pi*vv*zeta2 + (80*vv**2*zeta2)/9. - 
     -     (20*vv*zeta3)/9. + 6*vv*G1(0) - 
     -     (0,2.2222222222222223)*pi*vv**2*G1(0) - 
     -     (20*vv*zeta2*G1(0))/9. - (0,0.3333333333333333)*pi*H1(0) - 
     -     (0,1.8888888888888888)*pi*vv*H1(0) + 
     -     (0,2.2222222222222223)*pi*vv**2*H1(0) + 
     -     (80*vv*zeta2*H1(0))/9.)/(uu + vv) + 
     -  ((20*vv*zeta2)/9. + (0,4.444444444444445)*pi*vv*zeta2 - 
     -     (100*vv**2*zeta2)/9. - (0,8.88888888888889)*pi*vv**2*zeta2 + 
     -     (80*vv**3*zeta2)/9. + (20*vv*zeta3)/9. - 
     -     (40*vv**2*zeta3)/9. + 
     -     (0,2.2222222222222223)*pi*vv**2*G1(0) - 
     -     (0,2.2222222222222223)*pi*vv**3*G1(0) + 
     -     (20*vv*zeta2*G1(0))/9. - (40*vv**2*zeta2*G1(0))/9. + 
     -     (0,2.2222222222222223)*pi*vv*H1(0) - 
     -     (0,4.444444444444445)*pi*vv**2*H1(0) + 
     -     (0,2.2222222222222223)*pi*vv**3*H1(0) - 
     -     (80*vv*zeta2*H1(0))/9. + (160*vv**2*zeta2*H1(0))/9.)/
     -   ((1 - vv)*(uu + vv)) + 
     -  ((0,4.444444444444445)*pi*vv**2*zeta2 - 
     -     (0,4.444444444444445)*pi*vv**3*zeta2 + (20*vv**2*zeta3)/9. - 
     -     (20*vv**3*zeta3)/9. + (20*vv**2*zeta2*G1(0))/9. - 
     -     (20*vv**3*zeta2*G1(0))/9. - (80*vv**2*zeta2*H1(0))/9. + 
     -     (80*vv**3*zeta2*H1(0))/9.)/((1 - vv)**2*(uu + vv)) + 
     -  (6*vv*G1(2) - 6*vv*H1(1))/(1 - uu) + 
     -  ((0,0.3333333333333333)*pi*G1(1) - 
     -     (0,0.3333333333333333)*pi*vv*G1(1) + G2(1,2)/3. - 
     -     (vv*G2(1,2))/3. - (G1(1)*H1(0))/3. + (vv*G1(1)*H1(0))/3. + 
     -     (G1(2)*H1(0))/3. - (vv*G1(2)*H1(0))/3. - (G1(1)*H1(1))/3. + 
     -     (vv*G1(1)*H1(1))/3.)/uu + 
     -  ((0,0.3333333333333333)*pi*G1(3) - 
     -     (0,0.3333333333333333)*pi*vv*G1(3) + G2(3,0)/3. - 
     -     (vv*G2(3,0))/3. + (G1(0)*H1(0))/3. - (vv*G1(0)*H1(0))/3. - 
     -     (G1(3)*H1(0))/3. + (vv*G1(3)*H1(0))/3. - H2(0,0)/3. + 
     -     (vv*H2(0,0))/3. + H2(1,0)/3. - (vv*H2(1,0))/3.)/(1 - uu - vv)
     -    + ((-20*zeta2)/9. - (0,4.444444444444445)*pi*zeta2 - 
     -     (100*uu*zeta2)/9. - (0,8.88888888888889)*pi*uu*zeta2 - 
     -     (80*uu**2*zeta2)/9. + (100*vv*zeta2)/9. + 
     -     (0,8.88888888888889)*pi*vv*zeta2 + (80*uu*vv*zeta2)/9. - 
     -     (80*vv**2*zeta2)/9. - (20*zeta3)/9. - (40*uu*zeta3)/9. + 
     -     (40*vv*zeta3)/9. + (0,2.2222222222222223)*pi*uu*G1(0) + 
     -     (0,2.2222222222222223)*pi*uu**2*G1(0) - 
     -     (0,2.2222222222222223)*pi*vv*G1(0) - 
     -     (0,2.2222222222222223)*pi*uu*vv*G1(0) + 
     -     (0,2.2222222222222223)*pi*vv**2*G1(0) - 
     -     (20*zeta2*G1(0))/9. - (40*uu*zeta2*G1(0))/9. + 
     -     (40*vv*zeta2*G1(0))/9. + 
     -     (0,2.2222222222222223)*pi*uu*G1(1) - 
     -     (0,2.2222222222222223)*pi*uu**2*G1(1) + 
     -     (0,2.2222222222222223)*pi*uu**2*G1(2) + 
     -     (20*uu*zeta2*G1(2))/9. - 
     -     (0,2.2222222222222223)*pi*uu*G1(3) - 
     -     (0,2.2222222222222223)*pi*uu**2*G1(3) + 
     -     (0,2.2222222222222223)*pi*uu*G2(0,1) + 
     -     (20*uu**2*G2(0,2))/9. - 
     -     (0,2.2222222222222223)*pi*uu*G2(0,3) + (20*uu*G2(1,2))/9. - 
     -     (20*uu**2*G2(1,2))/9. + (20*uu**2*G2(2,0))/9. - 
     -     (0,2.2222222222222223)*pi*uu*G2(2,1) + 
     -     (0,2.2222222222222223)*pi*uu*G2(2,3) - (20*uu*G2(3,0))/9. - 
     -     (20*uu**2*G2(3,0))/9. + (20*uu*G3(0,1,2))/9. - 
     -     (20*uu*G3(0,3,0))/9. - (20*uu*G3(2,1,2))/9. + 
     -     (20*uu*G3(2,3,0))/9. - (0,2.2222222222222223)*pi*H1(0) - 
     -     (0,4.444444444444445)*pi*uu*H1(0) - 
     -     (0,2.2222222222222223)*pi*uu**2*H1(0) + 
     -     (0,4.444444444444445)*pi*vv*H1(0) + 
     -     (0,2.2222222222222223)*pi*uu*vv*H1(0) - 
     -     (0,2.2222222222222223)*pi*vv**2*H1(0) + 
     -     (80*zeta2*H1(0))/9. + (160*uu*zeta2*H1(0))/9. - 
     -     (160*vv*zeta2*H1(0))/9. - 
     -     (0,2.2222222222222223)*pi*uu*G1(0)*H1(0) - 
     -     (20*uu**2*G1(0)*H1(0))/9. - (20*uu*G1(1)*H1(0))/9. + 
     -     (20*uu**2*G1(1)*H1(0))/9. + 
     -     (0,2.2222222222222223)*pi*uu*G1(2)*H1(0) - 
     -     (20*uu**2*G1(2)*H1(0))/9. + (20*uu*G1(3)*H1(0))/9. + 
     -     (20*uu**2*G1(3)*H1(0))/9. - (20*uu*G2(0,1)*H1(0))/9. + 
     -     (20*uu*G2(0,3)*H1(0))/9. + (20*uu*G2(2,1)*H1(0))/9. - 
     -     (20*uu*G2(2,3)*H1(0))/9. - 
     -     (0,2.2222222222222223)*pi*uu**2*H1(1) - 
     -     (20*uu*zeta2*H1(1))/9. - (20*uu**2*G1(0)*H1(1))/9. - 
     -     (20*uu*G1(1)*H1(1))/9. + (20*uu**2*G1(1)*H1(1))/9. - 
     -     (20*uu*G2(0,1)*H1(1))/9. + (20*uu*G2(2,1)*H1(1))/9. + 
     -     (20*uu*H2(0,0))/9. + (0,2.2222222222222223)*pi*uu*H2(0,0) + 
     -     (20*uu**2*H2(0,0))/9. + (20*uu*G1(0)*H2(0,0))/9. - 
     -     (20*uu*G1(2)*H2(0,0))/9. + (20*uu*H2(0,1))/9. + 
     -     (0,2.2222222222222223)*pi*uu*H2(0,1) + 
     -     (20*uu*G1(0)*H2(0,1))/9. - (20*uu*G1(2)*H2(0,1))/9. - 
     -     (0,2.2222222222222223)*pi*uu*H2(1,0) + 
     -     (20*uu**2*H2(1,0))/9. - (20*uu*H3(0,0,0))/9. - 
     -     (20*uu*H3(0,1,0))/9. + (20*uu*H3(1,0,0))/9. + 
     -     (20*uu*H3(1,0,1))/9.)/(1 - vv) + 
     -  ((0,4.444444444444445)*pi*uu*zeta2 + 
     -     (0,4.444444444444445)*pi*uu**2*zeta2 - 
     -     (0,4.444444444444445)*pi*vv*zeta2 - 
     -     (0,4.444444444444445)*pi*uu*vv*zeta2 + 
     -     (0,4.444444444444445)*pi*vv**2*zeta2 + (20*uu*zeta3)/9. + 
     -     (20*uu**2*zeta3)/9. - (20*vv*zeta3)/9. - 
     -     (20*uu*vv*zeta3)/9. + (20*vv**2*zeta3)/9. + 
     -     (20*uu*zeta2*G1(0))/9. + (20*uu**2*zeta2*G1(0))/9. - 
     -     (20*vv*zeta2*G1(0))/9. - (20*uu*vv*zeta2*G1(0))/9. + 
     -     (20*vv**2*zeta2*G1(0))/9. - (20*uu**2*zeta2*G1(2))/9. - 
     -     (0,2.2222222222222223)*pi*uu**2*G2(0,1) + 
     -     (0,2.2222222222222223)*pi*uu**2*G2(0,3) + 
     -     (0,2.2222222222222223)*pi*uu**2*G2(2,1) - 
     -     (0,2.2222222222222223)*pi*uu**2*G2(2,3) - 
     -     (20*uu**2*G3(0,1,2))/9. + (20*uu**2*G3(0,3,0))/9. + 
     -     (20*uu**2*G3(2,1,2))/9. - (20*uu**2*G3(2,3,0))/9. - 
     -     (80*uu*zeta2*H1(0))/9. - (80*uu**2*zeta2*H1(0))/9. + 
     -     (80*vv*zeta2*H1(0))/9. + (80*uu*vv*zeta2*H1(0))/9. - 
     -     (80*vv**2*zeta2*H1(0))/9. + 
     -     (0,2.2222222222222223)*pi*uu**2*G1(0)*H1(0) - 
     -     (0,2.2222222222222223)*pi*uu**2*G1(2)*H1(0) + 
     -     (20*uu**2*G2(0,1)*H1(0))/9. - (20*uu**2*G2(0,3)*H1(0))/9. - 
     -     (20*uu**2*G2(2,1)*H1(0))/9. + (20*uu**2*G2(2,3)*H1(0))/9. + 
     -     (20*uu**2*zeta2*H1(1))/9. + (20*uu**2*G2(0,1)*H1(1))/9. - 
     -     (20*uu**2*G2(2,1)*H1(1))/9. - 
     -     (0,2.2222222222222223)*pi*uu**2*H2(0,0) - 
     -     (20*uu**2*G1(0)*H2(0,0))/9. + (20*uu**2*G1(2)*H2(0,0))/9. - 
     -     (0,2.2222222222222223)*pi*uu**2*H2(0,1) - 
     -     (20*uu**2*G1(0)*H2(0,1))/9. + (20*uu**2*G1(2)*H2(0,1))/9. + 
     -     (0,2.2222222222222223)*pi*uu**2*H2(1,0) + 
     -     (20*uu**2*H3(0,0,0))/9. + (20*uu**2*H3(0,1,0))/9. - 
     -     (20*uu**2*H3(1,0,0))/9. - (20*uu**2*H3(1,0,1))/9.)/
     -   (1 - vv)**2 + (-51866 + (0,78948)*pi - (0,21060)*pi*vv + 
     -     456921*zeta2 - (0,111132)*pi*zeta2 + 162000*uu*zeta2 + 
     -     (0,17280)*pi*uu*zeta2 - 84240*uu**2*zeta2 - 
     -     (0,17280)*pi*uu**2*zeta2 - 34560*vv*zeta2 - 
     -     (0,17928)*pi*vv*zeta2 - 84240*uu*vv*zeta2 - 
     -     (0,17280)*pi*uu*vv*zeta2 - (0,8640)*pi*vv**2*zeta2 + 
     -     67806*zeta3 - (0,54216)*pi*zeta3 + 17928*uu*zeta3 + 
     -     8640*uu**2*zeta3 + 17928*vv*zeta3 + 17280*uu*vv*zeta3 + 
     -     8640*vv**2*zeta3 + 28944*G1(0) - (0,82296)*pi*G1(0) + 
     -     21060*uu*G1(0) - (0,38340)*pi*uu*G1(0) + 
     -     (0,21060)*pi*uu**2*G1(0) + (0,8640)*pi*vv*G1(0) + 
     -     (0,21060)*pi*uu*vv*G1(0) - 234792*zeta2*G1(0) - 
     -     (0,154548)*pi*zeta2*G1(0) + 35208*uu*zeta2*G1(0) - 
     -     8640*uu**2*zeta2*G1(0) - 8748*zeta3*G1(0) + 
     -     (0,35280)*pi*G1(1) + (0,21060)*pi*uu*G1(1) - 
     -     (0,21060)*pi*uu**2*G1(1) - (0,4320)*pi*vv*G1(1) - 
     -     (0,21060)*pi*uu*vv*G1(1) + 31320*zeta2*G1(1) + 
     -     (0,66528)*pi*zeta2*G1(1) + 35208*uu*zeta2*G1(1) - 
     -     8640*uu**2*zeta2*G1(1) - 30024*zeta3*G1(1) + 50004*G1(2) - 
     -     (0,82296)*pi*G1(2) - 21060*uu*G1(2) - 
     -     (0,29700)*pi*uu*G1(2) + (0,21060)*pi*uu**2*G1(2) - 
     -     21060*vv*G1(2) + (0,21060)*pi*uu*vv*G1(2) - 
     -     247644*zeta2*G1(2) - (0,154548)*pi*zeta2*G1(2) + 
     -     24516*zeta3*G1(2) + (0,35280)*pi*G1(3) + 
     -     (0,38340)*pi*uu*G1(3) - (0,21060)*pi*uu**2*G1(3) + 
     -     (0,4320)*pi*vv*G1(3) - (0,21060)*pi*uu*vv*G1(3) + 
     -     58752*zeta2*G1(3) + (0,68688)*pi*zeta2*G1(3) - 
     -     35856*uu*zeta2*G1(3) - 17280*uu**2*zeta2*G1(3) - 
     -     35856*vv*zeta2*G1(3) - 34560*uu*vv*zeta2*G1(3) - 
     -     17280*vv**2*zeta2*G1(3) - 17280*zeta3*G1(3) - 
     -     41148*G2(0,0) + (0,44712)*pi*G2(0,0) - 
     -     287226*zeta2*G2(0,0) - (0,24948)*pi*G2(0,1) - 
     -     (0,17604)*pi*uu*G2(0,1) + (0,4320)*pi*uu**2*G2(0,1) + 
     -     34560*zeta2*G2(0,1) - 41148*G2(0,2) + (0,44712)*pi*G2(0,2) - 
     -     29700*uu*G2(0,2) + 21060*uu**2*G2(0,2) + 
     -     21060*uu*vv*G2(0,2) - 285066*zeta2*G2(0,2) - 
     -     (0,20736)*pi*G2(0,3) + (0,17604)*pi*uu*G2(0,3) - 
     -     (0,4320)*pi*uu**2*G2(0,3) + 79488*zeta2*G2(0,3) + 
     -     (0,13284)*pi*G2(1,0) - (0,17604)*pi*uu*G2(1,0) + 
     -     (0,4320)*pi*uu**2*G2(1,0) + 73872*zeta2*G2(1,0) - 
     -     (0,25704)*pi*G2(1,1) + (0,35208)*pi*uu*G2(1,1) - 
     -     (0,8640)*pi*uu**2*G2(1,1) - 13824*zeta2*G2(1,1) + 
     -     35280*G2(1,2) - (0,6588)*pi*G2(1,2) + 21060*uu*G2(1,2) - 
     -     (0,17604)*pi*uu*G2(1,2) - 21060*uu**2*G2(1,2) + 
     -     (0,4320)*pi*uu**2*G2(1,2) - 4320*vv*G2(1,2) - 
     -     21060*uu*vv*G2(1,2) + 149256*zeta2*G2(1,2) - 41148*G2(2,0) + 
     -     (0,44712)*pi*G2(2,0) - 29700*uu*G2(2,0) + 
     -     21060*uu**2*G2(2,0) + 21060*uu*vv*G2(2,0) - 
     -     267786*zeta2*G2(2,0) - (0,7452)*pi*G2(2,1) - 
     -     (0,8964)*pi*uu*G2(2,1) - (0,4320)*pi*uu**2*G2(2,1) - 
     -     (0,8964)*pi*vv*G2(2,1) - (0,8640)*pi*uu*vv*G2(2,1) - 
     -     (0,4320)*pi*vv**2*G2(2,1) + 65664*zeta2*G2(2,1) - 
     -     41148*G2(2,2) + (0,44712)*pi*G2(2,2) - 
     -     273402*zeta2*G2(2,2) - (0,38232)*pi*G2(2,3) + 
     -     (0,8964)*pi*uu*G2(2,3) + (0,4320)*pi*uu**2*G2(2,3) + 
     -     (0,8964)*pi*vv*G2(2,3) + (0,8640)*pi*uu*vv*G2(2,3) + 
     -     (0,4320)*pi*vv**2*G2(2,3) + 17280*zeta2*G2(2,3) + 
     -     35280*G2(3,0) - (0,19872)*pi*G2(3,0) + 38340*uu*G2(3,0) + 
     -     (0,8964)*pi*uu*G2(3,0) - 21060*uu**2*G2(3,0) + 
     -     (0,4320)*pi*uu**2*G2(3,0) + 4320*vv*G2(3,0) + 
     -     (0,8964)*pi*vv*G2(3,0) - 21060*uu*vv*G2(3,0) + 
     -     (0,8640)*pi*uu*vv*G2(3,0) + (0,4320)*pi*vv**2*G2(3,0) + 
     -     162000*zeta2*G2(3,0) + (0,8964)*pi*uu*G2(3,2) + 
     -     (0,4320)*pi*uu**2*G2(3,2) + (0,8964)*pi*vv*G2(3,2) + 
     -     (0,8640)*pi*uu*vv*G2(3,2) + (0,4320)*pi*vv**2*G2(3,2) + 
     -     72792*zeta2*G2(3,2) + (0,864)*pi*G2(3,3) - 
     -     (0,17928)*pi*uu*G2(3,3) - (0,8640)*pi*uu**2*G2(3,3) - 
     -     (0,17928)*pi*vv*G2(3,3) - (0,17280)*pi*uu*vv*G2(3,3) - 
     -     (0,8640)*pi*vv**2*G2(3,3) - 27648*zeta2*G2(3,3) + 
     -     22356*G3(0,0,0) + (0,52488)*pi*G3(0,0,0) - 
     -     (0,9504)*pi*G3(0,0,1) + 22356*G3(0,0,2) + 
     -     (0,52488)*pi*G3(0,0,2) - (0,23328)*pi*G3(0,0,3) + 
     -     (0,7776)*pi*G3(0,1,0) - (0,17280)*pi*G3(0,1,1) - 
     -     24948*G3(0,1,2) - (0,7776)*pi*G3(0,1,2) - 
     -     17604*uu*G3(0,1,2) + 4320*uu**2*G3(0,1,2) + 
     -     22356*G3(0,2,0) + (0,52488)*pi*G3(0,2,0) - 
     -     (0,9936)*pi*G3(0,2,1) + 22356*G3(0,2,2) + 
     -     (0,52488)*pi*G3(0,2,2) - (0,27216)*pi*G3(0,2,3) - 
     -     20736*G3(0,3,0) - (0,23328)*pi*G3(0,3,0) + 
     -     17604*uu*G3(0,3,0) - 4320*uu**2*G3(0,3,0) - 
     -     (0,7776)*pi*G3(0,3,2) + (0,13824)*pi*G3(0,3,3) - 
     -     (0,7776)*pi*G3(1,0,0) + (0,8640)*pi*G3(1,0,1) + 
     -     13284*G3(1,0,2) - (0,23328)*pi*G3(1,0,2) - 
     -     17604*uu*G3(1,0,2) + 4320*uu**2*G3(1,0,2) + 
     -     (0,9720)*pi*G3(1,0,3) - (0,13824)*pi*G3(1,1,1) - 
     -     25704*G3(1,1,2) + (0,13824)*pi*G3(1,1,2) + 
     -     35208*uu*G3(1,1,2) - 8640*uu**2*G3(1,1,2) + 
     -     13284*G3(1,2,0) - (0,23328)*pi*G3(1,2,0) - 
     -     17604*uu*G3(1,2,0) + 4320*uu**2*G3(1,2,0) + 
     -     (0,22464)*pi*G3(1,2,1) - 19872*G3(1,2,2) - 
     -     (0,38880)*pi*G3(1,2,2) + (0,9720)*pi*G3(1,2,3) + 
     -     22356*G3(2,0,0) + (0,52488)*pi*G3(2,0,0) - 
     -     (0,27216)*pi*G3(2,0,1) + 22356*G3(2,0,2) + 
     -     (0,52488)*pi*G3(2,0,2) - (0,9936)*pi*G3(2,0,3) - 
     -     (0,7776)*pi*G3(2,1,0) + (0,13824)*pi*G3(2,1,1) - 
     -     7452*G3(2,1,2) - (0,23328)*pi*G3(2,1,2) - 
     -     8964*uu*G3(2,1,2) - 4320*uu**2*G3(2,1,2) - 
     -     8964*vv*G3(2,1,2) - 8640*uu*vv*G3(2,1,2) - 
     -     4320*vv**2*G3(2,1,2) + 22356*G3(2,2,0) + 
     -     (0,52488)*pi*G3(2,2,0) - (0,23328)*pi*G3(2,2,1) + 
     -     22356*G3(2,2,2) + (0,52488)*pi*G3(2,2,2) - 
     -     (0,9504)*pi*G3(2,2,3) - 38232*G3(2,3,0) - 
     -     (0,7776)*pi*G3(2,3,0) + 8964*uu*G3(2,3,0) + 
     -     4320*uu**2*G3(2,3,0) + 8964*vv*G3(2,3,0) + 
     -     8640*uu*vv*G3(2,3,0) + 4320*vv**2*G3(2,3,0) + 
     -     (0,7776)*pi*G3(2,3,2) - (0,17280)*pi*G3(2,3,3) - 
     -     19872*G3(3,0,0) - (0,38880)*pi*G3(3,0,0) + 
     -     (0,9720)*pi*G3(3,0,1) - (0,23328)*pi*G3(3,0,2) + 
     -     8964*uu*G3(3,0,2) + 4320*uu**2*G3(3,0,2) + 
     -     8964*vv*G3(3,0,2) + 8640*uu*vv*G3(3,0,2) + 
     -     4320*vv**2*G3(3,0,2) + (0,22464)*pi*G3(3,0,3) - 
     -     (0,23328)*pi*G3(3,2,0) + 8964*uu*G3(3,2,0) + 
     -     4320*uu**2*G3(3,2,0) + 8964*vv*G3(3,2,0) + 
     -     8640*uu*vv*G3(3,2,0) + 4320*vv**2*G3(3,2,0) + 
     -     (0,9720)*pi*G3(3,2,1) - (0,7776)*pi*G3(3,2,2) + 
     -     (0,8640)*pi*G3(3,2,3) + 864*G3(3,3,0) + 
     -     (0,13824)*pi*G3(3,3,0) - 17928*uu*G3(3,3,0) - 
     -     8640*uu**2*G3(3,3,0) - 17928*vv*G3(3,3,0) - 
     -     17280*uu*vv*G3(3,3,0) - 8640*vv**2*G3(3,3,0) - 
     -     (0,13824)*pi*G3(3,3,3) + 26244*G4(0,0,0,0) + 
     -     26244*G4(0,0,0,2) - 9504*G4(0,0,1,2) + 26244*G4(0,0,2,0) + 
     -     26244*G4(0,0,2,2) - 23328*G4(0,0,3,0) + 7776*G4(0,1,0,2) - 
     -     17280*G4(0,1,1,2) + 7776*G4(0,1,2,0) - 15552*G4(0,1,2,2) + 
     -     26244*G4(0,2,0,0) + 26244*G4(0,2,0,2) - 9936*G4(0,2,1,2) + 
     -     26244*G4(0,2,2,0) + 26244*G4(0,2,2,2) - 27216*G4(0,2,3,0) - 
     -     15552*G4(0,3,0,0) - 7776*G4(0,3,0,2) - 7776*G4(0,3,2,0) + 
     -     13824*G4(0,3,3,0) - 7776*G4(1,0,0,2) + 8640*G4(1,0,1,2) - 
     -     7776*G4(1,0,2,0) - 15552*G4(1,0,2,2) + 9720*G4(1,0,3,0) - 
     -     13824*G4(1,1,1,2) + 13824*G4(1,1,2,2) - 7776*G4(1,2,0,0) - 
     -     15552*G4(1,2,0,2) + 22464*G4(1,2,1,2) - 15552*G4(1,2,2,0) - 
     -     23328*G4(1,2,2,2) + 9720*G4(1,2,3,0) + 26244*G4(2,0,0,0) + 
     -     26244*G4(2,0,0,2) - 27216*G4(2,0,1,2) + 26244*G4(2,0,2,0) + 
     -     26244*G4(2,0,2,2) - 9936*G4(2,0,3,0) - 7776*G4(2,1,0,2) + 
     -     13824*G4(2,1,1,2) - 7776*G4(2,1,2,0) - 15552*G4(2,1,2,2) + 
     -     26244*G4(2,2,0,0) + 26244*G4(2,2,0,2) - 23328*G4(2,2,1,2) + 
     -     26244*G4(2,2,2,0) + 26244*G4(2,2,2,2) - 9504*G4(2,2,3,0) - 
     -     15552*G4(2,3,0,0) + 7776*G4(2,3,0,2) + 7776*G4(2,3,2,0) - 
     -     17280*G4(2,3,3,0) - 23328*G4(3,0,0,0) - 15552*G4(3,0,0,2) + 
     -     9720*G4(3,0,1,2) - 15552*G4(3,0,2,0) - 7776*G4(3,0,2,2) + 
     -     22464*G4(3,0,3,0) - 15552*G4(3,2,0,0) - 7776*G4(3,2,0,2) + 
     -     9720*G4(3,2,1,2) - 7776*G4(3,2,2,0) + 8640*G4(3,2,3,0) + 
     -     13824*G4(3,3,0,0) - 13824*G4(3,3,3,0) - 73680*H1(0) + 
     -     (0,126216)*pi*H1(0) + (0,46980)*pi*uu*H1(0) - 
     -     (0,21060)*pi*uu**2*H1(0) + 21060*vv*H1(0) - 
     -     (0,4320)*pi*vv*H1(0) - (0,21060)*pi*uu*vv*H1(0) + 
     -     268488*zeta2*H1(0) + (0,223236)*pi*zeta2*H1(0) - 
     -     70416*uu*zeta2*H1(0) + 17280*uu**2*zeta2*H1(0) + 
     -     51084*zeta3*H1(0) + 76428*G1(0)*H1(0) - 
     -     (0,64584)*pi*G1(0)*H1(0) + 29700*uu*G1(0)*H1(0) + 
     -     (0,17604)*pi*uu*G1(0)*H1(0) - 21060*uu**2*G1(0)*H1(0) - 
     -     (0,4320)*pi*uu**2*G1(0)*H1(0) - 21060*uu*vv*G1(0)*H1(0) + 
     -     449226*zeta2*G1(0)*H1(0) - 35280*G1(1)*H1(0) + 
     -     (0,19008)*pi*G1(1)*H1(0) - 21060*uu*G1(1)*H1(0) + 
     -     21060*uu**2*G1(1)*H1(0) + 4320*vv*G1(1)*H1(0) + 
     -     21060*uu*vv*G1(1)*H1(0) - 176904*zeta2*G1(1)*H1(0) + 
     -     76428*G1(2)*H1(0) - (0,82080)*pi*G1(2)*H1(0) + 
     -     29700*uu*G1(2)*H1(0) + (0,8964)*pi*uu*G1(2)*H1(0) - 
     -     21060*uu**2*G1(2)*H1(0) + (0,4320)*pi*uu**2*G1(2)*H1(0) + 
     -     (0,8964)*pi*vv*G1(2)*H1(0) - 21060*uu*vv*G1(2)*H1(0) + 
     -     (0,8640)*pi*uu*vv*G1(2)*H1(0) + 
     -     (0,4320)*pi*vv**2*G1(2)*H1(0) + 388098*zeta2*G1(2)*H1(0) - 
     -     35280*G1(3)*H1(0) + (0,19872)*pi*G1(3)*H1(0) - 
     -     38340*uu*G1(3)*H1(0) - (0,17928)*pi*uu*G1(3)*H1(0) + 
     -     21060*uu**2*G1(3)*H1(0) - (0,8640)*pi*uu**2*G1(3)*H1(0) - 
     -     4320*vv*G1(3)*H1(0) - (0,17928)*pi*vv*G1(3)*H1(0) + 
     -     21060*uu*vv*G1(3)*H1(0) - (0,17280)*pi*uu*vv*G1(3)*H1(0) - 
     -     (0,8640)*pi*vv**2*G1(3)*H1(0) - 203472*zeta2*G1(3)*H1(0) - 
     -     42228*G2(0,0)*H1(0) - (0,91368)*pi*G2(0,0)*H1(0) + 
     -     24948*G2(0,1)*H1(0) + (0,15336)*pi*G2(0,1)*H1(0) + 
     -     17604*uu*G2(0,1)*H1(0) - 4320*uu**2*G2(0,1)*H1(0) - 
     -     47304*G2(0,2)*H1(0) - (0,93096)*pi*G2(0,2)*H1(0) + 
     -     20736*G2(0,3)*H1(0) + (0,30240)*pi*G2(0,3)*H1(0) - 
     -     17604*uu*G2(0,3)*H1(0) + 4320*uu**2*G2(0,3)*H1(0) - 
     -     13284*G2(1,0)*H1(0) + (0,23328)*pi*G2(1,0)*H1(0) + 
     -     17604*uu*G2(1,0)*H1(0) - 4320*uu**2*G2(1,0)*H1(0) + 
     -     25704*G2(1,1)*H1(0) - 35208*uu*G2(1,1)*H1(0) + 
     -     8640*uu**2*G2(1,1)*H1(0) + 6588*G2(1,2)*H1(0) + 
     -     (0,39960)*pi*G2(1,2)*H1(0) + 17604*uu*G2(1,2)*H1(0) - 
     -     4320*uu**2*G2(1,2)*H1(0) - 47304*G2(2,0)*H1(0) - 
     -     (0,75816)*pi*G2(2,0)*H1(0) + 7452*G2(2,1)*H1(0) + 
     -     (0,16416)*pi*G2(2,1)*H1(0) + 8964*uu*G2(2,1)*H1(0) + 
     -     4320*uu**2*G2(2,1)*H1(0) + 8964*vv*G2(2,1)*H1(0) + 
     -     8640*uu*vv*G2(2,1)*H1(0) + 4320*vv**2*G2(2,1)*H1(0) - 
     -     42228*G2(2,2)*H1(0) - (0,77544)*pi*G2(2,2)*H1(0) + 
     -     38232*G2(2,3)*H1(0) - (0,1944)*pi*G2(2,3)*H1(0) - 
     -     8964*uu*G2(2,3)*H1(0) - 4320*uu**2*G2(2,3)*H1(0) - 
     -     8964*vv*G2(2,3)*H1(0) - 8640*uu*vv*G2(2,3)*H1(0) - 
     -     4320*vv**2*G2(2,3)*H1(0) + 19872*G2(3,0)*H1(0) + 
     -     (0,52704)*pi*G2(3,0)*H1(0) - 8964*uu*G2(3,0)*H1(0) - 
     -     4320*uu**2*G2(3,0)*H1(0) - 8964*vv*G2(3,0)*H1(0) - 
     -     8640*uu*vv*G2(3,0)*H1(0) - 4320*vv**2*G2(3,0)*H1(0) + 
     -     (0,22248)*pi*G2(3,2)*H1(0) - 8964*uu*G2(3,2)*H1(0) - 
     -     4320*uu**2*G2(3,2)*H1(0) - 8964*vv*G2(3,2)*H1(0) - 
     -     8640*uu*vv*G2(3,2)*H1(0) - 4320*vv**2*G2(3,2)*H1(0) - 
     -     864*G2(3,3)*H1(0) - (0,13824)*pi*G2(3,3)*H1(0) + 
     -     17928*uu*G2(3,3)*H1(0) + 8640*uu**2*G2(3,3)*H1(0) + 
     -     17928*vv*G2(3,3)*H1(0) + 17280*uu*vv*G2(3,3)*H1(0) + 
     -     8640*vv**2*G2(3,3)*H1(0) - 49572*G3(0,0,0)*H1(0) + 
     -     9504*G3(0,0,1)*H1(0) - 51300*G3(0,0,2)*H1(0) + 
     -     23328*G3(0,0,3)*H1(0) - 7776*G3(0,1,0)*H1(0) + 
     -     17280*G3(0,1,1)*H1(0) + 5832*G3(0,1,2)*H1(0) - 
     -     53460*G3(0,2,0)*H1(0) + 9936*G3(0,2,1)*H1(0) - 
     -     49572*G3(0,2,2)*H1(0) + 27216*G3(0,2,3)*H1(0) + 
     -     22464*G3(0,3,0)*H1(0) + 7776*G3(0,3,2)*H1(0) - 
     -     13824*G3(0,3,3)*H1(0) + 7776*G3(1,0,0)*H1(0) - 
     -     8640*G3(1,0,1)*H1(0) + 24192*G3(1,0,2)*H1(0) - 
     -     9720*G3(1,0,3)*H1(0) + 13824*G3(1,1,1)*H1(0) - 
     -     13824*G3(1,1,2)*H1(0) + 25272*G3(1,2,0)*H1(0) - 
     -     22464*G3(1,2,1)*H1(0) + 37152*G3(1,2,2)*H1(0) - 
     -     9720*G3(1,2,3)*H1(0) - 49572*G3(2,0,0)*H1(0) + 
     -     27216*G3(2,0,1)*H1(0) - 53460*G3(2,0,2)*H1(0) + 
     -     9936*G3(2,0,3)*H1(0) + 7776*G3(2,1,0)*H1(0) - 
     -     13824*G3(2,1,1)*H1(0) + 22464*G3(2,1,2)*H1(0) - 
     -     51300*G3(2,2,0)*H1(0) + 23328*G3(2,2,1)*H1(0) - 
     -     49572*G3(2,2,2)*H1(0) + 9504*G3(2,2,3)*H1(0) + 
     -     5832*G3(2,3,0)*H1(0) - 7776*G3(2,3,2)*H1(0) + 
     -     17280*G3(2,3,3)*H1(0) + 37152*G3(3,0,0)*H1(0) - 
     -     9720*G3(3,0,1)*H1(0) + 25272*G3(3,0,2)*H1(0) - 
     -     22464*G3(3,0,3)*H1(0) + 24192*G3(3,2,0)*H1(0) - 
     -     9720*G3(3,2,1)*H1(0) + 7776*G3(3,2,2)*H1(0) - 
     -     8640*G3(3,2,3)*H1(0) - 13824*G3(3,3,0)*H1(0) + 
     -     13824*G3(3,3,3)*H1(0) - 50004*H1(1) + (0,82296)*pi*H1(1) + 
     -     21060*uu*H1(1) + (0,29700)*pi*uu*H1(1) - 
     -     (0,21060)*pi*uu**2*H1(1) + 21060*vv*H1(1) - 
     -     (0,21060)*pi*uu*vv*H1(1) + 247428*zeta2*H1(1) + 
     -     (0,193428)*pi*zeta2*H1(1) - 10692*zeta3*H1(1) + 
     -     41148*G1(0)*H1(1) - (0,44712)*pi*G1(0)*H1(1) + 
     -     29700*uu*G1(0)*H1(1) - 21060*uu**2*G1(0)*H1(1) - 
     -     21060*uu*vv*G1(0)*H1(1) + 287226*zeta2*G1(0)*H1(1) - 
     -     35280*G1(1)*H1(1) + (0,6588)*pi*G1(1)*H1(1) - 
     -     21060*uu*G1(1)*H1(1) + (0,17604)*pi*uu*G1(1)*H1(1) + 
     -     21060*uu**2*G1(1)*H1(1) - (0,4320)*pi*uu**2*G1(1)*H1(1) + 
     -     4320*vv*G1(1)*H1(1) + 21060*uu*vv*G1(1)*H1(1) - 
     -     149256*zeta2*G1(1)*H1(1) + 41148*G1(2)*H1(1) - 
     -     (0,44712)*pi*G1(2)*H1(1) + 275562*zeta2*G1(2)*H1(1) - 
     -     (0,8964)*pi*uu*G1(3)*H1(1) - (0,4320)*pi*uu**2*G1(3)*H1(1) - 
     -     (0,8964)*pi*vv*G1(3)*H1(1) - (0,8640)*pi*uu*vv*G1(3)*H1(1) - 
     -     (0,4320)*pi*vv**2*G1(3)*H1(1) - 72792*zeta2*G1(3)*H1(1) - 
     -     22356*G2(0,0)*H1(1) - (0,52488)*pi*G2(0,0)*H1(1) + 
     -     24948*G2(0,1)*H1(1) + (0,7776)*pi*G2(0,1)*H1(1) + 
     -     17604*uu*G2(0,1)*H1(1) - 4320*uu**2*G2(0,1)*H1(1) - 
     -     22356*G2(0,2)*H1(1) - (0,52488)*pi*G2(0,2)*H1(1) + 
     -     (0,7776)*pi*G2(0,3)*H1(1) - 13284*G2(1,0)*H1(1) + 
     -     (0,23328)*pi*G2(1,0)*H1(1) + 17604*uu*G2(1,0)*H1(1) - 
     -     4320*uu**2*G2(1,0)*H1(1) + 25704*G2(1,1)*H1(1) - 
     -     (0,13824)*pi*G2(1,1)*H1(1) - 35208*uu*G2(1,1)*H1(1) + 
     -     8640*uu**2*G2(1,1)*H1(1) + 19872*G2(1,2)*H1(1) + 
     -     (0,38880)*pi*G2(1,2)*H1(1) - 22356*G2(2,0)*H1(1) - 
     -     (0,52488)*pi*G2(2,0)*H1(1) + 7452*G2(2,1)*H1(1) + 
     -     (0,23328)*pi*G2(2,1)*H1(1) + 8964*uu*G2(2,1)*H1(1) + 
     -     4320*uu**2*G2(2,1)*H1(1) + 8964*vv*G2(2,1)*H1(1) + 
     -     8640*uu*vv*G2(2,1)*H1(1) + 4320*vv**2*G2(2,1)*H1(1) - 
     -     22356*G2(2,2)*H1(1) - (0,52488)*pi*G2(2,2)*H1(1) - 
     -     (0,7776)*pi*G2(2,3)*H1(1) + (0,23328)*pi*G2(3,0)*H1(1) - 
     -     8964*uu*G2(3,0)*H1(1) - 4320*uu**2*G2(3,0)*H1(1) - 
     -     8964*vv*G2(3,0)*H1(1) - 8640*uu*vv*G2(3,0)*H1(1) - 
     -     4320*vv**2*G2(3,0)*H1(1) + (0,7776)*pi*G2(3,2)*H1(1) - 
     -     26244*G3(0,0,0)*H1(1) + 9504*G3(0,0,1)*H1(1) - 
     -     26244*G3(0,0,2)*H1(1) - 7776*G3(0,1,0)*H1(1) + 
     -     17280*G3(0,1,1)*H1(1) + 15552*G3(0,1,2)*H1(1) - 
     -     26244*G3(0,2,0)*H1(1) + 9936*G3(0,2,1)*H1(1) - 
     -     26244*G3(0,2,2)*H1(1) + 7776*G3(0,3,0)*H1(1) + 
     -     7776*G3(1,0,0)*H1(1) - 8640*G3(1,0,1)*H1(1) + 
     -     15552*G3(1,0,2)*H1(1) + 13824*G3(1,1,1)*H1(1) - 
     -     13824*G3(1,1,2)*H1(1) + 15552*G3(1,2,0)*H1(1) - 
     -     22464*G3(1,2,1)*H1(1) + 23328*G3(1,2,2)*H1(1) - 
     -     26244*G3(2,0,0)*H1(1) + 27216*G3(2,0,1)*H1(1) - 
     -     26244*G3(2,0,2)*H1(1) + 7776*G3(2,1,0)*H1(1) - 
     -     13824*G3(2,1,1)*H1(1) + 15552*G3(2,1,2)*H1(1) - 
     -     26244*G3(2,2,0)*H1(1) + 23328*G3(2,2,1)*H1(1) - 
     -     26244*G3(2,2,2)*H1(1) - 7776*G3(2,3,0)*H1(1) + 
     -     15552*G3(3,0,0)*H1(1) - 9720*G3(3,0,1)*H1(1) + 
     -     7776*G3(3,0,2)*H1(1) + 7776*G3(3,2,0)*H1(1) - 
     -     9720*G3(3,2,1)*H1(1) - 113004*H2(0,0) + 
     -     (0,84456)*pi*H2(0,0) - 38340*uu*H2(0,0) - 
     -     (0,26568)*pi*uu*H2(0,0) + 21060*uu**2*H2(0,0) - 
     -     4320*vv*H2(0,0) - (0,8964)*pi*vv*H2(0,0) + 
     -     21060*uu*vv*H2(0,0) - (0,8640)*pi*uu*vv*H2(0,0) - 
     -     (0,4320)*pi*vv**2*H2(0,0) - 652698*zeta2*H2(0,0) + 
     -     62100*G1(0)*H2(0,0) + (0,144072)*pi*G1(0)*H2(0,0) - 
     -     17604*uu*G1(0)*H2(0,0) + 4320*uu**2*G1(0)*H2(0,0) - 
     -     19008*G1(1)*H2(0,0) - (0,53784)*pi*G1(1)*H2(0,0) + 
     -     79596*G1(2)*H2(0,0) + (0,112968)*pi*G1(2)*H2(0,0) - 
     -     8964*uu*G1(2)*H2(0,0) - 4320*uu**2*G1(2)*H2(0,0) - 
     -     8964*vv*G1(2)*H2(0,0) - 8640*uu*vv*G1(2)*H2(0,0) - 
     -     4320*vv**2*G1(2)*H2(0,0) - 19872*G1(3)*H2(0,0) - 
     -     (0,66528)*pi*G1(3)*H2(0,0) + 17928*uu*G1(3)*H2(0,0) + 
     -     8640*uu**2*G1(3)*H2(0,0) + 17928*vv*G1(3)*H2(0,0) + 
     -     17280*uu*vv*G1(3)*H2(0,0) + 8640*vv**2*G1(3)*H2(0,0) + 
     -     86724*G2(0,0)*H2(0,0) - 13392*G2(0,1)*H2(0,0) + 
     -     90180*G2(0,2)*H2(0,0) - 29376*G2(0,3)*H2(0,0) - 
     -     25272*G2(1,0)*H2(0,0) - 38232*G2(1,2)*H2(0,0) + 
     -     72900*G2(2,0)*H2(0,0) - 15552*G2(2,1)*H2(0,0) + 
     -     72900*G2(2,2)*H2(0,0) + 3888*G2(2,3)*H2(0,0) - 
     -     50976*G2(3,0)*H2(0,0) - 24192*G2(3,2)*H2(0,0) + 
     -     13824*G2(3,3)*H2(0,0) - 41148*H2(0,1) + 
     -     (0,44712)*pi*H2(0,1) - 8640*uu*H2(0,1) - 
     -     (0,17604)*pi*uu*H2(0,1) + (0,4320)*pi*uu**2*H2(0,1) - 
     -     4320*vv*H2(0,1) - 360018*zeta2*H2(0,1) + 
     -     22356*G1(0)*H2(0,1) + (0,75816)*pi*G1(0)*H2(0,1) - 
     -     17604*uu*G1(0)*H2(0,1) + 4320*uu**2*G1(0)*H2(0,1) - 
     -     32292*G1(1)*H2(0,1) - (0,38880)*pi*G1(1)*H2(0,1) + 
     -     17604*uu*G1(1)*H2(0,1) - 4320*uu**2*G1(1)*H2(0,1) + 
     -     34776*G1(2)*H2(0,1) + (0,60264)*pi*G1(2)*H2(0,1) - 
     -     8964*uu*G1(2)*H2(0,1) - 4320*uu**2*G1(2)*H2(0,1) - 
     -     8964*vv*G1(2)*H2(0,1) - 8640*uu*vv*G1(2)*H2(0,1) - 
     -     4320*vv**2*G1(2)*H2(0,1) - (0,23328)*pi*G1(3)*H2(0,1) + 
     -     8964*uu*G1(3)*H2(0,1) + 4320*uu**2*G1(3)*H2(0,1) + 
     -     8964*vv*G1(3)*H2(0,1) + 8640*uu*vv*G1(3)*H2(0,1) + 
     -     4320*vv**2*G1(3)*H2(0,1) + 41796*G2(0,0)*H2(0,1) - 
     -     23112*G2(0,1)*H2(0,1) + 39636*G2(0,2)*H2(0,1) - 
     -     7776*G2(0,3)*H2(0,1) - 15552*G2(1,0)*H2(0,1) - 
     -     14688*G2(1,2)*H2(0,1) + 26244*G2(2,0)*H2(0,1) - 
     -     8640*G2(2,1)*H2(0,1) + 26244*G2(2,2)*H2(0,1) + 
     -     7776*G2(2,3)*H2(0,1) - 15552*G2(3,0)*H2(0,1) + 
     -     1944*G2(3,2)*H2(0,1) - 88164*H2(1,0) + 
     -     (0,81864)*pi*H2(1,0) - 29700*uu*H2(1,0) - 
     -     (0,8964)*pi*uu*H2(1,0) + 21060*uu**2*H2(1,0) - 
     -     (0,4320)*pi*uu**2*H2(1,0) - (0,8964)*pi*vv*H2(1,0) + 
     -     21060*uu*vv*H2(1,0) - (0,8640)*pi*uu*vv*H2(1,0) - 
     -     (0,4320)*pi*vv**2*H2(1,0) - 467478*zeta2*H2(1,0) + 
     -     47196*G1(0)*H2(1,0) + (0,95256)*pi*G1(0)*H2(1,0) - 
     -     6588*G1(1)*H2(1,0) - (0,39960)*pi*G1(1)*H2(1,0) - 
     -     17604*uu*G1(1)*H2(1,0) + 4320*uu**2*G1(1)*H2(1,0) + 
     -     42120*G1(2)*H2(1,0) + (0,79704)*pi*G1(2)*H2(1,0) - 
     -     (0,22248)*pi*G1(3)*H2(1,0) + 8964*uu*G1(3)*H2(1,0) + 
     -     4320*uu**2*G1(3)*H2(1,0) + 8964*vv*G1(3)*H2(1,0) + 
     -     8640*uu*vv*G1(3)*H2(1,0) + 4320*vv**2*G1(3)*H2(1,0) + 
     -     51516*G2(0,0)*H2(1,0) - 5832*G2(0,1)*H2(1,0) + 
     -     51516*G2(0,2)*H2(1,0) - 7776*G2(0,3)*H2(1,0) - 
     -     25272*G2(1,0)*H2(1,0) + 13824*G2(1,1)*H2(1,0) - 
     -     36072*G2(1,2)*H2(1,0) + 55404*G2(2,0)*H2(1,0) - 
     -     22464*G2(2,1)*H2(1,0) + 49788*G2(2,2)*H2(1,0) + 
     -     7776*G2(2,3)*H2(1,0) - 24192*G2(3,0)*H2(1,0) - 
     -     8856*G2(3,2)*H2(1,0) - 41148*H2(1,1) + 
     -     (0,44712)*pi*H2(1,1) - 248346*zeta2*H2(1,1) + 
     -     22356*G1(0)*H2(1,1) + (0,52488)*pi*G1(0)*H2(1,1) - 
     -     19872*G1(1)*H2(1,1) - (0,38880)*pi*G1(1)*H2(1,1) + 
     -     22356*G1(2)*H2(1,1) + (0,52488)*pi*G1(2)*H2(1,1) - 
     -     (0,7776)*pi*G1(3)*H2(1,1) + 26244*G2(0,0)*H2(1,1) - 
     -     15552*G2(0,1)*H2(1,1) + 26244*G2(0,2)*H2(1,1) - 
     -     15552*G2(1,0)*H2(1,1) + 13824*G2(1,1)*H2(1,1) - 
     -     23328*G2(1,2)*H2(1,1) + 26244*G2(2,0)*H2(1,1) - 
     -     15552*G2(2,1)*H2(1,1) + 26244*G2(2,2)*H2(1,1) - 
     -     7776*G2(3,0)*H2(1,1) - 81972*H3(0,0,0) - 
     -     (0,210600)*pi*H3(0,0,0) + 26568*uu*H3(0,0,0) + 
     -     8964*vv*H3(0,0,0) + 8640*uu*vv*H3(0,0,0) + 
     -     4320*vv**2*H3(0,0,0) - 137700*G1(0)*H3(0,0,0) + 
     -     52056*G1(1)*H3(0,0,0) - 105516*G1(2)*H3(0,0,0) + 
     -     64800*G1(3)*H3(0,0,0) - 22356*H3(0,0,1) - 
     -     (0,99144)*pi*H3(0,0,1) + 8964*uu*H3(0,0,1) + 
     -     4320*uu**2*H3(0,0,1) + 8964*vv*H3(0,0,1) + 
     -     8640*uu*vv*H3(0,0,1) + 4320*vv**2*H3(0,0,1) - 
     -     57348*G1(0)*H3(0,0,1) + 28512*G1(1)*H3(0,0,1) - 
     -     40068*G1(2)*H3(0,0,1) + 15552*G1(3)*H3(0,0,1) - 
     -     47196*H3(0,1,0) - (0,117504)*pi*H3(0,1,0) + 
     -     17604*uu*H3(0,1,0) - 4320*uu**2*H3(0,1,0) - 
     -     75708*G1(0)*H3(0,1,0) + 36072*G1(1)*H3(0,1,0) - 
     -     57348*G1(2)*H3(0,1,0) + 24192*G1(3)*H3(0,1,0) - 
     -     22356*H3(0,1,1) - (0,60264)*pi*H3(0,1,1) - 
     -     34020*G1(0)*H3(0,1,1) + 23328*G1(1)*H3(0,1,1) - 
     -     34020*G1(2)*H3(0,1,1) + 7776*G1(3)*H3(0,1,1) - 
     -     84348*H3(1,0,0) - (0,138024)*pi*H3(1,0,0) + 
     -     8964*uu*H3(1,0,0) + 4320*uu**2*H3(1,0,0) + 
     -     8964*vv*H3(1,0,0) + 8640*uu*vv*H3(1,0,0) + 
     -     4320*vv**2*H3(1,0,0) - 94068*G1(0)*H3(1,0,0) + 
     -     38232*G1(1)*H3(1,0,0) - 76788*G1(2)*H3(1,0,0) + 
     -     24192*G1(3)*H3(1,0,0) - 34668*H3(1,0,1) - 
     -     (0,79704)*pi*H3(1,0,1) + 8964*uu*H3(1,0,1) + 
     -     4320*uu**2*H3(1,0,1) + 8964*vv*H3(1,0,1) + 
     -     8640*uu*vv*H3(1,0,1) + 4320*vv**2*H3(1,0,1) - 
     -     43740*G1(0)*H3(1,0,1) + 14688*G1(1)*H3(1,0,1) - 
     -     24300*G1(2)*H3(1,0,1) - 1944*G1(3)*H3(1,0,1) - 
     -     42228*H3(1,1,0) - (0,60264)*pi*H3(1,1,0) - 
     -     53460*G1(0)*H3(1,1,0) + 36072*G1(1)*H3(1,1,0) - 
     -     45684*G1(2)*H3(1,1,0) + 8856*G1(3)*H3(1,1,0) - 
     -     22356*H3(1,1,1) - (0,52488)*pi*H3(1,1,1) - 
     -     26244*G1(0)*H3(1,1,1) + 23328*G1(1)*H3(1,1,1) - 
     -     26244*G1(2)*H3(1,1,1) + 202500*H4(0,0,0,0) + 
     -     72900*H4(0,0,0,1) + 99900*H4(0,0,1,0) + 41796*H4(0,0,1,1) + 
     -     118260*H4(0,1,0,0) + 41796*H4(0,1,0,1) + 62316*H4(0,1,1,0) + 
     -     26244*H4(0,1,1,1) + 137700*H4(1,0,0,0) + 43956*H4(1,0,0,1) + 
     -     82836*H4(1,0,1,0) + 35964*H4(1,0,1,1) + 63396*H4(1,1,0,0) + 
     -     6804*H4(1,1,0,1) + 55404*H4(1,1,1,0) + 26244*H4(1,1,1,1))/
     -   3888.  

      a2finite = a2finitelogcoeff(1)*lnrat(-hmass**2,musq)**4
     & + a2finitelogcoeff(2)*lnrat(-hmass**2,musq)**3
     & + a2finitelogcoeff(3)*lnrat(-hmass**2,musq)**2
     & + a2finitelogcoeff(4)*lnrat(-hmass**2,musq)
     & + a2finitelogcoeff(5)

      a2finite = a2born*a2finite

   99 format(a35,e20.12)
      end
