!=================================================================================================

        double precision function YLP0qiqi(z,x,f,f0)
        implicit none
        double precision z,x,f,f0
        include 'paramqt.f'
        YLP0qiqi = -(CF*(z + dble(1))*PlusDist(-2,z,x,f,f0)) + CF*dble(dble(3)/dble(2))*PlusDist(-1,z,x,f,f0) + CF*dble(2)*PlusDist(0,z,x,f,f0)
        return
        end

        double precision function YLP0qg(z,x,f,f0)
        implicit none
        double precision z,x,f,f0
        include 'paramqt.f'
        YLP0qg = TF*(z*dble(-2) + dble(1) + z**2*dble(2))*PlusDist(-2,z,x,f,f0)
        return
        end

        double precision function YLP0gg(z,x,f,f0)
        implicit none
        double precision z,x,f,f0
        include 'paramqt.f'
        YLP0gg = (CA*(z**2 + z*dble(-2) + z**3*dble(-1) + dble(1))*dble(2)*PlusDist(-2,z,x,f,f0))/z + 
     -  beta0*dble(dble(1)/dble(2))*PlusDist(-1,z,x,f,f0) + CA*dble(2)*PlusDist(0,z,x,f,f0)
        return
        end

        double precision function YLP0gq(z,x,f,f0)
        implicit none
        double precision z,x,f,f0
        include 'paramqt.f'
        YLP0gq = (CF*(z**2 + z*dble(-2) + dble(2))*PlusDist(-2,z,x,f,f0))/z
        return
        end

!-----------------------------------------------

        double precision function YLP1qiqi(z,x,f,f0)
        implicit none 
        double precision z,x,f,f0
        include 'paramqt.f'
        YLP1qiqi = (CF*dble(dble(-1)/dble(72))*((z + dble(-1))*(z*dble(9)*
     -           (G1 + G1*z + beta0*dble(-8) + CF*dble(-40)*(z + dble(-1)) + beta0*z*dble(8) + CA*(z + dble(-1))*dble(24)) + 
     -          TF*dble(16)*(z**2*dble(-27) + dble(-10) + z*dble(9) + z**3*dble(28))) + 
     -       z*(CA*(z**2 + dble(1)) + (z**2 + dble(-1))*(CF + TF*dble(2)))*dble(36)*log(z)**2 + 
     -       z*dble(-12)*log(z)*(TF*z*dble(-24) + CF*z**2*dble(-12) + TF*dble(-6) + beta0*dble(-3) + beta0*z**2*dble(-3) + 
     -          CA*(z**2 + dble(-1))*dble(6) + CF*z*dble(12) + TF*z**2*dble(14) + TF*z**3*dble(16) + CF*dble(18) + 
     -          CF*(z**2 + dble(1))*dble(12)*log(-z + dble(1))))*PlusDist(-2,z,x,f,f0))/(z*(z + dble(-1))) + 
     -  CF*(CA*(Zeta3*dble(-72) + dble(6)) + beta0*(dble(3) + Zeta2*dble(24)) + CF*(Zeta2*dble(-72) + dble(9) + Zeta3*dble(144)))*
     -   dble(dble(1)/dble(24))*PlusDist(-1,z,x,f,f0) + CF*G1*dble(dble(1)/dble(4))*PlusDist(0,z,x,f,f0)
        return
        end

        double precision function YLP1qiqj(z,x,f,f0)
        implicit none
        double precision z,x,f,f0
        include 'paramqt.f'
        YLP1qiqj = (CF*TF*dble(dble(1)/dble(9))*(z**3*dble(-56) + z*dble(-18) + dble(20) + z**2*dble(54) + 
     -      z*dble(3)*(dble(3) + z**2*dble(8) + z*dble(15))*log(z) + z*dble(-9)*(z + dble(1))*log(z)**2)*PlusDist(-2,z,x,f,f0))/z
        return
        end

        double precision function YLP1qiqbi(z,x,f,f0)
        implicit none
        double precision z,x,f,f0
        include 'paramqt.f'
        YLP1qiqbi = (CF*dble(dble(1)/dble(18))*(TF*z**4*dble(-112) + CF*z**3*dble(-72) + CA*z*dble(-36) + CF*z*Zeta2*dble(-36) + 
     -      CF*z**3*Zeta2*dble(-36) + TF*z**3*dble(-4) + TF*z*dble(4) + CA*z*Zeta2*dble(18) + CA*z**3*Zeta2*dble(18) + 
     -      CA*z**3*dble(36) + TF*dble(40) + CF*z*dble(72) + TF*z**2*dble(72) + 
     -      z*(CA + CF*dble(-2))*(z**2 + dble(1))*dble(36)*Li2(-z) + 
     -      z*dble(-9)*(CA*(z**2 + dble(1)) + (TF*(z + dble(1))**2 - CF*(z**2 + dble(1)))*dble(2))*log(z)**2 + 
     -      z*dble(6)*log(z)*((z + dble(1))*(CA*dble(-3)*(z + dble(1)) + CF*(z + dble(1))*dble(6) + 
     -            TF*(dble(3) + z**2*dble(8) + z*dble(15))) + (CA + CF*dble(-2))*(z**2 + dble(1))*dble(6)*log(z + dble(1))))*
     -    PlusDist(-2,z,x,f,f0))/(z*(z + dble(1)))
        return
        end

        double precision function YLP1qg(z,x,f,f0)
        implicit none
        double precision z,x,f,f0
        include 'paramqt.f'
        YLP1qg = (TF*dble(dble(-1)/dble(18))*(CA*z**2*dble(-450) + CF*z**3*dble(-180) + CF*z*dble(-126) + CF*z**2*Zeta2*dble(-72) + 
     -      CA*dble(-40) + CA*z*dble(36) + CF*z*Zeta2*dble(36) + CA*z**2*Zeta2*dble(72) + CF*z**3*Zeta2*dble(72) + 
     -      CF*z**2*dble(261) + CA*z**3*dble(436) + CA*z*(dble(1) + z*dble(2) + z**2*dble(2))*dble(36)*Li2(-z) + 
     -      CA*z**3*dble(-264)*log(z) + CA*z**2*dble(-144)*log(z) + CF*z**3*dble(-72)*log(z) + CF*z*dble(-27)*log(z) + 
     -      CA*z*dble(-18)*log(z) + CF*z**2*dble(36)*log(z) + CF*z**3*dble(-36)*log(z)**2 + CF*z*dble(-9)*log(z)**2 + 
     -      CA*z*dble(18)*log(z)**2 + CF*z**2*dble(18)*log(z)**2 + CA*z**2*dble(36)*log(z)**2 + 
     -      z*dble(-36)*((CA - CF)*z*(z + dble(-1))*dble(2) - CF*(z*dble(-2) + dble(1) + z**2*dble(2))*log(z))*log(-z + dble(1)) + 
     -      (CA - CF)*z*(z*dble(-2) + dble(1) + z**2*dble(2))*dble(18)*log(-z + dble(1))**2 + 
     -      CA*z*dble(36)*log(z)*log(z + dble(1)) + CA*z**2*dble(72)*log(z)*log(z + dble(1)) + 
     -      CA*z**3*dble(72)*log(z)*log(z + dble(1)))*PlusDist(-2,z,x,f,f0))/z
        return
        end

        double precision function YLP1gg(z,x,f,f0)
        implicit none
        double precision z,x,f,f0
        include 'paramqt.f'
        YLP1gg = (dble(dble(1)/dble(36))*((z + dble(-1))*(CA*dble(-3)*
     -           (G1*(z + z**2 + z**4 + dble(-1))*dble(3) + beta0*(dble(-26) + z*dble(-8) + z**3*dble(8) + z**4*dble(26))) + 
     -          CF*nf*TF*(z*dble(-11) + z**2*dble(-6) + dble(1) + z**4*dble(5) + z**3*dble(11))*dble(48) + 
     -          CA**2*dble(2)*(dble(-277) + z**3*dble(-65) + Zeta2*(z + z**2 + dble(1))**2*dble(36) + z*dble(65) + z**4*dble(277))
     -          ) + CA**2*(z + dble(-1))*(z + z**2 + dble(1))**2*dble(144)*Li2(-z) + 
     -       z*dble(-72)*(CF*nf*TF*(z + dble(-1))*(z + dble(1))**2 + CA**2*(z + z**2*dble(-1) + dble(1))**2)*log(z)**2 + 
     -       dble(12)*log(z)*(CA**2*(z + dble(1))*(-z + z**2 + dble(1))**2*dble(12)*log(-z + dble(1)) + 
     -          (z + dble(-1))*(-(z*(z + dble(1))*(beta0*CA*dble(-3)*(z + dble(1)) + CF*nf*TF*(dble(3) + z*dble(5))*dble(6) + 
     -                  CA**2*dble(4)*(dble(9) + z**2*dble(11)))) + CA**2*(z + z**2 + dble(1))**2*dble(12)*log(z + dble(1)))))*
     -     PlusDist(-2,z,x,f,f0))/(z*(z**2 + dble(-1))) + 
     -  (beta0*CA - CF*nf*TF + CA**2*(dble(-1) + Zeta3*dble(3)))*PlusDist(-1,z,x,f,f0) + 
     -  CA*G1*dble(dble(1)/dble(4))*PlusDist(0,z,x,f,f0)
        return
        end

        double precision function YLP1gq(z,x,f,f0)
        implicit none
        double precision z,x,f,f0
        include 'paramqt.f'
        YLP1gq = (CF*dble(dble(1)/dble(18))*(CA*dble(-202) + CA*z**2*dble(-102) + CF*z**2*dble(-63) + beta0*z*dble(-60) + CF*z*dble(-45) + 
     -      beta0*z**2*dble(48) + beta0*dble(60) + CA*z*Zeta2*dble(72) + CA*z**3*dble(88) + CA*z*dble(258) + 
     -      CA*(z**2 + dble(2) + z*dble(2))*dble(36)*Li2(-z) + CA*z*dble(-216)*log(z) + CA*z**2*dble(-90)*log(z) + 
     -      CA*z**3*dble(-48)*log(z) + CF*z*dble(36)*log(z) + CF*z**2*dble(63)*log(z) + CF*z*dble(-18)*log(z)**2 + 
     -      CF*z**2*dble(9)*log(z)**2 + CA*z**2*dble(18)*log(z)**2 + CA*z*dble(36)*log(z)**2 + 
     -      dble(-18)*(-(beta0*z**2) + beta0*dble(-2) + CA*z**2*dble(-2) + beta0*z*dble(2) + 
     -         CF*(z*dble(-6) + z**2*dble(5) + dble(6)) + CA*dble(2)*(z**2 + z*dble(-2) + dble(2))*log(z))*log(-z + dble(1)) + 
     -      (CA - CF)*(z**2 + z*dble(-2) + dble(2))*dble(18)*log(-z + dble(1))**2 + CA*z**2*dble(36)*log(z)*log(z + dble(1)) + 
     -      CA*dble(72)*log(z)*log(z + dble(1)) + CA*z*dble(72)*log(z)*log(z + dble(1)))*PlusDist(-2,z,x,f,f0))/z
        return
        end

!===============================================

        double precision function YLI1qiqi(z,x,f,f0)
        implicit none
        double precision z,x,f,f0
        include 'paramqt.f'
        YLI1qiqi = CF*dble(-2)*(z + dble(-1))*PlusDist(-2,z,x,f,f0) - CF*Zeta2*PlusDist(-1,z,x,f,f0)
        return
        end

        double precision function YLI1qg(z,x,f,f0)
        implicit none
        double precision z,x,f,f0
        include 'paramqt.f'
        YLI1qg = TF*z*dble(-4)*(z + dble(-1))*PlusDist(-2,z,x,f,f0)
        return
        end

        double precision function YLI1gg(z,x,f,f0)
        implicit none
        double precision z,x,f,f0
        include 'paramqt.f'
        YLI1gg = -(CA*Zeta2*PlusDist(-1,z,x,f,f0))
        return
        end

        double precision function YLI1gq(z,x,f,f0)
        implicit none
        double precision z,x,f,f0
        include 'paramqt.f'
        YLI1gq = CF*z*dble(2)*PlusDist(-2,z,x,f,f0)
        return
        end

        double precision function YLIi1gg(z,x,f,f0)
        implicit none
        double precision z,x,f,f0
        include 'paramqt.f'
        YLIi1gg = (CA*(z + dble(-1))*dble(4)*PlusDist(-2,z,x,f,f0))/z
        return
        end

        double precision function YLIi1gq(z,x,f,f0)
        implicit none
        double precision z,x,f,f0
        include 'paramqt.f'
        YLIi1gq = (CF*(z + dble(-1))*dble(4)*PlusDist(-2,z,x,f,f0))/z
        return
        end

!-----------------------------------------------

        double precision function YLI2qiqi(z,x,f,f0)
        implicit none
        double precision z,x,f,f0
        include 'paramqt.f'
        YLI2qiqi = (CF*dble(dble(1)/dble(54))*(CF*z**2*dble(-2376) + TF*z**2*dble(-2376) + CA*z*dble(-1600) + CF*z*Zeta3*dble(-1080) + 
     -       CF*z**3*Zeta3*dble(-1080) + TF*z*Zeta2*dble(-720) + TF*z**3*Zeta2*dble(-720) + TF*dble(-688) + 
     -       CA*z**2*Zeta2*dble(-648) + TF*z**4*dble(-544) + CF*z*Zeta2*dble(-324) + CF*z**3*Zeta2*dble(-324) + 
     -       CA*z**3*Zeta3*dble(-324) + nf*TF*z**3*dble(-152) + nf*TF*z**2*dble(-144) + CA*z**3*dble(16) + TF*Zeta2*dble(288) + 
     -       TF*z**4*Zeta2*dble(288) + nf*TF*z*dble(296) + CA*z*Zeta2*dble(324) + CA*z**3*Zeta2*dble(324) + 
     -       CF*z**2*Zeta2*dble(648) + TF*z**2*Zeta2*dble(864) + CF*z*dble(1188) + CF*z**3*dble(1188) + CA*z*Zeta3*dble(1188) + 
     -       CA*z**2*dble(1584) + TF*z**3*dble(1660) + TF*z*dble(1948) + CA*z*dble(-432)*Li3(z) + CA*z**3*dble(-432)*Li3(z) + 
     -       CF*z*dble(1080)*Li3(z) + CF*z**3*dble(1080)*Li3(z) + (CA - CF)*z*(z**2 + dble(1))*dble(216)*Li3(-z + dble(1)) + 
     -       CF*z**3*dble(-1728)*log(z) + TF*z**3*dble(-1488)*log(z) + CF*z*dble(-540)*log(z) + TF*z*dble(-504)*log(z) + 
     -       CA*z**2*dble(-432)*log(z) + nf*TF*z*dble(-120)*log(z) + nf*TF*z**3*dble(-120)*log(z) + CA*z*dble(348)*log(z) + 
     -       TF*z**4*dble(768)*log(z) + CA*z**3*dble(996)*log(z) + TF*z**2*dble(1224)*log(z) + CF*z**2*dble(1404)*log(z) + 
     -       CF*z*dble(-162)*log(z)**2 + TF*z**4*dble(-144)*log(z)**2 + CF*z**2*dble(-108)*log(z)**2 + 
     -       nf*TF*z*dble(-36)*log(z)**2 + nf*TF*z**3*dble(-36)*log(z)**2 + CA*z**3*dble(-9)*log(z)**2 + TF*z*dble(54)*log(z)**2 + 
     -       TF*z**3*dble(90)*log(z)**2 + CA*z*dble(99)*log(z)**2 + CA*z**2*dble(108)*log(z)**2 + CF*z**3*dble(108)*log(z)**2 + 
     -       TF*z*dble(-36)*log(z)**3 + CF*z*dble(-18)*log(z)**3 + CA*z*dble(18)*log(z)**3 + CA*z**3*dble(18)*log(z)**3 + 
     -       CF*z**3*dble(18)*log(z)**3 + TF*z**3*dble(36)*log(z)**3 + CA*z*Zeta2*dble(-216)*log(-z + dble(1)) + 
     -       CA*z**3*Zeta2*dble(-216)*log(-z + dble(1)) + CA*z**2*dble(-108)*log(-z + dble(1)) + 
     -       CF*z**3*dble(-108)*log(-z + dble(1)) + CF*z**2*dble(108)*log(-z + dble(1)) + CA*z**3*dble(108)*log(-z + dble(1)) + 
     -       CF*z*Zeta2*dble(216)*log(-z + dble(1)) + CF*z**3*Zeta2*dble(216)*log(-z + dble(1)) + 
     -       CF*z**2*dble(-1296)*log(z)*log(-z + dble(1)) + TF*z**2*dble(-864)*log(z)*log(-z + dble(1)) + 
     -       TF*dble(-288)*log(z)*log(-z + dble(1)) + TF*z**4*dble(-288)*log(z)*log(-z + dble(1)) + 
     -       CA*z*dble(-216)*log(z)*log(-z + dble(1)) + CA*z**3*dble(-216)*log(z)*log(-z + dble(1)) + 
     -       CA*z**2*dble(432)*log(z)*log(-z + dble(1)) + CF*z*dble(648)*log(z)*log(-z + dble(1)) + 
     -       CF*z**3*dble(648)*log(z)*log(-z + dble(1)) + TF*z*dble(720)*log(z)*log(-z + dble(1)) + 
     -       TF*z**3*dble(720)*log(z)*log(-z + dble(1)) + CF*z*dble(-108)*log(z)**2*log(-z + dble(1)) + 
     -       CF*z**3*dble(-108)*log(z)**2*log(-z + dble(1)) + CF*z*dble(-324)*log(z)*log(-z + dble(1))**2 + 
     -       CF*z**3*dble(-324)*log(z)*log(-z + dble(1))**2 + CA*z*dble(216)*log(z)*log(-z + dble(1))**2 + 
     -       CA*z**3*dble(216)*log(z)*log(-z + dble(1))**2 + 
     -       dble(-72)*Li2(z)*((z + dble(-1))**2*(z*(CA + CF*dble(-2))*dble(3) + TF*(z*dble(-2) + dble(4) + z**2*dble(4))) + 
     -          z*dble(-3)*(CA + CF*dble(-3))*(z**2 + dble(1))*log(z) + (CA - CF)*z*dble(-3)*(z**2 + dble(1))*log(-z + dble(1))))*
     -     PlusDist(-2,z,x,f,f0))/(z*(z + dble(-1))) + 
     -  CF*(nf*TF*dble(8)*(dble(-328) + Zeta3*dble(252) + Zeta2*dble(270)) + CF*Zeta4*dble(810) + 
     -     CA*dble(2)*(Zeta2*dble(-3618) + Zeta3*dble(-2772) + Zeta4*dble(1620) + dble(4856)))*dble(dble(1)/dble(648))*
     -   PlusDist(-1,z,x,f,f0) + CF*(nf*TF*dble(56) + CA*(dble(-202) + Zeta3*dble(189)))*dble(dble(4)/dble(27))*PlusDist(0,z,x,f,f0)
        return
        end

        double precision function YLI2qiqj(z,x,f,f0)
        implicit none
        double precision z,x,f,f0
        include 'paramqt.f'
        YLI2qiqj = (CF*TF*dble(dble(1)/dble(27))*((z + dble(-1))*dble(2)*
     -       (dble(-172) + z**2*dble(-136) + Zeta2*(-z + dble(2) + z**2*dble(2))*dble(36) + z*dble(143)) + 
     -      dble(-72)*(z**2*dble(-3) + dble(-2) + z**3*dble(2) + z*dble(3))*Li2(z) + 
     -      z*dble(-9)*(dble(3) + z*dble(3) + z**2*dble(8))*log(z)**2 + z*(z + dble(1))*dble(18)*log(z)**3 + 
     -      dble(-12)*log(z)*(z*(z**2*dble(-32) + dble(-21) + z*dble(30)) + 
     -         (z**2*dble(-3) + dble(-2) + z**3*dble(2) + z*dble(3))*dble(6)*log(-z + dble(1))))*PlusDist(-2,z,x,f,f0))/z
        return
        end

        double precision function YLI2qiqbi(z,x,f,f0)
        implicit none
        double precision z,x,f,f0
        include 'paramqt.f'
        YLI2qiqbi = (CF*dble(dble(1)/dble(27))*(CF*z**3*dble(-810) + CA*z*dble(-405) + CF*z*Zeta2*dble(-324) + CF*z*Zeta3*dble(-324) + 
     -      CF*z**3*Zeta3*dble(-324) + TF*z*dble(-286) + TF*z**4*dble(-272) + CF*z**2*Zeta2*dble(-216) + TF*Zeta2*dble(-144) + 
     -      TF*z**2*dble(-72) + TF*z**3*Zeta2*dble(-72) + CA*z**3*Zeta2*dble(-54) + TF*z*Zeta2*dble(72) + CA*z**2*Zeta2*dble(108) + 
     -      CF*z**3*Zeta2*dble(108) + TF*z**4*Zeta2*dble(144) + CA*z*Zeta2*dble(162) + CA*z*Zeta3*dble(162) + 
     -      CA*z**3*Zeta3*dble(162) + TF*z**3*dble(286) + TF*dble(344) + CA*z**3*dble(405) + CF*z*dble(810) + 
     -      CA*z*dble(-324)*Li3(-z) + CA*z**3*dble(-324)*Li3(-z) + CF*z*dble(648)*Li3(-z) + CF*z**3*dble(648)*Li3(-z) + 
     -      CA*z*dble(-216)*Li3(z) + CA*z**3*dble(-216)*Li3(z) + CF*z*dble(432)*Li3(z) + CF*z**3*dble(432)*Li3(z) + 
     -      CA*z*dble(-216)*Li3(1/(z + dble(1))) + CA*z**3*dble(-216)*Li3(1/(z + dble(1))) + CF*z*dble(432)*Li3(1/(z + dble(1))) + 
     -      CF*z**3*dble(432)*Li3(1/(z + dble(1))) + CA*z**2*dble(-378)*log(z) + CA*z**3*dble(-297)*log(z) + 
     -      TF*z**2*dble(-108)*log(z) + CA*z*dble(-81)*log(z) + TF*z**3*dble(24)*log(z) + CF*z*dble(162)*log(z) + 
     -      TF*z*dble(252)*log(z) + TF*z**4*dble(384)*log(z) + CF*z**3*dble(594)*log(z) + CF*z**2*dble(756)*log(z) + 
     -      TF*z**3*dble(-99)*log(z)**2 + TF*z**4*dble(-72)*log(z)**2 + TF*z**2*dble(-54)*log(z)**2 + TF*z*dble(-27)*log(z)**2 + 
     -      CF*z*dble(-18)*log(z)**3 + CF*z**3*dble(-18)*log(z)**3 + CA*z*dble(9)*log(z)**3 + CA*z**3*dble(9)*log(z)**3 + 
     -      TF*z*dble(18)*log(z)**3 + TF*z**3*dble(18)*log(z)**3 + TF*z**2*dble(36)*log(z)**3 + 
     -      z*(CA + CF*dble(-2))*dble(108)*Li2(-z)*((z + dble(1))**2 + (z**2 + dble(1))*log(z)) + 
     -      dble(-36)*Li2(z)*((z**2 + dble(-1))*(CA*z*dble(-3) + TF*z*dble(-2) + TF*dble(4) + TF*z**2*dble(4) + CF*z*dble(6)) + 
     -         z*dble(-3)*(CA + CF*dble(-2))*(z**2 + dble(1))*log(z)) + CF*z**3*dble(-216)*log(z)*log(-z + dble(1)) + 
     -      TF*z**4*dble(-144)*log(z)*log(-z + dble(1)) + CA*z*dble(-108)*log(z)*log(-z + dble(1)) + 
     -      TF*z*dble(-72)*log(z)*log(-z + dble(1)) + TF*z**3*dble(72)*log(z)*log(-z + dble(1)) + 
     -      CA*z**3*dble(108)*log(z)*log(-z + dble(1)) + TF*dble(144)*log(z)*log(-z + dble(1)) + 
     -      CF*z*dble(216)*log(z)*log(-z + dble(1)) + CA*z*Zeta2*dble(-108)*log(z + dble(1)) + 
     -      CA*z**3*Zeta2*dble(-108)*log(z + dble(1)) + CF*z*Zeta2*dble(216)*log(z + dble(1)) + 
     -      CF*z**3*Zeta2*dble(216)*log(z + dble(1)) + CF*z**2*dble(-432)*log(z)*log(z + dble(1)) + 
     -      CF*z*dble(-216)*log(z)*log(z + dble(1)) + CF*z**3*dble(-216)*log(z)*log(z + dble(1)) + 
     -      CA*z*dble(108)*log(z)*log(z + dble(1)) + CA*z**3*dble(108)*log(z)*log(z + dble(1)) + 
     -      CA*z**2*dble(216)*log(z)*log(z + dble(1)) + CA*z*dble(-54)*log(z)**2*log(z + dble(1)) + 
     -      CA*z**3*dble(-54)*log(z)**2*log(z + dble(1)) + CF*z*dble(108)*log(z)**2*log(z + dble(1)) + 
     -      CF*z**3*dble(108)*log(z)**2*log(z + dble(1)) + CF*z*dble(-72)*log(z + dble(1))**3 + 
     -      CF*z**3*dble(-72)*log(z + dble(1))**3 + CA*z*dble(36)*log(z + dble(1))**3 + CA*z**3*dble(36)*log(z + dble(1))**3)*
     -    PlusDist(-2,z,x,f,f0))/(z*(z + dble(1)))
        return
        end

        double precision function YLI2qg(z,x,f,f0)
        implicit none
        double precision z,x,f,f0
        include 'paramqt.f'
        YLI2qg = (TF*dble(dble(1)/dble(54))*(CF*z**3*dble(-3888) + CF*z**2*Zeta3*dble(-3456) + CA*z**2*Zeta2*dble(-1296) + 
     -      CA*z**3*Zeta3*dble(-1296) + CA*z*dble(-1260) + CA*z**3*dble(-1192) + CF*z*dble(-702) + CF*z**2*Zeta2*dble(-648) + 
     -      CA*z*Zeta3*dble(-648) + CA*Zeta2*dble(-288) + CA*z*Zeta2*dble(432) + CF*z**3*Zeta2*dble(648) + CA*dble(688) + 
     -      CA*z**2*dble(1548) + CA*z**3*Zeta2*dble(1584) + CF*z*Zeta3*dble(1728) + CF*z**3*Zeta3*dble(3456) + CF*z**2*dble(4050) + 
     -      CA*z*dble(648)*Li3(-z) + CA*z**2*dble(1296)*Li3(-z) + CA*z**3*dble(1296)*Li3(-z) + CF*z**3*dble(-432)*Li3(z) + 
     -      CF*z*dble(-216)*Li3(z) + CF*z**2*dble(432)*Li3(z) + CA*z**2*dble(1728)*Li3(z) + CA*z**2*dble(-432)*Li3(-z + dble(1)) + 
     -      CF*z**3*dble(-432)*Li3(-z + dble(1)) + CF*z*dble(-216)*Li3(-z + dble(1)) + CA*z*dble(216)*Li3(-z + dble(1)) + 
     -      CF*z**2*dble(432)*Li3(-z + dble(1)) + CA*z**3*dble(432)*Li3(-z + dble(1)) + CA*z*dble(432)*Li3(1/(z + dble(1))) + 
     -      CA*z**2*dble(864)*Li3(1/(z + dble(1))) + CA*z**3*dble(864)*Li3(1/(z + dble(1))) + CA*z**2*dble(-720)*log(z) + 
     -      CF*z**3*dble(-432)*log(z) + CF*z*dble(432)*log(z) + CA*z*dble(504)*log(z) + CF*z**2*dble(810)*log(z) + 
     -      CA*z**3*dble(1632)*log(z) + CA*z**3*dble(-792)*log(z)**2 + CF*z**3*dble(-216)*log(z)**2 + CA*z*dble(-54)*log(z)**2 + 
     -      CF*z*dble(27)*log(z)**2 + CA*z**2*dble(216)*log(z)**2 + CF*z**2*dble(324)*log(z)**2 + CF*z**3*dble(-72)*log(z)**3 + 
     -      CF*z*dble(-18)*log(z)**3 + CA*z*dble(36)*log(z)**3 + CF*z**2*dble(36)*log(z)**3 + CA*z**2*dble(72)*log(z)**3 + 
     -      CA*z*dble(-216)*Li2(-z)*(z*dble(-2)*(z + dble(1)) + (dble(1) + z*dble(2) + z**2*dble(2))*log(z)) + 
     -      CA*z**3*dble(-432)*log(-z + dble(1)) + CF*z**2*Zeta2*dble(-432)*log(-z + dble(1)) + 
     -      CA*z**3*Zeta2*dble(-432)*log(-z + dble(1)) + CF*z**2*dble(-324)*log(-z + dble(1)) + 
     -      CA*z*Zeta2*dble(-216)*log(-z + dble(1)) + CF*z*Zeta2*dble(216)*log(-z + dble(1)) + 
     -      CA*z**2*dble(324)*log(-z + dble(1)) + CF*z**3*dble(432)*log(-z + dble(1)) + CA*z**2*Zeta2*dble(432)*log(-z + dble(1)) + 
     -      CF*z**3*Zeta2*dble(432)*log(-z + dble(1)) + CA*z**3*dble(-1584)*log(z)*log(-z + dble(1)) + 
     -      CA*z*dble(-432)*log(z)*log(-z + dble(1)) + CF*z**2*dble(-432)*log(z)*log(-z + dble(1)) + 
     -      CA*dble(288)*log(z)*log(-z + dble(1)) + CF*z**3*dble(432)*log(z)*log(-z + dble(1)) + 
     -      CA*z**2*dble(1728)*log(z)*log(-z + dble(1)) + CF*z**2*dble(-216)*log(z)**2*log(-z + dble(1)) + 
     -      CF*z*dble(108)*log(z)**2*log(-z + dble(1)) + CF*z**3*dble(216)*log(z)**2*log(-z + dble(1)) + 
     -      CA*z**2*dble(-216)*log(-z + dble(1))**2 + CF*z**3*dble(-216)*log(-z + dble(1))**2 + 
     -      CF*z**2*dble(216)*log(-z + dble(1))**2 + CA*z**3*dble(216)*log(-z + dble(1))**2 + 
     -      CA*z**2*dble(-432)*log(z)*log(-z + dble(1))**2 + CF*z**3*dble(-216)*log(z)*log(-z + dble(1))**2 + 
     -      CF*z*dble(-108)*log(z)*log(-z + dble(1))**2 + CA*z*dble(216)*log(z)*log(-z + dble(1))**2 + 
     -      CF*z**2*dble(216)*log(z)*log(-z + dble(1))**2 + CA*z**3*dble(432)*log(z)*log(-z + dble(1))**2 + 
     -      CA*z**2*dble(-72)*log(-z + dble(1))**3 + CF*z**3*dble(-72)*log(-z + dble(1))**3 + CF*z*dble(-36)*log(-z + dble(1))**3 + 
     -      CA*z*dble(36)*log(-z + dble(1))**3 + CF*z**2*dble(72)*log(-z + dble(1))**3 + CA*z**3*dble(72)*log(-z + dble(1))**3 + 
     -      dble(72)*Li2(z)*(CA*dble(2)*(z**3*dble(-11) + z*dble(-3) + dble(2) + z**2*dble(12)) + 
     -         z*(CF + CA*z*dble(-4) + CF*z*dble(-2) + CF*z**2*dble(2))*dble(3)*log(z) + 
     -         (CA - CF)*z*(z*dble(-2) + dble(1) + z**2*dble(2))*dble(3)*log(-z + dble(1))) + 
     -      CA*z*Zeta2*dble(216)*log(z + dble(1)) + CA*z**2*Zeta2*dble(432)*log(z + dble(1)) + 
     -      CA*z**3*Zeta2*dble(432)*log(z + dble(1)) + CA*z**2*dble(432)*log(z)*log(z + dble(1)) + 
     -      CA*z**3*dble(432)*log(z)*log(z + dble(1)) + CA*z*dble(108)*log(z)**2*log(z + dble(1)) + 
     -      CA*z**2*dble(216)*log(z)**2*log(z + dble(1)) + CA*z**3*dble(216)*log(z)**2*log(z + dble(1)) + 
     -      CA*z**2*dble(-144)*log(z + dble(1))**3 + CA*z**3*dble(-144)*log(z + dble(1))**3 + CA*z*dble(-72)*log(z + dble(1))**3)*
     -    PlusDist(-2,z,x,f,f0))/z
        return
        end

        double precision function YLI2gg(z,x,f,f0)
        implicit none
        double precision z,x,f,f0
        include 'paramqt.f'
        YLI2gg = (dble(dble(-1)/dble(27))*(CA**2*z**5*dble(-3340) + CA**2*dble(-3160) + CA**2*z*Zeta3*dble(-1836) + 
     -       CA**2*z**4*Zeta3*dble(-1728) + CF*nf*TF*z**2*dble(-1656) + CF*nf*TF*z**3*dble(-1656) + CA**2*z*Zeta2*dble(-864) + 
     -       CA**2*z**4*Zeta2*dble(-864) + CA*nf*TF*z*dble(-664) + CA*nf*TF*z**4*dble(-440) + CA**2*z**3*dble(-364) + 
     -       CF*nf*TF*dble(-72) + CF*nf*TF*z**5*dble(-72) + CA*nf*TF*z**2*dble(-44) + CA**2*z**2*Zeta2*dble(72) + 
     -       CA**2*z**3*Zeta2*dble(72) + CA*nf*TF*z**3*dble(108) + CA**2*z**2*dble(120) + CA*nf*TF*dble(484) + 
     -       CA*nf*TF*z**5*dble(556) + CA**2*Zeta2*dble(792) + CA**2*z**5*Zeta2*dble(792) + CA**2*z**2*Zeta3*dble(972) + 
     -       CA**2*z**3*Zeta3*dble(1080) + CA**2*z**5*Zeta3*dble(1080) + CF*nf*TF*z*dble(1728) + CF*nf*TF*z**4*dble(1728) + 
     -       CA**2*Zeta3*dble(1728) + CA**2*z**4*dble(3040) + CA**2*z*dble(3704) + CA**2*dble(-648)*Li3(-z) + 
     -       CA**2*z*dble(-648)*Li3(-z) + CA**2*z**2*dble(-648)*Li3(-z) + CA**2*z**3*dble(648)*Li3(-z) + 
     -       CA**2*z**4*dble(648)*Li3(-z) + CA**2*z**5*dble(648)*Li3(-z) + CA**2*dble(-1080)*Li3(z) + 
     -       CA**2*z**2*dble(-1080)*Li3(z) + CA**2*z**3*dble(-216)*Li3(z) + CA**2*z**5*dble(-216)*Li3(z) + 
     -       CA**2*z*dble(216)*Li3(z) + CA**2*z**4*dble(1080)*Li3(z) + CA**2*dble(-432)*Li3(1/(z + dble(1))) + 
     -       CA**2*z*dble(-432)*Li3(1/(z + dble(1))) + CA**2*z**2*dble(-432)*Li3(1/(z + dble(1))) + 
     -       CA**2*z**3*dble(432)*Li3(1/(z + dble(1))) + CA**2*z**4*dble(432)*Li3(1/(z + dble(1))) + 
     -       CA**2*z**5*dble(432)*Li3(1/(z + dble(1))) + CA**2*z*dble(-2103)*log(z) + CF*nf*TF*z**3*dble(-648)*log(z) + 
     -       CF*nf*TF*z**4*dble(-648)*log(z) + CA**2*z**2*dble(-447)*log(z) + CA*nf*TF*z**3*dble(-156)*log(z) + 
     -       CA*nf*TF*z**4*dble(-120)*log(z) + CA*nf*TF*z**2*dble(120)*log(z) + CA*nf*TF*z*dble(156)*log(z) + 
     -       CA**2*z**4*dble(447)*log(z) + CA**2*z**3*dble(495)*log(z) + CF*nf*TF*z*dble(648)*log(z) + 
     -       CF*nf*TF*z**2*dble(648)*log(z) + CA**2*z**5*dble(1608)*log(z) + 
     -       CA**2*dble(-216)*(z + dble(-1))*(z + z**2 + dble(1))**2*Li2(-z)*log(z) + CA**2*z**5*dble(-396)*log(z)**2 + 
     -       CF*nf*TF*z**3*dble(-162)*log(z)**2 + CA**2*z**2*dble(-99)*log(z)**2 + CF*nf*TF*z**4*dble(-54)*log(z)**2 + 
     -       CA*nf*TF*z**3*dble(-36)*log(z)**2 + CA*nf*TF*z**4*dble(-36)*log(z)**2 + CA*nf*TF*z*dble(36)*log(z)**2 + 
     -       CA*nf*TF*z**2*dble(36)*log(z)**2 + CF*nf*TF*z**2*dble(54)*log(z)**2 + CA**2*z**4*dble(99)*log(z)**2 + 
     -       CF*nf*TF*z*dble(162)*log(z)**2 + CA**2*z**3*dble(171)*log(z)**2 + CA**2*z*dble(225)*log(z)**2 + 
     -       CA**2*z**2*dble(-72)*log(z)**3 + CA**2*z*dble(-36)*log(z)**3 + CF*nf*TF*z**3*dble(-36)*log(z)**3 + 
     -       CF*nf*TF*z**4*dble(-36)*log(z)**3 + CA**2*z**5*dble(-36)*log(z)**3 + CF*nf*TF*z*dble(36)*log(z)**3 + 
     -       CF*nf*TF*z**2*dble(36)*log(z)**3 + CA**2*z**3*dble(36)*log(z)**3 + CA**2*z**4*dble(72)*log(z)**3 + 
     -       CA**2*dble(72)*Li2(z)*(dble(-1)*(z + dble(-1))**2*(z*dble(10) + z**2*dble(10) + dble(11) + z**3*dble(11)) + 
     -          dble(3)*(-z + z**3 + z**5 + z**4*dble(-3) + dble(3) + z**2*dble(3))*log(z)) + 
     -       CA*nf*TF*z**2*dble(-36)*log(-z + dble(1)) + CA**2*z**4*dble(-18)*log(-z + dble(1)) + 
     -       CA**2*z**2*dble(18)*log(-z + dble(1)) + CA*nf*TF*z**4*dble(36)*log(-z + dble(1)) + 
     -       CA**2*dble(-792)*log(z)*log(-z + dble(1)) + CA**2*z**5*dble(-792)*log(z)*log(-z + dble(1)) + 
     -       CA**2*z**2*dble(-72)*log(z)*log(-z + dble(1)) + CA**2*z**3*dble(-72)*log(z)*log(-z + dble(1)) + 
     -       CA**2*z*dble(864)*log(z)*log(-z + dble(1)) + CA**2*z**4*dble(864)*log(z)*log(-z + dble(1)) + 
     -       CA**2*z*dble(-108)*log(z)**2*log(-z + dble(1)) + CA**2*z**4*dble(-108)*log(z)**2*log(-z + dble(1)) + 
     -       CA**2*dble(108)*log(z)**2*log(-z + dble(1)) + CA**2*z**2*dble(108)*log(z)**2*log(-z + dble(1)) + 
     -       CA**2*z**3*dble(108)*log(z)**2*log(-z + dble(1)) + CA**2*z**5*dble(108)*log(z)**2*log(-z + dble(1)) + 
     -       CA**2*z*dble(-108)*log(z)*log(-z + dble(1))**2 + CA**2*z**4*dble(-108)*log(z)*log(-z + dble(1))**2 + 
     -       CA**2*dble(108)*log(z)*log(-z + dble(1))**2 + CA**2*z**2*dble(108)*log(z)*log(-z + dble(1))**2 + 
     -       CA**2*z**3*dble(108)*log(z)*log(-z + dble(1))**2 + CA**2*z**5*dble(108)*log(z)*log(-z + dble(1))**2 + 
     -       CA**2*Zeta2*dble(-216)*log(z + dble(1)) + CA**2*z*Zeta2*dble(-216)*log(z + dble(1)) + 
     -       CA**2*z**2*Zeta2*dble(-216)*log(z + dble(1)) + CA**2*z**3*Zeta2*dble(216)*log(z + dble(1)) + 
     -       CA**2*z**4*Zeta2*dble(216)*log(z + dble(1)) + CA**2*z**5*Zeta2*dble(216)*log(z + dble(1)) + 
     -       CA**2*dble(-108)*log(z)**2*log(z + dble(1)) + CA**2*z*dble(-108)*log(z)**2*log(z + dble(1)) + 
     -       CA**2*z**2*dble(-108)*log(z)**2*log(z + dble(1)) + CA**2*z**3*dble(108)*log(z)**2*log(z + dble(1)) + 
     -       CA**2*z**4*dble(108)*log(z)**2*log(z + dble(1)) + CA**2*z**5*dble(108)*log(z)**2*log(z + dble(1)) + 
     -       CA**2*z**3*dble(-72)*log(z + dble(1))**3 + CA**2*z**4*dble(-72)*log(z + dble(1))**3 + 
     -       CA**2*z**5*dble(-72)*log(z + dble(1))**3 + CA**2*dble(72)*log(z + dble(1))**3 + CA**2*z*dble(72)*log(z + dble(1))**3 + 
     -       CA**2*z**2*dble(72)*log(z + dble(1))**3)*PlusDist(-2,z,x,f,f0))/(z*(z**2 + dble(-1))) + 
     -  CA*(nf*TF*dble(8)*(dble(-328) + Zeta3*dble(252) + Zeta2*dble(270)) + 
     -     CA*(Zeta2*dble(-7236) + Zeta3*dble(-5544) + Zeta4*dble(4050) + dble(9712)))*dble(dble(1)/dble(648))*PlusDist(-1,z,x,f,f0)
     -    + CA*(nf*TF*dble(56) + CA*(dble(-202) + Zeta3*dble(189)))*dble(dble(4)/dble(27))*PlusDist(0,z,x,f,f0)
        return
        end

        double precision function YLI2gq(z,x,f,f0)
        implicit none
        double precision z,x,f,f0
        include 'paramqt.f'
        YLI2gq = (CF*dble(dble(1)/dble(54))*(CA*dble(-6320) + CA*z*Zeta3*dble(-2160) + CA*z**2*dble(-2144) + CA*z*Zeta2*dble(-1728) + 
     -      nf*TF*z*dble(-896) + CA*z**3*Zeta2*dble(-288) + CF*z**2*dble(-54) + nf*TF*z**2*dble(208) + CA*z**2*Zeta2*dble(324) + 
     -      CF*z*dble(540) + nf*TF*dble(896) + CA*z**3*dble(1216) + CA*Zeta2*dble(1584) + CA*z**2*Zeta3*dble(1728) + 
     -      CA*Zeta3*dble(3456) + CA*z*dble(6328) + CA*dble(-1296)*Li3(-z) + CA*z*dble(-1296)*Li3(-z) + 
     -      CA*z**2*dble(-648)*Li3(-z) + CA*dble(-2160)*Li3(z) + CA*z**2*dble(-1080)*Li3(z) + CA*z*dble(432)*Li3(z) + 
     -      CA*dble(-864)*Li3(1/(z + dble(1))) + CA*z*dble(-864)*Li3(1/(z + dble(1))) + CA*z**2*dble(-432)*Li3(1/(z + dble(1))) + 
     -      CA*z*dble(-2988)*log(z) + CA*z**3*dble(-1056)*log(z) + CF*z*dble(-810)*log(z) + CA*z**2*dble(72)*log(z) + 
     -      CF*z**2*dble(270)*log(z) + CF*z*dble(-108)*log(z)**2 + CF*z**2*dble(-81)*log(z)**2 + CA*z**3*dble(144)*log(z)**2 + 
     -      CA*z**2*dble(162)*log(z)**2 + CA*z*dble(648)*log(z)**2 + CA*z*dble(-72)*log(z)**3 + CA*z**2*dble(-36)*log(z)**3 + 
     -      CF*z**2*dble(-18)*log(z)**3 + CF*z*dble(36)*log(z)**3 + 
     -      CA*dble(216)*Li2(-z)*(z**2 + (z**2 + dble(2) + z*dble(2))*log(z)) + 
     -      CA*dble(72)*Li2(z)*(dble(-22) + z**2*dble(-6) + z**3*dble(4) + z*dble(24) + 
     -         dble(3)*(z*dble(-2) + z**2*dble(3) + dble(6))*log(z)) + CA*dble(-1824)*log(-z + dble(1)) + 
     -      CF*z*dble(-1728)*log(-z + dble(1)) + CA*z**2*dble(-516)*log(-z + dble(1)) + nf*TF*z*dble(-480)*log(-z + dble(1)) + 
     -      nf*TF*z**2*dble(96)*log(-z + dble(1)) + nf*TF*dble(480)*log(-z + dble(1)) + CF*z**2*dble(540)*log(-z + dble(1)) + 
     -      CF*dble(1728)*log(-z + dble(1)) + CA*z*dble(1824)*log(-z + dble(1)) + CA*dble(-1584)*log(z)*log(-z + dble(1)) + 
     -      CA*z**2*dble(-648)*log(z)*log(-z + dble(1)) + CA*z**3*dble(288)*log(z)*log(-z + dble(1)) + 
     -      CA*z*dble(1728)*log(z)*log(-z + dble(1)) + CA*z*dble(-216)*log(z)**2*log(-z + dble(1)) + 
     -      CA*z**2*dble(108)*log(z)**2*log(-z + dble(1)) + CA*dble(216)*log(z)**2*log(-z + dble(1)) + 
     -      CA*dble(-396)*log(-z + dble(1))**2 + CF*z*dble(-324)*log(-z + dble(1))**2 + nf*TF*z*dble(-144)*log(-z + dble(1))**2 + 
     -      CA*z**2*dble(-90)*log(-z + dble(1))**2 + CF*z**2*dble(54)*log(-z + dble(1))**2 + 
     -      nf*TF*z**2*dble(72)*log(-z + dble(1))**2 + nf*TF*dble(144)*log(-z + dble(1))**2 + CF*dble(324)*log(-z + dble(1))**2 + 
     -      CA*z*dble(396)*log(-z + dble(1))**2 + CA*z*dble(-216)*log(z)*log(-z + dble(1))**2 + 
     -      CA*z**2*dble(108)*log(z)*log(-z + dble(1))**2 + CA*dble(216)*log(z)*log(-z + dble(1))**2 + 
     -      CA*dble(-72)*log(-z + dble(1))**3 + CF*z*dble(-72)*log(-z + dble(1))**3 + CA*z**2*dble(-36)*log(-z + dble(1))**3 + 
     -      CF*z**2*dble(36)*log(-z + dble(1))**3 + CF*dble(72)*log(-z + dble(1))**3 + CA*z*dble(72)*log(-z + dble(1))**3 + 
     -      CA*Zeta2*dble(-432)*log(z + dble(1)) + CA*z*Zeta2*dble(-432)*log(z + dble(1)) + 
     -      CA*z**2*Zeta2*dble(-216)*log(z + dble(1)) + CA*z**2*dble(216)*log(z)*log(z + dble(1)) + 
     -      CA*dble(-216)*log(z)**2*log(z + dble(1)) + CA*z*dble(-216)*log(z)**2*log(z + dble(1)) + 
     -      CA*z**2*dble(-108)*log(z)**2*log(z + dble(1)) + CA*z**2*dble(72)*log(z + dble(1))**3 + 
     -      CA*dble(144)*log(z + dble(1))**3 + CA*z*dble(144)*log(z + dble(1))**3)*PlusDist(-2,z,x,f,f0))/z
        return
        end

!===============================================

        double precision function YLP0qiqiYLP0qiqi(z,x,f,f0)
        implicit none
        double precision z,x,f,f0
        include 'paramqt.f'
        YLP0qiqiYLP0qiqi = (CF**2*(z*dble(-4) + z**2*dble(-1) + dble(5) + log(z) + z**2*dble(3)*log(z) + dble(-4)*(z**2 + dble(-1))*log(-z + dble(1)))*
     -     PlusDist(-2,z,x,f,f0))/(z + dble(-1)) + CF**2*(Zeta2*dble(-48) + dble(27))*dble(dble(1)/dble(12))*
     -   PlusDist(-1,z,x,f,f0) + CF**2*dble(6)*PlusDist(0,z,x,f,f0) + CF**2*dble(8)*PlusDist(1,z,x,f,f0)
        return
        end

        double precision function YLP0qgYLP0gq(z,x,f,f0)
        implicit none
        double precision z,x,f,f0
        include 'paramqt.f'
        YLP0qgYLP0gq = (CF*TF*dble(dble(1)/dble(3))*(z**3*dble(-4) + z**2*dble(-3) + z*dble(3) + dble(4) + z*(z + dble(1))*dble(6)*log(z))*
     -    PlusDist(-2,z,x,f,f0))/z
        return
        end

        double precision function YLP0qgYLP0gg(z,x,f,f0)
        implicit none
        double precision z,x,f,f0
        include 'paramqt.f'
        YLP0qgYLP0gg = (TF*dble(dble(1)/dble(6))*(CA*z**3*dble(-62) + beta0*z**2*dble(-6) + beta0*z*dble(3) + CA*z*dble(6) + beta0*z**3*dble(6) + 
     -      CA*dble(8) + CA*z**2*dble(48) + CA*z*(dble(1) + z*dble(4))*dble(12)*log(z) + 
     -      CA*z*(z*dble(-2) + dble(1) + z**2*dble(2))*dble(12)*log(-z + dble(1)))*PlusDist(-2,z,x,f,f0))/z
        return
        end

        double precision function YLP0qiqiYLP0qg(z,x,f,f0)
        implicit none
        double precision z,x,f,f0
        include 'paramqt.f'
        YLP0qiqiYLP0qg = CF*TF*dble(dble(1)/dble(2))*(dble(-1) + z*dble(4) + (z**2*dble(-8) + dble(-2) + z*dble(4))*log(z) + 
     -    (z*dble(-8) + dble(4) + z**2*dble(8))*log(-z + dble(1)))*PlusDist(-2,z,x,f,f0)
        return
        end

        double precision function YLP0ggYLP0gg(z,x,f,f0)
        implicit none
        double precision z,x,f,f0
        include 'paramqt.f'
        YLP0ggYLP0gg = (CA*dble(dble(-2)/dble(3))*(dble(-1)*(z + dble(-1))*
     -        (beta0*(z**2 + z*dble(-2) + z**3*dble(-1) + dble(1))*dble(3) + 
     -          CA*dble(2)*(dble(-11) + z**2*dble(-9) + z*dble(9) + z**3*dble(11))) + 
     -       CA*dble(-6)*(z**4 + z**3*dble(-4) + dble(1) + z**2*dble(3))*log(z) + 
     -       CA*(z**4 + z*dble(-3) + z**3*dble(-2) + dble(1) + z**2*dble(3))*dble(12)*log(-z + dble(1)))*PlusDist(-2,z,x,f,f0))/
     -   (z*(z + dble(-1))) + (CA**2*Zeta2*dble(-4) + beta0**2*dble(dble(1)/dble(4)))*PlusDist(-1,z,x,f,f0) + 
     -  beta0*CA*dble(2)*PlusDist(0,z,x,f,f0) + CA**2*dble(8)*PlusDist(1,z,x,f,f0)
        return
        end

        double precision function YLP0gqYLP0qg(z,x,f,f0)
        implicit none
        double precision z,x,f,f0
        include 'paramqt.f'
        double precision YLP0qgYLP0gq
        YLP0gqYLP0qg = YLP0qgYLP0gq(z,x,f,f0)
        return
        end

        double precision function YLP0gqYLP0qiqi(z,x,f,f0)
        implicit none
        double precision z,x,f,f0
        include 'paramqt.f'
        YLP0gqYLP0qiqi = (CF**2*dble(dble(1)/dble(2))*(-(z*(z + dble(-4) + (z + dble(-2))*dble(2)*log(z))) + 
     -      (z**2 + z*dble(-2) + dble(2))*dble(4)*log(-z + dble(1)))*PlusDist(-2,z,x,f,f0))/z
        return
        end

        double precision function YLP0ggYLP0gq(z,x,f,f0)
        implicit none
        double precision z,x,f,f0
        include 'paramqt.f'
        YLP0ggYLP0gq = (CF*dble(dble(1)/dble(6))*(CA*dble(-62) + beta0*z*dble(-6) + beta0*z**2*dble(3) + beta0*dble(6) + CA*z**2*dble(6) + 
     -      CA*z**3*dble(8) + CA*z*dble(48) + CA*dble(-24)*(z + z**2 + dble(1))*log(z) + 
     -      CA*(z**2 + z*dble(-2) + dble(2))*dble(12)*log(-z + dble(1)))*PlusDist(-2,z,x,f,f0))/z
        return
        end

!===============================================

        double precision function YLI1qiqiYLP0qiqi(z,x,f,f0)
        implicit none
        double precision z,x,f,f0
        include 'paramqt.f'
        YLI1qiqiYLP0qiqi = CF**2*dble(dble(1)/dble(6))*(dble(-6) + z*dble(6) + Zeta2*dble(6) + z*Zeta2*dble(6) + (z + dble(-1))*dble(12)*log(z) + 
     -     dble(-24)*(z + dble(-1))*log(-z + dble(1)))*PlusDist(-2,z,x,f,f0) + 
     -  CF**2*Zeta2*dble(dble(-3)/dble(2))*PlusDist(-1,z,x,f,f0) + CF**2*Zeta2*dble(-2)*PlusDist(0,z,x,f,f0)
        return
        end

        double precision function YLI1qgYLP0gq(z,x,f,f0)
        implicit none
        double precision z,x,f,f0
        include 'paramqt.f'
        YLI1qgYLP0gq = (CF*TF*dble(dble(4)/dble(3))*(z*dble(-3) + dble(1) + z**3*dble(2) + z**2*dble(-3)*log(z))*PlusDist(-2,z,x,f,f0))/z
        return
        end

        double precision function YLI1qgYLP0gg(z,x,f,f0)
        implicit none
        double precision z,x,f,f0
        include 'paramqt.f'
        YLI1qgYLP0gg = (TF*dble(dble(-2)/dble(3))*(dble(-1)*(z + dble(-1))*(beta0*z**2*dble(-3) + CA*(dble(-2) + z*dble(4) + z**2*dble(34))) + 
     -      CA*z**2*dble(24)*log(z) + CA*z**2*(z + dble(-1))*dble(12)*log(-z + dble(1)))*PlusDist(-2,z,x,f,f0))/z
        return
        end

        double precision function YLI1qiqiYLP0qg(z,x,f,f0)
        implicit none
        double precision z,x,f,f0
        include 'paramqt.f'
        YLI1qiqiYLP0qg = CF*TF*dble(dble(-1)/dble(6))*(dble(-12)*(z + z**2 + dble(-2)) + Zeta2*(z*dble(-2) + dble(1) + z**2*dble(2))*dble(6) + 
     -    (dble(1) + z*dble(2))*dble(12)*log(z))*PlusDist(-2,z,x,f,f0)
        return
        end

        double precision function YLI1ggYLP0gg(z,x,f,f0)
        implicit none
        double precision z,x,f,f0
        include 'paramqt.f'
        YLI1ggYLP0gg = (CA**2*Zeta2*dble(2)*(z**3 + dble(-1) + z**2*dble(-1) + z*dble(2))*PlusDist(-2,z,x,f,f0))/z + 
     -  beta0*CA*Zeta2*dble(dble(-1)/dble(2))*PlusDist(-1,z,x,f,f0) + CA**2*Zeta2*dble(-2)*PlusDist(0,z,x,f,f0)
        return
        end

        double precision function YLI1gqYLP0qg(z,x,f,f0)
        implicit none
        double precision z,x,f,f0
        include 'paramqt.f'
        YLI1gqYLP0qg = CF*TF*dble(2)*(z + z**2*dble(-2) + dble(1) + z*dble(2)*log(z))*PlusDist(-2,z,x,f,f0)
        return
        end

        double precision function YLI1gqYLP0qiqi(z,x,f,f0)
        implicit none
        double precision z,x,f,f0
        include 'paramqt.f'
        YLI1gqYLP0qiqi = CF**2*(z + dble(2) + z*dble(-2)*log(z) + z*dble(4)*log(-z + dble(1)))*PlusDist(-2,z,x,f,f0)
        return
        end

        double precision function YLI1ggYLP0gq(z,x,f,f0)
        implicit none
        double precision z,x,f,f0
        include 'paramqt.f'
        YLI1ggYLP0gq = -((CA*CF*Zeta2*(z**2 + z*dble(-2) + dble(2))*PlusDist(-2,z,x,f,f0))/z)
        return
        end

!=================================================================================================

        double precision function Cqq1qg(i1,i2,x1,x2,z1,z2,L,Lmu) ! include both qg->qq and gq->qq
        implicit none
        integer i1,i2
        double precision x1,x2,z1,z2,L,Lmu
        include 'paramqt.f'
        double precision YLP0qiqi,YLP0qg,YLI1qiqi,YLI1qg
        double precision fq,fq0,fg,fg0
        Cqq1qg = 0d0
        ! qg
        fq=PDFbuf1(i1,1)
        fq0=PDFbuf1(i1,0)
        fg=PDFbuf2(0,1)
        fg0=PDFbuf2(0,0)
        Cqq1qg = Cqq1qg + (YLI1qg(z2,x2,fg,fg0) + (L + Lmu)*dble(2)*YLP0qg(z2,x2,fg,fg0))*PlusDist(-1,z1,x1,fq,fq0)
        ! gq
        fq=PDFbuf2(i2,1)
        fq0=PDFbuf2(i2,0)
        fg=PDFbuf1(0,1)
        fg0=PDFbuf1(0,0)
        Cqq1qg = Cqq1qg + (YLI1qg(z1,x1,fg,fg0) + (L + Lmu)*dble(2)*YLP0qg(z1,x1,fg,fg0))*PlusDist(-1,z2,x2,fq,fq0)
        return
        end

        double precision function Cqq1qiqi(i1,i2,x1,x2,z1,z2,L,Lmu)
        implicit none
        integer i1,i2
        double precision x1,x2,z1,z2,L,Lmu
        include 'paramqt.f'
        double precision YLP0qiqi,YLP0qg,YLI1qiqi,YLI1qg
        double precision fqa,fqa0,fqb,fqb0
        fqa=PDFbuf1(i1,1)
        fqa0=PDFbuf1(i1,0)
        fqb=PDFbuf2(i2,1)
        fqb0=PDFbuf2(i2,0)
        Cqq1qiqi = YLI1qiqi(z2,x2,fqb,fqb0)*PlusDist(-1,z1,x1,fqa,fqa0) +
     -  (L + Lmu)*dble(2)*YLP0qiqi(z2,x2,fqb,fqb0)*PlusDist(-1,z1,x1,fqa,fqa0) +
     -  dble(dble(1)/dble(2))*(dble(2)*YLI1qiqi(z1,x1,fqa,fqa0) + (L + Lmu)*dble(4)*YLP0qiqi(z1,x1,fqa,fqa0) +
     -     (-(CF*G0*L**2) + gH0qq2*Lmu*dble(-2) + CF*d1*L*dble(2) + gH0qq*L*dble(2) + gH0qq*Lmu*dble(2) + cH1qq*dble(4))*
     -      PlusDist(-1,z1,x1,fqa,fqa0))*PlusDist(-1,z2,x2,fqb,fqb0)
        return
        end

        double precision function Cqq2qiqj(i1,i2,x1,x2,z1,z2,L,Lmu) ! include both qiqj->qiqi and qjqi->qiqi
        implicit none
        integer i1,i2
        double precision x1,x2,z1,z2,L,Lmu
        include 'paramqt.f'
        double precision YLP0qiqi,YLP0qg,YLP1qiqi,YLP1qiqj,YLP1qiqbi,YLP1qg
        double precision YLI1qiqi,YLI1qg,YLI2qiqi,YLI2qiqj,YLI2qiqbi,YLI2qg
        double precision YLP0qiqiYLP0qiqi,YLP0qgYLP0gq,YLP0qgYLP0gg,YLP0qiqiYLP0qg
        double precision YLI1qiqiYLP0qiqi,YLI1qgYLP0gq,YLI1qgYLP0gg,YLI1qiqiYLP0qg
        double precision fq,fq0,fqj,fqj0
        integer id
        Cqq2qiqj = 0d0
        ! qiqj
        fq=PDFbuf1(i1,1)
        fq0=PDFbuf1(i1,0)
        fqj=0d0
        fqj0=0d0
        id=int(nf)
        do while (id.gt.0)
                if ((id.ne.i2).and.(id.ne.-i2)) then
                        fqj=fqj+PDFbuf2(id,1)+PDFbuf2(-id,1)
                        fqj0=fqj0+PDFbuf2(id,0)+PDFbuf2(-id,0)
                end if
                id=id-1
        end do
        Cqq2qiqj = Cqq2qiqj + (YLI2qiqj(z2,x2,fqj,fqj0) + (L + Lmu)**2*dble(2)*YLP0qgYLP0gq(z2,x2,fqj,fqj0) + 
     -    (L + Lmu)*dble(2)*(YLI1qgYLP0gq(z2,x2,fqj,fqj0) + dble(2)*YLP1qiqj(z2,x2,fqj,fqj0)))*PlusDist(-1,z1,x1,fq,fq0)
        ! qjqi
        fq=PDFbuf2(i2,1)
        fq0=PDFbuf2(i2,0)
        fqj=0d0
        fqj0=0d0
        id=int(nf)
        do while (id.gt.0)
                if ((id.ne.i1).and.(id.ne.-i1)) then
                        fqj=fqj+PDFbuf1(id,1)+PDFbuf1(-id,1)
                        fqj0=fqj0+PDFbuf1(id,0)+PDFbuf1(-id,0)
                end if
                id=id-1
        end do
        Cqq2qiqj = Cqq2qiqj + (YLI2qiqj(z1,x1,fqj,fqj0) + (L + Lmu)**2*dble(2)*YLP0qgYLP0gq(z1,x1,fqj,fqj0) +
     -    (L + Lmu)*dble(2)*(YLI1qgYLP0gq(z1,x1,fqj,fqj0) + dble(2)*YLP1qiqj(z1,x1,fqj,fqj0)))*PlusDist(-1,z2,x2,fq,fq0)
        return
        end

        double precision function Cqq2qiqbi(i1,i2,x1,x2,z1,z2,L,Lmu) ! include both qiqbi->qiqi and qbiqi->qiqi
        implicit none
        integer i1,i2
        double precision x1,x2,z1,z2,L,Lmu
        include 'paramqt.f'
        double precision YLP0qiqi,YLP0qg,YLP1qiqi,YLP1qiqj,YLP1qiqbi,YLP1qg
        double precision YLI1qiqi,YLI1qg,YLI2qiqi,YLI2qiqj,YLI2qiqbi,YLI2qg
        double precision YLP0qiqiYLP0qiqi,YLP0qgYLP0gq,YLP0qgYLP0gg,YLP0qiqiYLP0qg
        double precision YLI1qiqiYLP0qiqi,YLI1qgYLP0gq,YLI1qgYLP0gg,YLI1qiqiYLP0qg
        double precision fq,fq0,fqb,fqb0
        integer id
        Cqq2qiqbi = 0d0
        ! qiqbi
        fq=PDFbuf1(i1,1)
        fq0=PDFbuf1(i1,0)
        fqb=PDFbuf2(-i2,1)
        fqb0=PDFbuf2(-i2,0)
        Cqq2qiqbi = Cqq2qiqbi + (YLI2qiqbi(z2,x2,fqb,fqb0) + (L + Lmu)**2*dble(2)*YLP0qgYLP0gq(z2,x2,fqb,fqb0) + 
     -    (L + Lmu)*dble(2)*(YLI1qgYLP0gq(z2,x2,fqb,fqb0) + dble(2)*YLP1qiqbi(z2,x2,fqb,fqb0)))*PlusDist(-1,z1,x1,fq,fq0)
        ! qbiqi
        fq=PDFbuf2(i2,1)
        fq0=PDFbuf2(i2,0)
        fqb=PDFbuf1(-i1,1)
        fqb0=PDFbuf1(-i1,0)
        Cqq2qiqbi = Cqq2qiqbi + (YLI2qiqbi(z1,x1,fqb,fqb0) + (L + Lmu)**2*dble(2)*YLP0qgYLP0gq(z1,x1,fqb,fqb0) +
     -    (L + Lmu)*dble(2)*(YLI1qgYLP0gq(z1,x1,fqb,fqb0) + dble(2)*YLP1qiqbi(z1,x1,fqb,fqb0)))*PlusDist(-1,z2,x2,fq,fq0)
        return
        end

        double precision function Cqq2qg(i1,i2,x1,x2,z1,z2,L,Lmu) ! include both qg->qq and gq->qq
        implicit none
        integer i1,i2
        double precision x1,x2,z1,z2,L,Lmu
        include 'paramqt.f'
        double precision YLP0qiqi,YLP0qg,YLP1qiqi,YLP1qiqj,YLP1qiqbi,YLP1qg
        double precision YLI1qiqi,YLI1qg,YLI2qiqi,YLI2qiqj,YLI2qiqbi,YLI2qg
        double precision YLP0qiqiYLP0qiqi,YLP0qgYLP0gq,YLP0qgYLP0gg,YLP0qiqiYLP0qg
        double precision YLI1qiqiYLP0qiqi,YLI1qgYLP0gq,YLI1qgYLP0gg,YLI1qiqiYLP0qg
        double precision fq,fq0,fg,fg0
        Cqq2qg = 0d0
        ! qg
        fq=PDFbuf1(i1,1)
        fq0=PDFbuf1(i1,0)
        fg=PDFbuf2(0,1)
        fg0=PDFbuf2(0,0)
        Cqq2qg = Cqq2qg + YLI1qg(z2,x2,fg,fg0)*YLI1qiqi(z1,x1,fq,fq0) - CF*d1*Lmu*YLI1qg(z2,x2,fg,fg0)*PlusDist(-1,z1,x1,fq,fq0) +
     -  (-(Lmu*(CF*G0*Lmu + gH0qq2*dble(2))) + cH1qq*dble(4))*dble(dble(1)/dble(2))*YLI1qg(z2,x2,fg,fg0)*
     -   PlusDist(-1,z1,x1,fq,fq0) + YLI2qg(z2,x2,fg,fg0)*PlusDist(-1,z1,x1,fq,fq0) + 
     -  CF*G0*(Zeta3*dble(-4) + (L + Lmu)**3*dble(-1))*YLP0qg(z2,x2,fg,fg0)*PlusDist(-1,z1,x1,fq,fq0) +
     -  (L + Lmu)**2*dble(dble(1)/dble(2))*((-(CF*G0*YLI1qg(z2,x2,fg,fg0)) +
     -        dble(4)*(YLP0qgYLP0gg(z2,x2,fg,fg0) + YLP0qiqiYLP0qg(z2,x2,fg,fg0)))*PlusDist(-1,z1,x1,fq,fq0) +
     -     YLP0qg(z2,x2,fg,fg0)*(dble(8)*YLP0qiqi(z1,x1,fq,fq0) +
     -        dble(2)*(-beta0 + gH0qq*dble(2) + CF*(d1 + G0*Lmu)*dble(2))*PlusDist(-1,z1,x1,fq,fq0))) + 
     -  (L + Lmu)*dble(dble(-1)/dble(2))*(dble(-4)*YLI1qiqi(z1,x1,fq,fq0)*YLP0qg(z2,x2,fg,fg0) +
     -     dble(-2)*((dble(2)*YLI1qgYLP0gg(z2,x2,fg,fg0) + dble(2)*YLI1qiqiYLP0qg(z2,x2,fg,fg0) - CF*G0*Lmu**2*YLP0qg(z2,x2,fg,fg0) +
     -           CF*d1*Lmu*dble(-2)*YLP0qg(z2,x2,fg,fg0) + gH0qq2*Lmu*dble(-2)*YLP0qg(z2,x2,fg,fg0) +
     -           cH1qq*dble(4)*YLP0qg(z2,x2,fg,fg0) + dble(4)*YLP1qg(z2,x2,fg,fg0))*PlusDist(-1,z1,x1,fq,fq0) +
     -        YLI1qg(z2,x2,fg,fg0)*(dble(2)*YLP0qiqi(z1,x1,fq,fq0) + (-beta0 + CF*d1 + gH0qq + CF*G0*Lmu)*PlusDist(-1,z1,x1,fq,fq0))))
        ! gq
        fq=PDFbuf2(i2,1)
        fq0=PDFbuf2(i2,0)
        fg=PDFbuf1(0,1)
        fg0=PDFbuf1(0,0)
        Cqq2qg = Cqq2qg + YLI1qg(z1,x1,fg,fg0)*YLI1qiqi(z2,x2,fq,fq0) - CF*d1*Lmu*YLI1qg(z1,x1,fg,fg0)*PlusDist(-1,z2,x2,fq,fq0) +
     -  (-(Lmu*(CF*G0*Lmu + gH0qq2*dble(2))) + cH1qq*dble(4))*dble(dble(1)/dble(2))*YLI1qg(z1,x1,fg,fg0)*
     -   PlusDist(-1,z2,x2,fq,fq0) + YLI2qg(z1,x1,fg,fg0)*PlusDist(-1,z2,x2,fq,fq0) + 
     -  CF*G0*(Zeta3*dble(-4) + (L + Lmu)**3*dble(-1))*YLP0qg(z1,x1,fg,fg0)*PlusDist(-1,z2,x2,fq,fq0) +
     -  (L + Lmu)**2*dble(dble(1)/dble(2))*((-(CF*G0*YLI1qg(z1,x1,fg,fg0)) +
     -        dble(4)*(YLP0qgYLP0gg(z1,x1,fg,fg0) + YLP0qiqiYLP0qg(z1,x1,fg,fg0)))*PlusDist(-1,z2,x2,fq,fq0) +
     -     YLP0qg(z1,x1,fg,fg0)*(dble(8)*YLP0qiqi(z2,x2,fq,fq0) +
     -        dble(2)*(-beta0 + gH0qq*dble(2) + CF*(d1 + G0*Lmu)*dble(2))*PlusDist(-1,z2,x2,fq,fq0))) + 
     -  (L + Lmu)*dble(dble(-1)/dble(2))*(dble(-4)*YLI1qiqi(z2,x2,fq,fq0)*YLP0qg(z1,x1,fg,fg0) +
     -     dble(-2)*((dble(2)*YLI1qgYLP0gg(z1,x1,fg,fg0) + dble(2)*YLI1qiqiYLP0qg(z1,x1,fg,fg0) - CF*G0*Lmu**2*YLP0qg(z1,x1,fg,fg0) +
     -           CF*d1*Lmu*dble(-2)*YLP0qg(z1,x1,fg,fg0) + gH0qq2*Lmu*dble(-2)*YLP0qg(z1,x1,fg,fg0) +
     -           cH1qq*dble(4)*YLP0qg(z1,x1,fg,fg0) + dble(4)*YLP1qg(z1,x1,fg,fg0))*PlusDist(-1,z2,x2,fq,fq0) +
     -        YLI1qg(z1,x1,fg,fg0)*(dble(2)*YLP0qiqi(z2,x2,fq,fq0) + (-beta0 + CF*d1 + gH0qq + CF*G0*Lmu)*PlusDist(-1,z2,x2,fq,fq0))))
        return
        end

        double precision function Cqq2gg(i1,i2,x1,x2,z1,z2,L,Lmu)
        implicit none
        integer i1,i2
        double precision x1,x2,z1,z2,L,Lmu
        include 'paramqt.f'
        double precision YLP0qiqi,YLP0qg,YLP1qiqi,YLP1qiqj,YLP1qiqbi,YLP1qg
        double precision YLI1qiqi,YLI1qg,YLI2qiqi,YLI2qiqj,YLI2qiqbi,YLI2qg
        double precision YLP0qiqiYLP0qiqi,YLP0qgYLP0gq,YLP0qgYLP0gg,YLP0qiqiYLP0qg
        double precision YLI1qiqiYLP0qiqi,YLI1qgYLP0gq,YLI1qgYLP0gg,YLI1qiqiYLP0qg
        double precision fga,fga0,fgb,fgb0
        fga=PDFbuf1(0,1)
        fga0=PDFbuf1(0,0)
        fgb=PDFbuf2(0,1)
        fgb0=PDFbuf2(0,0)
        Cqq2gg = (YLI1qg(z1,x1,fga,fga0) + (L + Lmu)*dble(2)*YLP0qg(z1,x1,fga,fga0))*(YLI1qg(z2,x2,fgb,fgb0) + (L + Lmu)*dble(2)*YLP0qg(z2,x2,fgb,fgb0)) 
        return
        end

        double precision function Cqq2qiqi(i1,i2,x1,x2,z1,z2,L,Lmu)
        implicit none
        integer i1,i2
        double precision x1,x2,z1,z2,L,Lmu
        include 'paramqt.f'
        double precision YLP0qiqi,YLP0qg,YLP1qiqi,YLP1qiqj,YLP1qiqbi,YLP1qg
        double precision YLI1qiqi,YLI1qg,YLI2qiqi,YLI2qiqj,YLI2qiqbi,YLI2qg
        double precision YLP0qiqiYLP0qiqi,YLP0qgYLP0gq,YLP0qgYLP0gg,YLP0qiqiYLP0qg
        double precision YLI1qiqiYLP0qiqi,YLI1qgYLP0gq,YLI1qgYLP0gg,YLI1qiqiYLP0qg
        double precision fqa,fqa0,fqb,fqb0
        fqa=PDFbuf1(i1,1)
        fqa0=PDFbuf1(i1,0)
        fqb=PDFbuf2(i2,1)
        fqb0=PDFbuf2(i2,0)
        Cqq2qiqi = YLI1qiqi(z1,x1,fqa,fqa0)*YLI1qiqi(z2,x2,fqb,fqb0) - CF*d1*Lmu*YLI1qiqi(z2,x2,fqb,fqb0)*PlusDist(-1,z1,x1,fqa,fqa0) + 
     -  (-(Lmu*(CF*G0*Lmu + gH0qq2*dble(2))) + cH1qq*dble(4))*dble(dble(1)/dble(2))*YLI1qiqi(z2,x2,fqb,fqb0)*
     -   PlusDist(-1,z1,x1,fqa,fqa0) + YLI2qiqi(z2,x2,fqb,fqb0)*PlusDist(-1,z1,x1,fqa,fqa0) - 
     -  CF*d1*Lmu*YLI1qiqi(z1,x1,fqa,fqa0)*PlusDist(-1,z2,x2,fqb,fqb0) + 
     -  (-(Lmu*(CF*G0*Lmu + gH0qq2*dble(2))) + cH1qq*dble(4))*dble(dble(1)/dble(2))*YLI1qiqi(z1,x1,fqa,fqa0)*
     -   PlusDist(-1,z2,x2,fqb,fqb0) + YLI2qiqi(z1,x1,fqa,fqa0)*PlusDist(-1,z2,x2,fqb,fqb0) - 
     -  CF*d2*Lmu*PlusDist(-1,z1,x1,fqa,fqa0)*PlusDist(-1,z2,x2,fqb,fqb0) + 
     -  cH1qq*(beta0 + gH0qq2)*Lmu*dble(-2)*PlusDist(-1,z1,x1,fqa,fqa0)*PlusDist(-1,z2,x2,fqb,fqb0) + 
     -  CF**2*d1**2*Lmu**2*dble(dble(1)/dble(2))*PlusDist(-1,z1,x1,fqa,fqa0)*PlusDist(-1,z2,x2,fqb,fqb0) + 
     -  Lmu**2*(gH0qq2*(beta0 + gH0qq2) + CF*cH1qq*G0*dble(-2))*dble(dble(1)/dble(2))*PlusDist(-1,z1,x1,fqa,fqa0)*
     -   PlusDist(-1,z2,x2,fqb,fqb0) + CF*d1*Lmu*(cH1qq*dble(-4) + Lmu*(CF*G0*Lmu + gH0qq2*dble(2)))*dble(dble(1)/dble(2))*
     -   PlusDist(-1,z1,x1,fqa,fqa0)*PlusDist(-1,z2,x2,fqb,fqb0) + 
     -  (-(Lmu*(CF*G1*Lmu + gH1qq2*dble(2))) + cH2qq*dble(4))*dble(dble(1)/dble(2))*PlusDist(-1,z1,x1,fqa,fqa0)*
     -   PlusDist(-1,z2,x2,fqb,fqb0) + CF*G0*Lmu**3*(beta0 + gH0qq2*dble(3))*dble(dble(1)/dble(6))*PlusDist(-1,z1,x1,fqa,fqa0)*
     -   PlusDist(-1,z2,x2,fqb,fqb0) + CF**2*G0**2*Lmu**4*dble(dble(1)/dble(8))*PlusDist(-1,z1,x1,fqa,fqa0)*
     -   PlusDist(-1,z2,x2,fqb,fqb0) + CF**2*G0**2*((L + Lmu)**4 + (L + Lmu)*Zeta3*dble(16))*dble(dble(1)/dble(8))*
     -   PlusDist(-1,z1,x1,fqa,fqa0)*PlusDist(-1,z2,x2,fqb,fqb0) + 
     -  CF*G0*(L**3 + Lmu**3 + L**2*Lmu*dble(3) + L*Lmu**2*dble(3) + Zeta3*dble(4))*dble(dble(-1)/dble(6))*
     -   (dble(6)*YLP0qiqi(z2,x2,fqb,fqb0)*PlusDist(-1,z1,x1,fqa,fqa0) + 
     -     (dble(6)*YLP0qiqi(z1,x1,fqa,fqa0) + (beta0*dble(-2) + gH0qq*dble(3) + CF*(d1 + G0*Lmu)*dble(3))*
     -         PlusDist(-1,z1,x1,fqa,fqa0))*PlusDist(-1,z2,x2,fqb,fqb0)) + 
     -  (L + Lmu)*dble(dble(-1)/dble(4))*(CF*d1*dble(-4)*YLI1qiqi(z2,x2,fqb,fqb0)*PlusDist(-1,z1,x1,fqa,fqa0) + 
     -     CF*G0*Lmu*dble(-4)*YLI1qiqi(z2,x2,fqb,fqb0)*PlusDist(-1,z1,x1,fqa,fqa0) + 
     -     dble(-2)*YLI1qiqi(z2,x2,fqb,fqb0)*(dble(4)*YLP0qiqi(z1,x1,fqa,fqa0) + gH0qq*PlusDist(-1,z1,x1,fqa,fqa0)) + 
     -     CF*d1*dble(-4)*YLI1qiqi(z1,x1,fqa,fqa0)*PlusDist(-1,z2,x2,fqb,fqb0) + 
     -     CF*G0*Lmu*dble(-4)*YLI1qiqi(z1,x1,fqa,fqa0)*PlusDist(-1,z2,x2,fqb,fqb0) + 
     -     CF*d2*dble(-4)*PlusDist(-1,z1,x1,fqa,fqa0)*PlusDist(-1,z2,x2,fqb,fqb0) + 
     -     beta0*CF*d1*Lmu*dble(-4)*PlusDist(-1,z1,x1,fqa,fqa0)*PlusDist(-1,z2,x2,fqb,fqb0) + 
     -     CF*G1*Lmu*dble(-4)*PlusDist(-1,z1,x1,fqa,fqa0)*PlusDist(-1,z2,x2,fqb,fqb0) + 
     -     CF*d1*dble(2)*(cH1qq*dble(-4) + Lmu*(CF*G0*Lmu + gH0qq2*dble(2)))*PlusDist(-1,z1,x1,fqa,fqa0)*
     -      PlusDist(-1,z2,x2,fqb,fqb0) + CF*G0*Lmu*dble(2)*(cH1qq*dble(-4) + Lmu*(CF*G0*Lmu + gH0qq2*dble(2)))*
     -      PlusDist(-1,z1,x1,fqa,fqa0)*PlusDist(-1,z2,x2,fqb,fqb0) + 
     -     CF**2*d1**2*Lmu*dble(4)*PlusDist(-1,z1,x1,fqa,fqa0)*PlusDist(-1,z2,x2,fqb,fqb0) + 
     -     CF**2*d1*G0*Lmu**2*dble(4)*PlusDist(-1,z1,x1,fqa,fqa0)*PlusDist(-1,z2,x2,fqb,fqb0) + 
     -     CF*d1*Lmu*dble(2)*(dble(4)*YLP0qiqi(z1,x1,fqa,fqa0) + gH0qq*PlusDist(-1,z1,x1,fqa,fqa0))*PlusDist(-1,z2,x2,fqb,fqb0) + 
     -     dble(-1)*(-(Lmu*(CF*G0*Lmu + gH0qq2*dble(2))) + cH1qq*dble(4))*
     -      (dble(4)*YLP0qiqi(z1,x1,fqa,fqa0) + gH0qq*PlusDist(-1,z1,x1,fqa,fqa0))*PlusDist(-1,z2,x2,fqb,fqb0) + 
     -     dble(-2)*(dble(4)*YLI1qgYLP0gq(z1,x1,fqa,fqa0) + (gH0qq + beta0*dble(-2))*YLI1qiqi(z1,x1,fqa,fqa0) + 
     -        dble(4)*YLI1qiqiYLP0qiqi(z1,x1,fqa,fqa0) + dble(8)*YLP1qiqi(z1,x1,fqa,fqa0) + gH1qq*PlusDist(-1,z1,x1,fqa,fqa0))*
     -      PlusDist(-1,z2,x2,fqb,fqb0) + dble(-2)*YLI1qiqi(z1,x1,fqa,fqa0)*
     -      (dble(4)*YLP0qiqi(z2,x2,fqb,fqb0) + gH0qq*PlusDist(-1,z2,x2,fqb,fqb0)) + 
     -     CF*d1*Lmu*dble(2)*PlusDist(-1,z1,x1,fqa,fqa0)*(dble(4)*YLP0qiqi(z2,x2,fqb,fqb0) + gH0qq*PlusDist(-1,z2,x2,fqb,fqb0)) + 
     -     dble(-1)*(-(Lmu*(CF*G0*Lmu + gH0qq2*dble(2))) + cH1qq*dble(4))*PlusDist(-1,z1,x1,fqa,fqa0)*
     -      (dble(4)*YLP0qiqi(z2,x2,fqb,fqb0) + gH0qq*PlusDist(-1,z2,x2,fqb,fqb0)) + 
     -     dble(-2)*PlusDist(-1,z1,x1,fqa,fqa0)*(dble(4)*YLI1qgYLP0gq(z2,x2,fqb,fqb0) + 
     -        (gH0qq + beta0*dble(-2))*YLI1qiqi(z2,x2,fqb,fqb0) + dble(4)*YLI1qiqiYLP0qiqi(z2,x2,fqb,fqb0) + 
     -        dble(8)*YLP1qiqi(z2,x2,fqb,fqb0) + gH1qq*PlusDist(-1,z2,x2,fqb,fqb0))) + 
     -  (L + Lmu)**2*dble(dble(1)/dble(4))*(CF*G0*dble(-2)*YLI1qiqi(z2,x2,fqb,fqb0)*PlusDist(-1,z1,x1,fqa,fqa0) + 
     -     dble(8)*YLP0qgYLP0gq(z2,x2,fqb,fqb0)*PlusDist(-1,z1,x1,fqa,fqa0) + 
     -     beta0*dble(-4)*YLP0qiqi(z2,x2,fqb,fqb0)*PlusDist(-1,z1,x1,fqa,fqa0) + 
     -     CF*d1*dble(8)*YLP0qiqi(z2,x2,fqb,fqb0)*PlusDist(-1,z1,x1,fqa,fqa0) + 
     -     gH0qq*dble(8)*YLP0qiqi(z2,x2,fqb,fqb0)*PlusDist(-1,z1,x1,fqa,fqa0) + 
     -     CF*G0*Lmu*dble(8)*YLP0qiqi(z2,x2,fqb,fqb0)*PlusDist(-1,z1,x1,fqa,fqa0) + 
     -     dble(8)*YLP0qiqiYLP0qiqi(z2,x2,fqb,fqb0)*PlusDist(-1,z1,x1,fqa,fqa0) + 
     -     CF*G0*dble(-2)*YLI1qiqi(z1,x1,fqa,fqa0)*PlusDist(-1,z2,x2,fqb,fqb0) + 
     -     dble(8)*YLP0qgYLP0gq(z1,x1,fqa,fqa0)*PlusDist(-1,z2,x2,fqb,fqb0) + 
     -     dble(8)*YLP0qiqiYLP0qiqi(z1,x1,fqa,fqa0)*PlusDist(-1,z2,x2,fqb,fqb0) + 
     -     beta0*CF*d1*dble(-4)*PlusDist(-1,z1,x1,fqa,fqa0)*PlusDist(-1,z2,x2,fqb,fqb0) + 
     -     CF*cH1qq*G0*dble(-4)*PlusDist(-1,z1,x1,fqa,fqa0)*PlusDist(-1,z2,x2,fqb,fqb0) + 
     -     CF*G1*dble(-2)*PlusDist(-1,z1,x1,fqa,fqa0)*PlusDist(-1,z2,x2,fqb,fqb0) + 
     -     beta0*gH0qq*dble(-2)*PlusDist(-1,z1,x1,fqa,fqa0)*PlusDist(-1,z2,x2,fqb,fqb0) + 
     -     beta0*CF*G0*Lmu*dble(-2)*PlusDist(-1,z1,x1,fqa,fqa0)*PlusDist(-1,z2,x2,fqb,fqb0) + 
     -     CF**2*d1**2*dble(2)*PlusDist(-1,z1,x1,fqa,fqa0)*PlusDist(-1,z2,x2,fqb,fqb0) + 
     -     gH0qq**2*dble(2)*PlusDist(-1,z1,x1,fqa,fqa0)*PlusDist(-1,z2,x2,fqb,fqb0) + 
     -     CF*G0*gH0qq2*Lmu*dble(2)*PlusDist(-1,z1,x1,fqa,fqa0)*PlusDist(-1,z2,x2,fqb,fqb0) + 
     -     CF**2*G0**2*Lmu**2*dble(3)*PlusDist(-1,z1,x1,fqa,fqa0)*PlusDist(-1,z2,x2,fqb,fqb0) + 
     -     CF*d1*gH0qq*dble(4)*PlusDist(-1,z1,x1,fqa,fqa0)*PlusDist(-1,z2,x2,fqb,fqb0) + 
     -     CF*G0*gH0qq*Lmu*dble(4)*PlusDist(-1,z1,x1,fqa,fqa0)*PlusDist(-1,z2,x2,fqb,fqb0) + 
     -     CF**2*d1*G0*Lmu*dble(6)*PlusDist(-1,z1,x1,fqa,fqa0)*PlusDist(-1,z2,x2,fqb,fqb0) + 
     -     dble(4)*YLP0qiqi(z1,x1,fqa,fqa0)*(dble(4)*YLP0qiqi(z2,x2,fqb,fqb0) + 
     -        (-beta0 + gH0qq*dble(2) + CF*(d1 + G0*Lmu)*dble(2))*PlusDist(-1,z2,x2,fqb,fqb0)))
        return
        end

!=================================================================================================

        double precision function Cgg1gq(x1,x2,z1,z2,L,Lmu) ! include both gq->qq and qg->qq
        implicit none
        double precision x1,x2,z1,z2,L,Lmu
        include 'paramqt.f'
        double precision YLP0gg,YLP0gq,YLI1gg,YLI1gq
        double precision fg,fg0,fq,fq0
        integer id
        Cgg1gq = 0d0
        ! gq
        fg=PDFbuf1(0,1)
        fg0=PDFbuf1(0,0)
        fq=0d0
        fq0=0d0
        id=int(nf)
        do while (id.gt.0)
                fq=fq+PDFbuf2(id,1)+PDFbuf2(-id,1)
                fq0=fq0+PDFbuf2(id,0)+PDFbuf2(-id,0)
                id=id-1
        end do
        Cgg1gq = Cgg1gq + (YLI1gq(z2,x2,fq,fq0) + (L + Lmu)*dble(2)*YLP0gq(z2,x2,fq,fq0))*PlusDist(-1,z1,x1,fg,fg0)
        ! qg
        fg=PDFbuf2(0,1)
        fg0=PDFbuf2(0,0)
        fq=0d0
        fq0=0d0
        id=int(nf)
        do while (id.gt.0)
                fq=fq+PDFbuf1(id,1)+PDFbuf1(-id,1)
                fq0=fq0+PDFbuf1(id,0)+PDFbuf1(-id,0)
                id=id-1
        end do
        Cgg1gq = Cgg1gq + (YLI1gq(z1,x1,fq,fq0) + (L + Lmu)*dble(2)*YLP0gq(z1,x1,fq,fq0))*PlusDist(-1,z2,x2,fg,fg0)
        return
        end

        double precision function Cgg1gg(x1,x2,z1,z2,L,Lmu)
        implicit none
        double precision x1,x2,z1,z2,L,Lmu
        include 'paramqt.f'
        double precision YLP0gg,YLP0gq,YLI1gg,YLI1gq
        double precision fga,fga0,fgb,fgb0
        fga=PDFbuf1(0,1)
        fga0=PDFbuf1(0,0)
        fgb=PDFbuf2(0,1)
        fgb0=PDFbuf2(0,0)
        Cgg1gg = YLI1gg(z2,x2,fgb,fgb0)*PlusDist(-1,z1,x1,fga,fga0) + (L + Lmu)*dble(2)*YLP0gg(z2,x2,fgb,fgb0)*PlusDist(-1,z1,x1,fga,fga0) + 
     -  dble(dble(1)/dble(2))*(dble(2)*YLI1gg(z1,x1,fga,fga0) + (L + Lmu)*dble(4)*YLP0gg(z1,x1,fga,fga0) + 
     -     (-(CA*G0*L**2) + Lmu*gH0gg2*dble(-2) + CA*d1*L*dble(2) + gH0gg*L*dble(2) + gH0gg*Lmu*dble(2) + cH1gg*dble(4))*
     -      PlusDist(-1,z1,x1,fga,fga0))*PlusDist(-1,z2,x2,fgb,fgb0)
        return
        end

        double precision function Cgg2qq(x1,x2,z1,z2,L,Lmu)
        implicit none
        double precision x1,x2,z1,z2,L,Lmu
        include 'paramqt.f'
        double precision YLP0gg,YLP0gq,YLP1gg,YLP1gq
        double precision YLI1gg,YLI1gq,YLIi1gq,YLIi1gg,YLI2gg,YLI2gq
        double precision YLP0ggYLP0gg,YLP0gqYLP0qg,YLP0gqYLP0qiqi,YLP0ggYLP0gq
        double precision YLI1ggYLP0gg,YLI1gqYLP0qg,YLI1gqYLP0qiqi,YLI1ggYLP0gq
        double precision fqa,fqa0,fqb,fqb0
        integer id
        fqa=0d0
        fqa0=0d0
        fqb=0d0
        fqb0=0d0
        id=int(nf)
        do while (id.gt.0)
                fqa=fqa+PDFbuf1(id,1)+PDFbuf1(-id,1)
                fqa0=fqa0+PDFbuf1(id,0)+PDFbuf1(-id,0)
                fqb=fqb+PDFbuf2(id,1)+PDFbuf2(-id,1)
                fqb0=fqb0+PDFbuf2(id,0)+PDFbuf2(-id,0)
                id=id-1
        end do
        Cgg2qq = YLIi1gq(z1,x1,fqa,fqa0)*YLIi1gq(z2,x2,fqb,fqb0) + 
     -  YLI1gq(z1,x1,fqa,fqa0)*(YLI1gq(z2,x2,fqb,fqb0) + (L + Lmu)*dble(2)*YLP0gq(z2,x2,fqb,fqb0)) + 
     -  (L + Lmu)*dble(2)*YLP0gq(z1,x1,fqa,fqa0)*(YLI1gq(z2,x2,fqb,fqb0) + (L + Lmu)*dble(2)*YLP0gq(z2,x2,fqb,fqb0))
        return
        end

        double precision function Cgg2gq(x1,x2,z1,z2,L,Lmu) ! include both gq->qq and qg->qq
        implicit none
        double precision x1,x2,z1,z2,L,Lmu
        include 'paramqt.f'
        double precision YLP0gg,YLP0gq,YLP1gg,YLP1gq
        double precision YLI1gg,YLI1gq,YLIi1gq,YLIi1gg,YLI2gg,YLI2gq
        double precision YLP0ggYLP0gg,YLP0gqYLP0qg,YLP0gqYLP0qiqi,YLP0ggYLP0gq
        double precision YLI1ggYLP0gg,YLI1gqYLP0qg,YLI1gqYLP0qiqi,YLI1ggYLP0gq
        double precision fg,fg0,fq,fq0
        integer id
        Cgg2gq = 0d0
        ! gq
        fg=PDFbuf1(0,1)
        fg0=PDFbuf1(0,0)
        fq=0d0
        fq0=0d0
        id=int(nf)
        do while (id.gt.0)
                fq=fq+PDFbuf2(id,1)+PDFbuf2(-id,1)
                fq0=fq0+PDFbuf2(id,0)+PDFbuf2(-id,0)
                id=id-1
        end do
        Cgg2gq = Cgg2gq + YLI1gg(z1,x1,fg,fg0)*YLI1gq(z2,x2,fq,fq0) + YLIi1gg(z1,x1,fg,fg0)*YLIi1gq(z2,x2,fq,fq0) -
     -  CA*d1*Lmu*YLI1gq(z2,x2,fq,fq0)*PlusDist(-1,z1,x1,fg,fg0) +
     -  (-(Lmu*(CA*G0*Lmu + gH0gg2*dble(2))) + cH1gg*dble(4))*dble(dble(1)/dble(2))*YLI1gq(z2,x2,fq,fq0)*
     -   PlusDist(-1,z1,x1,fg,fg0) + YLI2gq(z2,x2,fq,fq0)*PlusDist(-1,z1,x1,fg,fg0) +
     -  CA*G0*(Zeta3*dble(-4) + (L + Lmu)**3*dble(-1))*YLP0gq(z2,x2,fq,fq0)*PlusDist(-1,z1,x1,fg,fg0) +
     -  (L + Lmu)**2*dble(dble(1)/dble(2))*(dble(8)*YLP0gg(z1,x1,fg,fg0)*YLP0gq(z2,x2,fq,fq0) +
     -     (-(CA*G0*YLI1gq(z2,x2,fq,fq0)) + dble(4)*YLP0ggYLP0gq(z2,x2,fq,fq0) + beta0*dble(-2)*YLP0gq(z2,x2,fq,fq0) +
     -        CA*d1*dble(4)*YLP0gq(z2,x2,fq,fq0) + gH0gg*dble(4)*YLP0gq(z2,x2,fq,fq0) + CA*G0*Lmu*dble(4)*YLP0gq(z2,x2,fq,fq0) +
     -        dble(4)*YLP0gqYLP0qiqi(z2,x2,fq,fq0))*PlusDist(-1,z1,x1,fg,fg0)) +
     -  (L + Lmu)*dble(-1)*(dble(-2)*YLI1gg(z1,x1,fg,fg0)*YLP0gq(z2,x2,fq,fq0) +
     -     (dble(-2)*YLI1ggYLP0gq(z2,x2,fq,fq0) + dble(-2)*YLI1gqYLP0qiqi(z2,x2,fq,fq0) + CA*G0*Lmu**2*YLP0gq(z2,x2,fq,fq0) +
     -        cH1gg*dble(-4)*YLP0gq(z2,x2,fq,fq0) + CA*d1*Lmu*dble(2)*YLP0gq(z2,x2,fq,fq0) + gH0gg2*Lmu*dble(2)*YLP0gq(z2,x2,fq,fq0) +
     -        dble(-4)*YLP1gq(z2,x2,fq,fq0))*PlusDist(-1,z1,x1,fg,fg0) +
     -     dble(-1)*YLI1gq(z2,x2,fq,fq0)*(dble(2)*YLP0gg(z1,x1,fg,fg0) +
     -        (-beta0 + CA*d1 + gH0gg + CA*G0*Lmu)*PlusDist(-1,z1,x1,fg,fg0)))
        ! qg
        fg=PDFbuf2(0,1)
        fg0=PDFbuf2(0,0)
        fq=0d0
        fq0=0d0
        id=int(nf)
        do while (id.gt.0)
                fq=fq+PDFbuf1(id,1)+PDFbuf1(-id,1)
                fq0=fq0+PDFbuf1(id,0)+PDFbuf1(-id,0)
                id=id-1
        end do
        Cgg2gq = Cgg2gq + YLI1gg(z2,x2,fg,fg0)*YLI1gq(z1,x1,fq,fq0) + YLIi1gg(z2,x2,fg,fg0)*YLIi1gq(z1,x1,fq,fq0) -
     -  CA*d1*Lmu*YLI1gq(z1,x1,fq,fq0)*PlusDist(-1,z2,x2,fg,fg0) +
     -  (-(Lmu*(CA*G0*Lmu + gH0gg2*dble(2))) + cH1gg*dble(4))*dble(dble(1)/dble(2))*YLI1gq(z1,x1,fq,fq0)*
     -   PlusDist(-1,z2,x2,fg,fg0) + YLI2gq(z1,x1,fq,fq0)*PlusDist(-1,z2,x2,fg,fg0) +
     -  CA*G0*(Zeta3*dble(-4) + (L + Lmu)**3*dble(-1))*YLP0gq(z1,x1,fq,fq0)*PlusDist(-1,z2,x2,fg,fg0) +
     -  (L + Lmu)**2*dble(dble(1)/dble(2))*(dble(8)*YLP0gg(z2,x2,fg,fg0)*YLP0gq(z1,x1,fq,fq0) +
     -     (-(CA*G0*YLI1gq(z1,x1,fq,fq0)) + dble(4)*YLP0ggYLP0gq(z1,x1,fq,fq0) + beta0*dble(-2)*YLP0gq(z1,x1,fq,fq0) +
     -        CA*d1*dble(4)*YLP0gq(z1,x1,fq,fq0) + gH0gg*dble(4)*YLP0gq(z1,x1,fq,fq0) + CA*G0*Lmu*dble(4)*YLP0gq(z1,x1,fq,fq0) +
     -        dble(4)*YLP0gqYLP0qiqi(z1,x1,fq,fq0))*PlusDist(-1,z2,x2,fg,fg0)) +
     -  (L + Lmu)*dble(-1)*(dble(-2)*YLI1gg(z2,x2,fg,fg0)*YLP0gq(z1,x1,fq,fq0) +
     -     (dble(-2)*YLI1ggYLP0gq(z1,x1,fq,fq0) + dble(-2)*YLI1gqYLP0qiqi(z1,x1,fq,fq0) + CA*G0*Lmu**2*YLP0gq(z1,x1,fq,fq0) +
     -        cH1gg*dble(-4)*YLP0gq(z1,x1,fq,fq0) + CA*d1*Lmu*dble(2)*YLP0gq(z1,x1,fq,fq0) + gH0gg2*Lmu*dble(2)*YLP0gq(z1,x1,fq,fq0) +
     -        dble(-4)*YLP1gq(z1,x1,fq,fq0))*PlusDist(-1,z2,x2,fg,fg0) +
     -     dble(-1)*YLI1gq(z1,x1,fq,fq0)*(dble(2)*YLP0gg(z2,x2,fg,fg0) +
     -        (-beta0 + CA*d1 + gH0gg + CA*G0*Lmu)*PlusDist(-1,z2,x2,fg,fg0)))
        return
        end

        double precision function Cgg2gg(x1,x2,z1,z2,L,Lmu)
        implicit none
        double precision x1,x2,z1,z2,L,Lmu
        include 'paramqt.f'
        double precision YLP0gg,YLP0gq,YLP1gg,YLP1gq
        double precision YLI1gg,YLI1gq,YLIi1gq,YLIi1gg,YLI2gg,YLI2gq
        double precision YLP0ggYLP0gg,YLP0gqYLP0qg,YLP0gqYLP0qiqi,YLP0ggYLP0gq
        double precision YLI1ggYLP0gg,YLI1gqYLP0qg,YLI1gqYLP0qiqi,YLI1ggYLP0gq
        double precision fga,fga0,fgb,fgb0
        fga=PDFbuf1(0,1)
        fga0=PDFbuf1(0,0)
        fgb=PDFbuf2(0,1)
        fgb0=PDFbuf2(0,0)
        Cgg2gg = YLI1gg(z1,x1,fga,fga0)*YLI1gg(z2,x2,fgb,fgb0) + YLIi1gg(z1,x1,fga,fga0)*YLIi1gg(z2,x2,fgb,fgb0) - 
     -  CA*d1*Lmu*YLI1gg(z2,x2,fgb,fgb0)*PlusDist(-1,z1,x1,fga,fga0) + 
     -  (-(Lmu*(CA*G0*Lmu + gH0gg2*dble(2))) + cH1gg*dble(4))*dble(dble(1)/dble(2))*YLI1gg(z2,x2,fgb,fgb0)*
     -   PlusDist(-1,z1,x1,fga,fga0) + YLI2gg(z2,x2,fgb,fgb0)*PlusDist(-1,z1,x1,fga,fga0) - 
     -  CA*d1*Lmu*YLI1gg(z1,x1,fga,fga0)*PlusDist(-1,z2,x2,fgb,fgb0) + 
     -  (-(Lmu*(CA*G0*Lmu + gH0gg2*dble(2))) + cH1gg*dble(4))*dble(dble(1)/dble(2))*YLI1gg(z1,x1,fga,fga0)*
     -   PlusDist(-1,z2,x2,fgb,fgb0) + YLI2gg(z1,x1,fga,fga0)*PlusDist(-1,z2,x2,fgb,fgb0) - 
     -  CA*d2*Lmu*PlusDist(-1,z1,x1,fga,fga0)*PlusDist(-1,z2,x2,fgb,fgb0) + 
     -  cH1gg*(beta0 + gH0gg2)*Lmu*dble(-2)*PlusDist(-1,z1,x1,fga,fga0)*PlusDist(-1,z2,x2,fgb,fgb0) + 
     -  CA**2*d1**2*Lmu**2*dble(dble(1)/dble(2))*PlusDist(-1,z1,x1,fga,fga0)*PlusDist(-1,z2,x2,fgb,fgb0) + 
     -  Lmu**2*(gH0gg2*(beta0 + gH0gg2) + CA*cH1gg*G0*dble(-2))*dble(dble(1)/dble(2))*PlusDist(-1,z1,x1,fga,fga0)*
     -   PlusDist(-1,z2,x2,fgb,fgb0) + CA*d1*Lmu*(cH1gg*dble(-4) + Lmu*(CA*G0*Lmu + gH0gg2*dble(2)))*dble(dble(1)/dble(2))*
     -   PlusDist(-1,z1,x1,fga,fga0)*PlusDist(-1,z2,x2,fgb,fgb0) + 
     -  (-(Lmu*(CA*G1*Lmu + gH1gg2*dble(2))) + cH2gg*dble(4))*dble(dble(1)/dble(2))*PlusDist(-1,z1,x1,fga,fga0)*
     -   PlusDist(-1,z2,x2,fgb,fgb0) + CA*G0*Lmu**3*(beta0 + gH0gg2*dble(3))*dble(dble(1)/dble(6))*PlusDist(-1,z1,x1,fga,fga0)*
     -   PlusDist(-1,z2,x2,fgb,fgb0) + CA**2*G0**2*Lmu**4*dble(dble(1)/dble(8))*PlusDist(-1,z1,x1,fga,fga0)*
     -   PlusDist(-1,z2,x2,fgb,fgb0) + CA**2*G0**2*((L + Lmu)**4 + (L + Lmu)*Zeta3*dble(16))*dble(dble(1)/dble(8))*
     -   PlusDist(-1,z1,x1,fga,fga0)*PlusDist(-1,z2,x2,fgb,fgb0) + 
     -  CA*G0*(L**3 + Lmu**3 + L**2*Lmu*dble(3) + L*Lmu**2*dble(3) + Zeta3*dble(4))*dble(dble(-1)/dble(6))*
     -   (dble(6)*YLP0gg(z2,x2,fgb,fgb0)*PlusDist(-1,z1,x1,fga,fga0) + 
     -     (dble(6)*YLP0gg(z1,x1,fga,fga0) + (beta0*dble(-2) + gH0gg*dble(3) + CA*(d1 + G0*Lmu)*dble(3))*
     -         PlusDist(-1,z1,x1,fga,fga0))*PlusDist(-1,z2,x2,fgb,fgb0)) + 
     -  (L + Lmu)*dble(dble(-1)/dble(4))*(CA*d1*dble(-4)*YLI1gg(z2,x2,fgb,fgb0)*PlusDist(-1,z1,x1,fga,fga0) + 
     -     CA*G0*Lmu*dble(-4)*YLI1gg(z2,x2,fgb,fgb0)*PlusDist(-1,z1,x1,fga,fga0) + 
     -     dble(-2)*YLI1gg(z2,x2,fgb,fgb0)*(dble(4)*YLP0gg(z1,x1,fga,fga0) + gH0gg*PlusDist(-1,z1,x1,fga,fga0)) + 
     -     CA*d1*dble(-4)*YLI1gg(z1,x1,fga,fga0)*PlusDist(-1,z2,x2,fgb,fgb0) + 
     -     CA*G0*Lmu*dble(-4)*YLI1gg(z1,x1,fga,fga0)*PlusDist(-1,z2,x2,fgb,fgb0) + 
     -     CA*d2*dble(-4)*PlusDist(-1,z1,x1,fga,fga0)*PlusDist(-1,z2,x2,fgb,fgb0) + 
     -     beta0*CA*d1*Lmu*dble(-4)*PlusDist(-1,z1,x1,fga,fga0)*PlusDist(-1,z2,x2,fgb,fgb0) + 
     -     CA*G1*Lmu*dble(-4)*PlusDist(-1,z1,x1,fga,fga0)*PlusDist(-1,z2,x2,fgb,fgb0) + 
     -     CA*d1*dble(2)*(cH1gg*dble(-4) + Lmu*(CA*G0*Lmu + gH0gg2*dble(2)))*PlusDist(-1,z1,x1,fga,fga0)*
     -      PlusDist(-1,z2,x2,fgb,fgb0) + CA*G0*Lmu*dble(2)*(cH1gg*dble(-4) + Lmu*(CA*G0*Lmu + gH0gg2*dble(2)))*
     -      PlusDist(-1,z1,x1,fga,fga0)*PlusDist(-1,z2,x2,fgb,fgb0) + 
     -     CA**2*d1**2*Lmu*dble(4)*PlusDist(-1,z1,x1,fga,fga0)*PlusDist(-1,z2,x2,fgb,fgb0) + 
     -     CA**2*d1*G0*Lmu**2*dble(4)*PlusDist(-1,z1,x1,fga,fga0)*PlusDist(-1,z2,x2,fgb,fgb0) + 
     -     CA*d1*Lmu*dble(2)*(dble(4)*YLP0gg(z1,x1,fga,fga0) + gH0gg*PlusDist(-1,z1,x1,fga,fga0))*PlusDist(-1,z2,x2,fgb,fgb0) + 
     -     dble(-1)*(-(Lmu*(CA*G0*Lmu + gH0gg2*dble(2))) + cH1gg*dble(4))*
     -      (dble(4)*YLP0gg(z1,x1,fga,fga0) + gH0gg*PlusDist(-1,z1,x1,fga,fga0))*PlusDist(-1,z2,x2,fgb,fgb0) + 
     -     dble(2)*((-gH0gg + beta0*dble(2))*YLI1gg(z1,x1,fga,fga0) + dble(-4)*YLI1ggYLP0gg(z1,x1,fga,fga0) + 
     -        nf*dble(-8)*YLI1gqYLP0qg(z1,x1,fga,fga0) + dble(-8)*YLP1gg(z1,x1,fga,fga0) - gH1gg*PlusDist(-1,z1,x1,fga,fga0))*
     -      PlusDist(-1,z2,x2,fgb,fgb0) + dble(-2)*YLI1gg(z1,x1,fga,fga0)*
     -      (dble(4)*YLP0gg(z2,x2,fgb,fgb0) + gH0gg*PlusDist(-1,z2,x2,fgb,fgb0)) + 
     -     CA*d1*Lmu*dble(2)*PlusDist(-1,z1,x1,fga,fga0)*(dble(4)*YLP0gg(z2,x2,fgb,fgb0) + gH0gg*PlusDist(-1,z2,x2,fgb,fgb0)) + 
     -     dble(-1)*(-(Lmu*(CA*G0*Lmu + gH0gg2*dble(2))) + cH1gg*dble(4))*PlusDist(-1,z1,x1,fga,fga0)*
     -      (dble(4)*YLP0gg(z2,x2,fgb,fgb0) + gH0gg*PlusDist(-1,z2,x2,fgb,fgb0)) + 
     -     dble(2)*PlusDist(-1,z1,x1,fga,fga0)*((-gH0gg + beta0*dble(2))*YLI1gg(z2,x2,fgb,fgb0) + 
     -        dble(-4)*YLI1ggYLP0gg(z2,x2,fgb,fgb0) + nf*dble(-8)*YLI1gqYLP0qg(z2,x2,fgb,fgb0) + dble(-8)*YLP1gg(z2,x2,fgb,fgb0) - 
     -        gH1gg*PlusDist(-1,z2,x2,fgb,fgb0))) + (L + Lmu)**2*dble(dble(1)/dble(4))*
     -   (CA*G0*dble(-2)*YLI1gg(z2,x2,fgb,fgb0)*PlusDist(-1,z1,x1,fga,fga0) + 
     -     beta0*dble(-4)*YLP0gg(z2,x2,fgb,fgb0)*PlusDist(-1,z1,x1,fga,fga0) + 
     -     CA*d1*dble(8)*YLP0gg(z2,x2,fgb,fgb0)*PlusDist(-1,z1,x1,fga,fga0) + 
     -     gH0gg*dble(8)*YLP0gg(z2,x2,fgb,fgb0)*PlusDist(-1,z1,x1,fga,fga0) + 
     -     CA*G0*Lmu*dble(8)*YLP0gg(z2,x2,fgb,fgb0)*PlusDist(-1,z1,x1,fga,fga0) + 
     -     dble(8)*YLP0ggYLP0gg(z2,x2,fgb,fgb0)*PlusDist(-1,z1,x1,fga,fga0) + 
     -     nf*dble(16)*YLP0gqYLP0qg(z2,x2,fgb,fgb0)*PlusDist(-1,z1,x1,fga,fga0) + 
     -     CA*G0*dble(-2)*YLI1gg(z1,x1,fga,fga0)*PlusDist(-1,z2,x2,fgb,fgb0) + 
     -     dble(8)*YLP0ggYLP0gg(z1,x1,fga,fga0)*PlusDist(-1,z2,x2,fgb,fgb0) + 
     -     nf*dble(16)*YLP0gqYLP0qg(z1,x1,fga,fga0)*PlusDist(-1,z2,x2,fgb,fgb0) + 
     -     beta0*CA*d1*dble(-4)*PlusDist(-1,z1,x1,fga,fga0)*PlusDist(-1,z2,x2,fgb,fgb0) + 
     -     CA*cH1gg*G0*dble(-4)*PlusDist(-1,z1,x1,fga,fga0)*PlusDist(-1,z2,x2,fgb,fgb0) + 
     -     CA*G1*dble(-2)*PlusDist(-1,z1,x1,fga,fga0)*PlusDist(-1,z2,x2,fgb,fgb0) + 
     -     beta0*gH0gg*dble(-2)*PlusDist(-1,z1,x1,fga,fga0)*PlusDist(-1,z2,x2,fgb,fgb0) + 
     -     beta0*CA*G0*Lmu*dble(-2)*PlusDist(-1,z1,x1,fga,fga0)*PlusDist(-1,z2,x2,fgb,fgb0) + 
     -     CA**2*d1**2*dble(2)*PlusDist(-1,z1,x1,fga,fga0)*PlusDist(-1,z2,x2,fgb,fgb0) + 
     -     gH0gg**2*dble(2)*PlusDist(-1,z1,x1,fga,fga0)*PlusDist(-1,z2,x2,fgb,fgb0) + 
     -     CA*G0*gH0gg2*Lmu*dble(2)*PlusDist(-1,z1,x1,fga,fga0)*PlusDist(-1,z2,x2,fgb,fgb0) + 
     -     CA**2*G0**2*Lmu**2*dble(3)*PlusDist(-1,z1,x1,fga,fga0)*PlusDist(-1,z2,x2,fgb,fgb0) + 
     -     CA*d1*gH0gg*dble(4)*PlusDist(-1,z1,x1,fga,fga0)*PlusDist(-1,z2,x2,fgb,fgb0) + 
     -     CA*G0*gH0gg*Lmu*dble(4)*PlusDist(-1,z1,x1,fga,fga0)*PlusDist(-1,z2,x2,fgb,fgb0) + 
     -     CA**2*d1*G0*Lmu*dble(6)*PlusDist(-1,z1,x1,fga,fga0)*PlusDist(-1,z2,x2,fgb,fgb0) + 
     -     dble(4)*YLP0gg(z1,x1,fga,fga0)*(dble(4)*YLP0gg(z2,x2,fgb,fgb0) + 
     -        (-beta0 + gH0gg*dble(2) + CA*(d1 + G0*Lmu)*dble(2))*PlusDist(-1,z2,x2,fgb,fgb0)))
        return
        end

!=================================================================================================

ccc     "odr" specifies the order of the coefficient in alphaS
ccc     "i1,i2" denotes the incoming parton flavors
ccc     "x1,x2" records the Bjorken scales
ccc     "z1,z2" are two random variables from 0 to 1 for convolution
ccc     "qTsq" is the qT cutoff
ccc     "pmod" flags the process type: 0 for DY/Higgs, 1 for diphoton
 
        double precision function Kqq(odr,i1,i2,x1,x2,z1,z2,as,mur,muf,mom,qTsq,pmod)
        implicit none
        integer odr,pmod
        integer i1,i2
        double precision x1,x2,z1,z2,as,mur,muf
        double precision mom(1:4,1:4)
        double precision qTsq,Qsq,L,Lmu,Lrf
        double precision ck0,ck1,ck2
        double precision Cqq1qiqi,Cqq1qg,Cqq2gg,Cqq2qiqj,Cqq2qiqbi,Cqq2qg,Cqq2qiqi
        include 'paramqt.f'

        call UpdateParam(odr,mom,pmod)

        Qsq = (mom(1,4)+mom(2,4))**2-(mom(1,3)+mom(2,3))**2
        ck0 = PDFbuf1(i1,0)*PDFbuf2(i2,0)
        L = log(qTsq/Qsq)
        Lmu = log(Qsq/muf/muf)
        Lrf = 0d0
        if ( mur.ne.muf ) Lrf = log(mur/muf)
        Kqq = 1d0
        if      ( odr.eq.1 ) then ! 1 + O(alphaS)
                ck1 = Cqq1qg(i1,i2,x1,x2,z1,z2,L,Lmu)+Cqq1qiqi(i1,i2,x1,x2,z1,z2,L,Lmu)
                Kqq = ck1/ck0*as/4d0/Pi + 1d0
        else if ( odr.eq.-1 ) then ! O(alphaS) alone
                ck1 = Cqq1qg(i1,i2,x1,x2,z1,z2,L,Lmu)+Cqq1qiqi(i1,i2,x1,x2,z1,z2,L,Lmu)
                Kqq = ck1/ck0*as/4d0/Pi
        else if ( odr.eq.2 ) then ! 1 + O(alphaS) + O(alphaS^2)
                ck1 = Cqq1qg(i1,i2,x1,x2,z1,z2,L,Lmu)+Cqq1qiqi(i1,i2,x1,x2,z1,z2,L,Lmu)
                ck2 = Cqq2gg(i1,i2,x1,x2,z1,z2,L,Lmu)+Cqq2qiqj(i1,i2,x1,x2,z1,z2,L,Lmu)+Cqq2qiqbi(i1,i2,x1,x2,z1,z2,L,Lmu)
     -               +Cqq2qg(i1,i2,x1,x2,z1,z2,L,Lmu)+Cqq2qiqi(i1,i2,x1,x2,z1,z2,L,Lmu)
                Kqq = (ck2+2d0*beta0*Lrf*ck1)/ck0
                Kqq = Kqq*as/4d0/Pi + ck1/ck0
                Kqq = Kqq*as/4d0/Pi + 1d0
        else if ( odr.eq.-2 ) then ! O(alphaS^2) alone
                if ( mur.ne.muf ) ck1 = Cqq1qg(i1,i2,x1,x2,z1,z2,L,Lmu)+Cqq1qiqi(i1,i2,x1,x2,z1,z2,L,Lmu)
                ck2 = Cqq2gg(i1,i2,x1,x2,z1,z2,L,Lmu)+Cqq2qiqj(i1,i2,x1,x2,z1,z2,L,Lmu)+Cqq2qiqbi(i1,i2,x1,x2,z1,z2,L,Lmu)
     -               +Cqq2qg(i1,i2,x1,x2,z1,z2,L,Lmu)+Cqq2qiqi(i1,i2,x1,x2,z1,z2,L,Lmu)
                Kqq = (ck2+2d0*beta0*Lrf*ck1)/ck0
                Kqq = Kqq*(as/4d0/Pi)**2
        end if
        
        Kqq=Kqq*ck0
        
        return
        end

        double precision function Kgg(odr,x1,x2,z1,z2,as,mur,muf,mt,mom,qTsq,pmod)
        implicit none
        integer odr,pmod
        double precision x1,x2,z1,z2,as,mur,muf,mt
        double precision mom(1:4,1:4)
        double precision qTsq,Qsq,L,Lmu,Lrf
        double precision ck0,ck1,ck2
        double precision h1t,h2t
        double precision Cgg1gq,Cgg1gg,Cgg2qq,Cgg2gq,Cgg2gg
        include 'paramqt.f'
        call LoadAllPDF(x1,x2,z1,z2,muf)
        call UpdateParam(odr,mom,pmod)
        Qsq = (mom(1,4)+mom(2,4))**2-(mom(1,3)+mom(2,3))**2
        ck0 = PDFbuf1(0,0)*PDFbuf2(0,0)
        L = log(qTsq/Qsq)
        Lmu = log(Qsq/muf/muf)
        Lrf = 0d0
        if ( mur.ne.muf ) Lrf = log(mur/muf)
        Kgg = 1d0
        if      ( odr.eq.1 ) then
                h1t = 5d0*CA-3d0*CF
                ck1 = Cgg1gq(x1,x2,z1,z2,L,Lmu)+Cgg1gg(x1,x2,z1,z2,L,Lmu)
                Kgg = (ck1/ck0+2d0*beta0*Lrf) + 2d0*h1t
                Kgg = Kgg*as/4d0/Pi + 1d0
        else if ( odr.eq.-1 ) then
                h1t = 5d0*CA-3d0*CF
                ck1 = Cgg1gq(x1,x2,z1,z2,L,Lmu)+Cgg1gg(x1,x2,z1,z2,L,Lmu)
                Kgg = (ck1/ck0+2d0*beta0*Lrf) + 2d0*h1t
                Kgg = Kgg*as/4d0/Pi
        else if ( odr.eq.2 ) then
                h1t = 5d0*CA-3d0*CF
                h2t = (beta0*h1t-beta1)*2d0*log(mur/mt)
     -                -100d0*CA*CF/3d0-5d0*CA*TF/6d0-4d0*CF*TF/3d0-47d0*CA*nf*TF/9d0-5d0*CF*nf*TF+1063d0*CA**2/36d0+27d0*CF**2/2d0
                ck1 = Cgg1gq(x1,x2,z1,z2,L,Lmu)+Cgg1gg(x1,x2,z1,z2,L,Lmu)
                ck2 = Cgg2qq(x1,x2,z1,z2,L,Lmu)+Cgg2gq(x1,x2,z1,z2,L,Lmu)+Cgg2gg(x1,x2,z1,z2,L,Lmu)
                Kgg = ((ck2+4d0*beta0*Lrf*ck1)/ck0+4d0*((beta0*Lrf)**2)+6d0*beta1*Lrf) + 2d0*h2t+h1t**2
                Kgg = Kgg*as/4d0/Pi + (1d0+2d0*h1t*as/4d0/Pi)*(ck1/ck0+2d0*beta0*Lrf) + 2d0*h1t
                Kgg = Kgg*as/4d0/Pi + 1d0
        else if ( odr.eq.-2 ) then
                h1t = 5d0*CA-3d0*CF
                h2t = (beta0*h1t-beta1)*2d0*log(mur/mt)
     -                -100d0*CA*CF/3d0-5d0*CA*TF/6d0-4d0*CF*TF/3d0-47d0*CA*nf*TF/9d0-5d0*CF*nf*TF+1063d0*CA**2/36d0+27d0*CF**2/2d0
                ck1 = Cgg1gq(x1,x2,z1,z2,L,Lmu)+Cgg1gg(x1,x2,z1,z2,L,Lmu)
                ck2 = Cgg2qq(x1,x2,z1,z2,L,Lmu)+Cgg2gq(x1,x2,z1,z2,L,Lmu)+Cgg2gg(x1,x2,z1,z2,L,Lmu)
                Kgg = ((ck2+4d0*beta0*Lrf*ck1)/ck0+4d0*((beta0*Lrf)**2)+6d0*beta1*Lrf) + 2d0*h2t+h1t**2
                Kgg = Kgg + 2d0*h1t*(ck1/ck0+2d0*beta0*Lrf)
                Kgg = Kgg*(as/4d0/Pi)**2
        end if
        
        Kgg=Kgg*ck0
        
        return
        end


