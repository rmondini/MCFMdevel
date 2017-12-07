! expect f=PDF(x/z)/z and f0=PDF(x)
! assuming x<z<1, i.e an integral of \int_1^x dz
! also work for \int_1^0 dz though (0,x) part is wasted

        double precision function PlusDist(n,z,x,f,f0)
        implicit none
        integer n
        double precision z,x,f,f0
        double precision fz,jc,lnnmz,lnnmx
!        if ( z.lt.x ) then
!                if ( n.eq.-1 ) then
!                        PlusDist = f0
!                else if ( n.gt.-1 ) then
!                        PlusDist = f0*log(1d0-x)**(n+1)/(n+1d0)
!                else
!                        PlusDist = 0d0
!                end if
!        }
!        else
!                if ( n.eq.-1 ) then
!                        PlusDist = f0
!                else if ( n.gt.-1 ) then
!                        PlusDist = f0*log(1d0-x)**(n+1)/(n+1d0) + (f/z-f0)*log(1d0-z)**n/(1d0-z)
!                else
!                        PlusDist = f/z
!                end if
!        end if
!        return
        if ( z.lt.x ) then
                PlusDist = 0d0
                return
        end if
        fz = f/z;
        if ( n.lt.-1 ) then
                PlusDist = fz
                return
        end if
        jc = 1d0/(1d0-x)
        if ( n.eq.-1 ) then
                PlusDist = f0*jc
        else
                lnnmz = 1d0
                if ( n.eq.1 ) then
                        lnnmz = log(1d0-z)
                else if ( n.gt.1 ) then
                        lnnmz = log(1d0-z)**n
                end if
                lnnmx = log(1d0-x)
                if ( n.gt.0 ) lnnmx = lnnmx**(n+1)
                PlusDist = (fz-f0)*lnnmz/(1d0-z)+f0*lnnmx/(n+1d0)*jc;
        end if
        return
        end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        subroutine LoadAllPDF(x1,x2,z1,z2,mu)
        implicit none
        double precision x1,x2,z1,z2,mu
        include 'paramqt.f'
        double precision x1old,x2old,z1old,z2old,muold
!        save x1old,x2old,z1old,z2old,muold
!        data x1old /0d0/
!        data x2old /0d0/
!        data z1old /0d0/
!        data z2old /0d0/
!        data muold /0d0/
        integer ih1,ih2
        double precision fx1(-5:5),fx2(-5:5),fx1z(-5:5),fx2z(-5:5)
        common/density/ih1,ih2

!        if ( x1.eq.x1old .and. x2.eq.x2old .and. z1.eq.z1old .and. z2.eq.z2old .and. mu.eq.muold ) return

        call fdist_select(ih1,1,x1,mu,fx1)
        call fdist_select(ih1,1,x1/z1,mu,fx1z)
        call fdist_select(ih2,2,x2,mu,fx2)
        call fdist_select(ih2,2,x2/z2,mu,fx2z)
        
        PDFbuf1(:,0)=fx1(:)
        PDFbuf1(:,1)=fx1z(:)
        PDFbuf2(:,0)=fx2(:)
        PDFbuf2(:,1)=fx2z(:)
        
        return
        end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        subroutine CalcParam(N_f)
        implicit none
        integer N_f
        include 'paramqt.f'
        nf = dble(N_f)
        sumQsq = 0d0
        if (N_f.ge.6) then
                sumQsq = sumQsq + 3d0*4d0/9d0
        else
                sumQsq = sumQsq + 2d0*4d0/9d0
        end if
        if (N_f.ge.5) then
                sumQsq = sumQsq + 3d0/9d0
        else
                sumQsq = sumQsq + 2d0/9d0
        end if
        beta0 = nf*TF*dble(dble(-4)/dble(3)) + CA*dble(dble(11)/dble(3))
        beta1 = -(nf*TF*(CF*dble(4) + CA*dble(dble(20)/dble(3)))) + CA**2*dble(dble(34)/dble(3))
        beta2 = nf*TF*dble(2)*(CF**2 + CA*CF*dble(dble(-205)/dble(18)) + CA**2*dble(dble(-1415)/dble(54))) + 
     -  nf**2*TF**2*dble(4)*(CF*dble(dble(11)/dble(9)) + CA*dble(dble(79)/dble(54))) + CA**3*dble(dble(2857)/dble(54))
        G0 = dble(4)
        G1 = dble(4)*(nf*TF*dble(dble(-20)/dble(9)) + CA*(Zeta2*dble(-2) + dble(dble(67)/dble(9))))
        G2 = dble(4)*(CF*nf*TF*(Zeta3*dble(16) + dble(dble(-55)/dble(3))) + nf**2*TF**2*dble(dble(-16)/dble(27)) + 
     -    CA*nf*TF*(Zeta3*dble(dble(-56)/dble(3)) + dble(dble(-418)/dble(27)) + Zeta2*dble(dble(80)/dble(9))) + 
     -    CA**2*(Zeta4*dble(22) + Zeta3*dble(dble(22)/dble(3)) + Zeta2*dble(dble(-268)/dble(9)) + dble(dble(245)/dble(6))))
        d1 = dble(0)
        d2 = nf*TF*dble(dble(-224)/dble(27)) + CA*(Zeta3*dble(-28) + dble(dble(808)/dble(27)))
        d3 = CF*nf*TF*(Zeta4*dble(32) + dble(dble(-3422)/dble(27)) + Zeta3*dble(dble(608)/dble(9))) + 
     -  CA*nf*TF*(Zeta4*dble(dble(-40)/dble(3)) + dble(dble(-125252)/dble(729)) + Zeta3*dble(dble(1808)/dble(27)) + 
     -     Zeta2*dble(dble(1648)/dble(81))) + nf**2*TF**2*(Zeta3*dble(dble(128)/dble(9)) + dble(dble(7424)/dble(729))) + 
     -  CA**2*(Zeta5*dble(192) + Zeta4*dble(dble(-154)/dble(3)) + Zeta3*dble(dble(-12328)/dble(27)) + 
     -     Zeta2*dble(dble(-6392)/dble(81)) + Zeta2*Zeta3*dble(dble(176)/dble(3)) + dble(dble(297029)/dble(729)))
        gH0qq = CF*dble(-6)
        gH1qq = CF*(CF*(Zeta3*dble(-48) + dble(-3) + Zeta2*dble(24)) + 
     -    CA*(Zeta3*dble(52) + Zeta2*dble(-22) + dble(dble(-961)/dble(27))) + 
     -    nf*TF*(Zeta2*dble(8) + dble(dble(260)/dble(27))))
        gH2qq = CF**3*(Zeta4*dble(-288) + Zeta3*dble(-136) + dble(-29) + Zeta2*dble(-36) + Zeta5*dble(480) + 
     -     Zeta2*Zeta3*dble(64)) + CA*CF**2*
     -   (Zeta5*dble(-240) + dble(dble(-151)/dble(2)) + Zeta3*dble(dble(-1688)/dble(3)) + Zeta2*Zeta3*dble(-32) + 
     -     Zeta4*dble(dble(988)/dble(3)) + Zeta2*dble(dble(820)/dble(3))) + 
     -  CF**2*nf*TF*(Zeta4*dble(dble(-560)/dble(3)) + Zeta2*dble(dble(-104)/dble(3)) + Zeta3*dble(dble(1024)/dble(9)) + 
     -     dble(dble(5906)/dble(27))) + CA**2*CF*(Zeta5*dble(-272) + Zeta4*dble(-166) + Zeta2*Zeta3*dble(dble(-176)/dble(3)) + 
     -     Zeta2*dble(dble(-14326)/dble(81)) + dble(dble(-139345)/dble(1458)) + Zeta3*dble(dble(7052)/dble(9))) + 
     -  CF*nf**2*TF**2*(Zeta2*dble(dble(-160)/dble(9)) + Zeta3*dble(dble(-64)/dble(27)) + dble(dble(19336)/dble(729))) + 
     -  CA*CF*nf*TF*(Zeta4*dble(88) + Zeta3*dble(dble(-3856)/dble(27)) + dble(dble(-34636)/dble(729)) + 
     -     Zeta2*dble(dble(10376)/dble(81)))
        gH0qq2 = gH0qq
        gH1qq2 = gH1qq
        gH2qq2 = gH2qq
        gH0gg = -dble(2)*beta0 + dble(0)
        gH1gg = -dble(4)*beta1 + CF*nf*TF*dble(-8) + CA*nf*TF*(Zeta2*dble(dble(-8)/dble(3)) + dble(dble(-208)/dble(27))) + 
     -  CA**2*(Zeta3*dble(4) + dble(dble(-160)/dble(27)) + Zeta2*dble(dble(22)/dble(3)))
        gH2gg = -dble(6)*beta2 + CF**2*nf*TF*dble(8) + CF*nf**2*TF**2*dble(dble(176)/dble(9)) + 
     -  CA*CF*nf*TF*(Zeta4*dble(-32) + Zeta2*dble(-8) + Zeta3*dble(dble(-608)/dble(9)) + 
     -     dble(dble(1178)/dble(27))) + CA**2*nf*TF*(Zeta2*dble(dble(-4792)/dble(81)) + Zeta4*dble(dble(328)/dble(3)) + 
     -     dble(dble(-167800)/dble(729)) + Zeta3*dble(dble(1424)/dble(27))) + 
     -  CA*nf**2*TF**2*(Zeta3*dble(dble(-448)/dble(27)) + Zeta2*dble(dble(160)/dble(27)) + dble(dble(24520)/dble(729))) + 
     -  CA**3*(Zeta5*dble(-32) + Zeta4*dble(dble(-638)/dble(3)) + 
     -     Zeta3*(Zeta2*dble(dble(-80)/dble(3)) + dble(dble(244)/dble(3))) + Zeta2*dble(dble(12218)/dble(81)) + 
     -     dble(dble(37045)/dble(729)))
        gH0gg2 = gH0gg + dble(2)*beta0
        gH1gg2 = gH1gg + dble(4)*beta1 
        gH2gg2 = gH2gg + dble(6)*beta2
        cH1qq = CF*(dble(-8) + Zeta2*dble(7))
        cH2qq = CF*(CF*(Zeta3*dble(-30) + Zeta2*dble(-83) + Zeta4*dble(dble(201)/dble(2)) + dble(dble(511)/dble(8))) + 
     -    CA*(Zeta4*dble(-8) + Zeta3*dble(dble(313)/dble(9)) + dble(dble(-51157)/dble(648)) + Zeta2*dble(dble(1061)/dble(18))) + 
     -    nf*TF*(Zeta3*dble(dble(4)/dble(9)) + Zeta2*dble(dble(-182)/dble(9)) + dble(dble(4085)/dble(162))))
        cH3qq = 0d0
        cH1gg = CA*Zeta2*dble(7)
        cH2gg = CF*nf*TF*(Zeta3*dble(16) + dble(dble(-67)/dble(3))) + 
     -  CA*nf*TF*(Zeta2*dble(dble(-50)/dble(3)) + Zeta3*dble(dble(-92)/dble(9)) + dble(dble(-1832)/dble(81))) + 
     -  CA**2*(Zeta3*dble(dble(-143)/dble(9)) + Zeta4*dble(dble(185)/dble(2)) + Zeta2*dble(dble(335)/dble(6)) + 
     -     dble(dble(5105)/dble(162)))
        cH3gg = 0d0
        return
        end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        subroutine UpdateParam(odr,mom,pmod)
        implicit none
        include 'paramqt.f'
        integer odr,pmod
        double precision mom(1:4,1:4)
        double precision s,t,u,X,Y,X2,Y2,X3,Y3,X4,Y4,V,W,EE,AA,DD,BB
        if (pmod.eq.0) return
        if (pmod.ne.0) then
                if (abs(odr).ge.1) then
                        t = (mom(1,4)+mom(3,4))**2-(mom(1,1)+mom(3,1))**2-(mom(1,2)+mom(3,2))**2-(mom(1,3)+mom(3,3))**2
                        u = (mom(1,4)+mom(4,4))**2-(mom(1,1)+mom(4,1))**2-(mom(1,2)+mom(4,2))**2-(mom(1,3)+mom(4,3))**2
                        s = -t-u
                        X = log(-t/s)
                        Y = log(-u/s)
                        cH1qq = (CF*(dble(6)*u*(dble(2)*t + dble(3)*u)*X + dble(6)*t*(dble(3)*t + dble(2)*u)*Y 
     -                          + dble(42)*(dble(-1) + Zeta2)*(t**2 + u**2) + dble(6)*(s**2 + t**2)*X**2 + dble(6)*(s**2 + u**2)*Y**2))
     -                          /(dble(6)*(t**2 + u**2))
                end if
                if (abs(odr).ge.2) then
                       X2 = Li2(-t/s)
                       Y2 = Li2(-u/s)
                       X3 = Li3(-t/s)
                       Y3 = Li3(-u/s)
                       X4 = Li4(-t/s)
                       Y4 = Li4(-u/s)
                       V = Li4(-t/u)
                       W = Li4(-u/t)
                       EE = ((t**2*Zeta2*dble(-1980) + u**2*Zeta2*dble(-1980) + 
     -      X**2*dble(-36)*(u**2 + t*u*dble(26) + t**2*dble(38)) + 
     -      Y**2*dble(-36)*(t**2 + t*u*dble(26) + u**2*dble(38)) + 
     -      X*dble(36)*((s**2 + t**2)*Zeta2*dble(24) - u*(t*dble(38) + u*dble(45))) + 
     -      Y*dble(36)*(t**2*(dble(-45) + Zeta2*dble(24)) + 
     -         t*u*dble(2)*(dble(-19) + Zeta2*dble(24)) + u**2*Zeta2*dble(48)) + 
     -      t**2*Zeta3*dble(72) + u**2*Zeta3*dble(72) + (s**2 + t**2)*X**3*dble(144) + 
     -      (s**2 + u**2)*Y**3*dble(144) + t*u*Zeta2*dble(1728) + t**2*dble(3401) + 
     -      u**2*dble(3401))*dble(dble(1)/dble(162)))/(t**2 + u**2)
                       AA = ((t**3*u**3*Zeta3*dble(-2880) + t**3*u**3*Zeta2*dble(-2220) + 
     -      t**2*u**4*X*Zeta3*dble(-1800) + t**4*u**2*Y*Zeta3*dble(-1800) + 
     -      t**4*u**2*X4*dble(-1440) + t**2*u**4*Y4*dble(-1440) + 
     -      t**4*u**2*X*Y*Zeta2*dble(-1440) + t**3*u**3*X*Y*Zeta2*dble(-1440) + 
     -      t**2*u**4*X*Y*Zeta2*dble(-1440) + t**3*u**3*X*Zeta3*dble(-1440) + 
     -      t**3*u**3*Y*Zeta3*dble(-1440) + t**4*u**2*Zeta3*dble(-1320) + 
     -      t**2*u**4*Zeta3*dble(-1320) + t**4*u**2*Zeta2*dble(-1110) + 
     -      t**2*u**4*Zeta2*dble(-1110) + t**5*u*X4*dble(-720) + 
     -      t**3*u**3*X**2*Y*dble(-720) + t**3*u**3*X*Y**2*dble(-720) + 
     -      t*u**5*Y4*dble(-720) + t**5*u*X*Y*Zeta2*dble(-720) + 
     -      t*u**5*X*Y*Zeta2*dble(-720) + t*u**5*X*Zeta3*dble(-720) + 
     -      t**5*u*Y*Zeta3*dble(-720) + t**2*u**4*X**2*Y*dble(-420) + 
     -      t**4*u**2*X*Y**2*dble(-420) + t**4*u**2*X**2*Y*dble(-360) + 
     -      t**2*u**4*X*Y**2*dble(-360) + t**4*u**2*X*Zeta3*dble(-360) + 
     -      t**2*u**4*Y*Zeta3*dble(-360) + t**4*u**2*X**3*Y*dble(-240) + 
     -      t**2*u**4*X*Y**3*dble(-240) + t**5*u*X3*dble(-120) + 
     -      t**5*u*X**3*Y*dble(-120) + t**3*u**3*X**3*Y*dble(-120) + 
     -      t**3*u**3*X*Y**3*dble(-120) + t*u**5*X*Y**3*dble(-120) + 
     -      t*u**5*Y3*dble(-120) + t**3*u**3*X*Y*dble(-60) + t**4*u**2*X*Y*dble(-30) + 
     -      t**2*u**4*X*Y*dble(-30) + 
     -      t*u**2*Y2*dble(-60)*(t*u**2*X + s**2*u*Zeta2*dble(12) + 
     -         Y*(u**3*dble(-2) + t*u**2*dble(10) + t**3*dble(13) + t**2*u*dble(24))) + 
     -      t**2*u*X2*dble(-60)*(t*(t*u*Y + s**2*Zeta2*dble(12)) + 
     -         X*(t**3*dble(-2) + t**2*u*dble(10) + u**3*dble(13) + t*u**2*dble(24))) + 
     -      t**5*u*X**4*dble(30) + t**3*u**3*X**4*dble(30) + t**3*u**3*Y**4*dble(30) + 
     -      t*u**5*Y**4*dble(30) + t**4*u**2*X**2*dble(45) + t**2*u**4*Y**2*dble(45) + 
     -      t**4*u**2*X**4*dble(60) + t**2*u**4*Y**4*dble(60) + 
     -      t**3*u**3*X**2*dble(90) + t**2*u**4*X**2*dble(90) + 
     -      t**4*u**2*Y**2*dble(90) + t**3*u**3*Y**2*dble(90) + t**5*u*Zeta3*dble(120) + 
     -      t*u**5*Zeta3*dble(120) + t**6*X**2*dble(135) + u**6*Y**2*dble(135) + 
     -      t**5*u*X**2*dble(180) + t*u**5*Y**2*dble(180) + t**5*u*X**2*Y**2*dble(180) + 
     -      t*u**5*X**2*Y**2*dble(180) + t**5*u*X*dble(270) + t**2*u**4*X*dble(270) + 
     -      t**4*u**2*Y*dble(270) + t*u**5*Y*dble(270) + t**2*u**4*X3*Y*dble(360) + 
     -      t**4*u**2*X*Y3*dble(360) + t**5*u*X**2*Zeta2*dble(360) + 
     -      t**3*u**3*X**2*Zeta2*dble(360) + t**3*u**3*Y**2*Zeta2*dble(360) + 
     -      t*u**5*Y**2*Zeta2*dble(360) + t**4*u**2*X**2*Y**2*dble(450) + 
     -      t**2*u**4*X**2*Y**2*dble(450) + t**3*u**3*X**2*Y**2*dble(540) + 
     -      t**4*u**2*X*dble(570) + t**3*u**3*X*dble(570) + t**3*u**3*Y*dble(570) + 
     -      t**2*u**4*Y*dble(570) + t**4*u**2*X3*dble(600) + t**2*u**4*Y3*dble(600) + 
     -      t**5*u*X*Zeta2*dble(600) + t*u**5*Y*Zeta2*dble(600) + 
     -      t**3*u**3*V*dble(720) + t*u**5*V*dble(720) + t**5*u*W*dble(720) + 
     -      t**3*u**3*W*dble(720) + t**2*u**4*X3*dble(720) + t*u**5*X4*dble(720) + 
     -      t**5*u*X3*Y*dble(720) + t**4*u**2*Y3*dble(720) + t*u**5*X*Y3*dble(720) + 
     -      t**5*u*Y4*dble(720) + t**4*u**2*X**2*Zeta2*dble(720) + 
     -      t**2*u**4*Y**2*Zeta2*dble(720) + t**2*u**4*X*Zeta2*dble(780) + 
     -      t**4*u**2*Y*Zeta2*dble(780) + t**4*u**2*X*Zeta2*dble(1320) + 
     -      t**2*u**4*Y*Zeta2*dble(1320) + t**2*u**4*V*dble(1440) + 
     -      t**4*u**2*W*dble(1440) + t**3*u**3*X3*dble(1440) + t**2*u**4*X4*dble(1440) + 
     -      t**3*u**3*X3*Y*dble(1440) + t**3*u**3*Y3*dble(1440) + 
     -      t**3*u**3*X*Y3*dble(1440) + t**4*u**2*Y4*dble(1440) + 
     -      t**3*u**3*X*Zeta2*dble(1440) + t**3*u**3*Y*Zeta2*dble(1440) + 
     -      t**4*u**2*X3*Y*dble(1800) + t**2*u**4*X*Y3*dble(1800) + 
     -      t**5*u*Zeta4*dble(2520) + t*u**5*Zeta4*dble(2520) + 
     -      t**4*u**2*Zeta4*dble(6030) + t**2*u**4*Zeta4*dble(6030) + 
     -      t**3*u**3*Zeta4*dble(7020))*dble(dble(4)/dble(45)))/(s**2*t*u*(t**2 + u**2))
                       DD = ((t**3*u**3*Zeta2*dble(-528840) + t**4*u**2*Zeta2*dble(-455400) + 
     -      t**2*u**4*Zeta2*dble(-455400) + t**2*u**4*X*dble(-300960) + 
     -      t**4*u**2*Y*dble(-300960) + t**3*u**3*X*dble(-286020) + 
     -      t**3*u**3*Y*dble(-286020) + t**4*u**2*X3*Y*dble(-259200) + 
     -      t**3*u**3*X3*Y*dble(-259200) + t**3*u**3*X*Y3*dble(-259200) + 
     -      t**2*u**4*X*Y3*dble(-259200) + t**3*u**3*X3*dble(-220320) + 
     -      t**3*u**3*Y3*dble(-220320) + t**4*u**2*X3*dble(-168480) + 
     -      t**2*u**4*Y3*dble(-168480) + t**3*u**3*X**2*dble(-165060) + 
     -      t**3*u**3*Y**2*dble(-165060) + t**4*u**2*X**2*dble(-156600) + 
     -      t**2*u**4*Y**2*dble(-156600) + t**2*u**4*V*dble(-155520) + 
     -      t**4*u**2*W*dble(-155520) + t**4*u**2*Zeta3*dble(-134640) + 
     -      t**2*u**4*Zeta3*dble(-134640) + t**5*u*Zeta2*dble(-132660) + 
     -      t*u**5*Zeta2*dble(-132660) + t**3*u**3*Zeta3*dble(-108720) + 
     -      t**2*u**4*X4*dble(-103680) + t**2*u**4*X3*Y*dble(-103680) + 
     -      t**4*u**2*X*Y3*dble(-103680) + t**4*u**2*Y4*dble(-103680) + 
     -      t*u**5*X*dble(-92340) + t**5*u*Y*dble(-92340) + t**5*u*X3*Y*dble(-90720) + 
     -      t*u**5*X*Y3*dble(-90720) + t**3*u**3*X**2*Y**2*dble(-84240) + 
     -      t**5*u*Zeta3*dble(-80280) + t*u**5*Zeta3*dble(-80280) + 
     -      t**3*u**3*V*dble(-77760) + t*u**5*V*dble(-77760) + t**5*u*W*dble(-77760) + 
     -      t**3*u**3*W*dble(-77760) + t**4*u**2*X*dble(-77400) + 
     -      t**2*u**4*Y*dble(-77400) + t**2*u**4*X**2*dble(-75960) + 
     -      t**4*u**2*Y**2*dble(-75960) + t**2*u**4*X3*dble(-64800) + 
     -      t**3*u**3*X4*dble(-64800) + t**4*u**2*X**2*Y**2*dble(-64800) + 
     -      t**2*u**4*X**2*Y**2*dble(-64800) + t**4*u**2*Y3*dble(-64800) + 
     -      t**3*u**3*Y4*dble(-64800) + t**2*u**4*X*X3*dble(-51840) + 
     -      t**4*u**2*Y*Y3*dble(-51840) + t*u**5*X4*dble(-38880) + 
     -      t**5*u*Y4*dble(-38880) + t**5*u*X3*dble(-32400) + t*u**5*Y3*dble(-32400) + 
     -      t**5*u*X**2*dble(-32040) + t*u**5*Y**2*dble(-32040) + 
     -      t*u**5*X*X3*dble(-25920) + t**5*u*Y*Y3*dble(-25920) + 
     -      t**5*u*X**2*Y**2*dble(-22680) + t*u**5*X**2*Y**2*dble(-22680) + 
     -      t**3*u**3*X*X3*dble(-12960) + t**2*u**4*X**3*Y*dble(-12960) + 
     -      t*u**5*X3*Y*dble(-12960) + t**4*u**2*X*Y**3*dble(-12960) + 
     -      t**5*u*X*Y3*dble(-12960) + t**3*u**3*Y*Y3*dble(-12960) + 
     -      t*u**5*X**3*Y*dble(-6480) + t**5*u*X*Y**3*dble(-6480) + 
     -      t*u**5*X**2*dble(-6300) + t**5*u*Y**2*dble(-6300) + t**5*u*X**4*dble(-540) + 
     -      t*u**5*Y**4*dble(-540) + 
     -      s**2*t**2*u*X2*dble(-2160)*
     -       (s*X**2*dble(6) + X*dble(-3)*(t*dble(5) + u*dble(16)) + 
     -         t*(Y*dble(-9) + Y**2*dble(3) + Zeta2*dble(24))) + 
     -      t**4*u**2*X**4*dble(1080) + t**2*u**4*Y**4*dble(1080) + 
     -      t**2*u**4*X**4*dble(2160) + t**4*u**2*Y**4*dble(2160) + 
     -      s**2*t*u**2*Y2*(u*Zeta2*dble(-24) + s*Y**2*dble(-6) + u*X**2*dble(-3) + 
     -         u*X*dble(9) + Y*dble(3)*(u*dble(5) + t*dble(16)))*dble(2160) + 
     -      t**3*u**3*X**4*dble(3780) + t**3*u**3*Y**4*dble(3780) + 
     -      t**3*u**3*X**3*Y*dble(6480) + t**3*u**3*X*Y**3*dble(6480) + 
     -      t*u**5*X**3*dble(7920) + t**5*u*Y**3*dble(7920) + t**6*X**2*dble(9720) + 
     -      t*u**5*X**2*Y*dble(9720) + u**6*Y**2*dble(9720) + t**5*u*X*Y**2*dble(9720) + 
     -      t**5*u*X**3*dble(12600) + t*u**5*Y**3*dble(12600) + 
     -      t**5*u*X*X3*dble(12960) + t**5*u*X**3*Y*dble(12960) + 
     -      t*u**5*X*Y**3*dble(12960) + t*u**5*Y*Y3*dble(12960) + 
     -      t*u**5*X**2*Zeta2*dble(12960) + t**5*u*Y**2*Zeta2*dble(12960) + 
     -      t*u**5*X3*dble(19440) + t**5*u*Y3*dble(19440) + t**5*u*X*Zeta2*dble(23760) + 
     -      t*u**5*Y*Zeta2*dble(23760) + t**4*u**2*X*X3*dble(25920) + 
     -      t**5*u*X4*dble(25920) + t**4*u**2*X4*dble(25920) + 
     -      t**5*u*X**2*Y*dble(25920) + t**4*u**2*X**3*Y*dble(25920) + 
     -      t*u**5*X*Y**2*dble(25920) + t**2*u**4*X*Y**3*dble(25920) + 
     -      t**2*u**4*Y*Y3*dble(25920) + t**2*u**4*Y4*dble(25920) + 
     -      t*u**5*Y4*dble(25920) + t*u**5*X*Zeta2*dble(28080) + 
     -      t**5*u*Y*Zeta2*dble(28080) + t**2*u**4*X**3*dble(36000) + 
     -      t**4*u**2*Y**3*dble(36000) + t**5*u*X**2*Zeta2*dble(38880) + 
     -      t*u**5*Y**2*Zeta2*dble(38880) + t**4*u**2*X**3*dble(45360) + 
     -      t**2*u**4*Y**3*dble(45360) + t**5*u*Zeta4*dble(55080) + 
     -      t*u**5*Zeta4*dble(55080) + t**2*u**4*X**2*Y*dble(58320) + 
     -      t**4*u**2*X*Y**2*dble(58320) + t**3*u**3*X**3*dble(60840) + 
     -      t**3*u**3*Y**3*dble(60840) + t**4*u**2*X*Y*dble(64800) + 
     -      t**2*u**4*X*Y*dble(64800) + t**5*u*X*Y*Zeta2*dble(77760) + 
     -      t*u**5*X*Y*Zeta2*dble(77760) + t**4*u**2*X*Zeta3*dble(77760) + 
     -      t**2*u**4*Y*Zeta3*dble(77760) + t**3*u**3*X*Y*dble(90720) + 
     -      t**4*u**2*X**2*Y*dble(90720) + t**2*u**4*X*Y**2*dble(90720) + 
     -      t**2*u**4*X**2*Zeta2*dble(90720) + t**4*u**2*Y**2*Zeta2*dble(90720) + 
     -      t**3*u**3*X**2*Y*dble(113400) + t**3*u**3*X*Y**2*dble(113400) + 
     -      t*u**5*X*Zeta3*dble(116640) + t**5*u*Y*Zeta3*dble(116640) + 
     -      t**4*u**2*X**2*Zeta2*dble(142560) + t**2*u**4*Y**2*Zeta2*dble(142560) + 
     -      t**3*u**3*X**2*Zeta2*dble(181440) + t**3*u**3*Y**2*Zeta2*dble(181440) + 
     -      t**4*u**2*X*Zeta2*dble(207360) + t**2*u**4*Y*Zeta2*dble(207360) + 
     -      t**4*u**2*Zeta4*dble(207360) + t**2*u**4*Zeta4*dble(207360) + 
     -      t**2*u**4*X*Zeta2*dble(216000) + t**4*u**2*Y*Zeta2*dble(216000) + 
     -      t**5*u*dble(217085) + t*u**5*dble(217085) + 
     -      t**4*u**2*X*Y*Zeta2*dble(259200) + t**2*u**4*X*Y*Zeta2*dble(259200) + 
     -      t**3*u**3*X*Zeta3*dble(272160) + t**3*u**3*Y*Zeta3*dble(272160) + 
     -      t**3*u**3*Zeta4*dble(304560) + t**2*u**4*X*Zeta3*dble(311040) + 
     -      t**4*u**2*Y*Zeta3*dble(311040) + t**3*u**3*X*Y*Zeta2*dble(362880) + 
     -      t**3*u**3*X*Zeta2*dble(371520) + t**3*u**3*Y*Zeta2*dble(371520) + 
     -      t**4*u**2*dble(434170) + t**3*u**3*dble(434170) + t**2*u**4*dble(434170))*
     -    dble(dble(-1)/dble(3240)))/(s**2*t*u*(t**2 + u**2))
                       BB = ((t**4*u**4*Zeta2*dble(-125280) + t**5*u**3*Zeta2*dble(-102240) + 
     -      t**3*u**5*Zeta2*dble(-102240) + t**4*u**4*X3*Y*dble(-72000) + 
     -      t**4*u**4*X*Y3*dble(-72000) + t**5*u**3*X3*Y*dble(-69120) + 
     -      t**3*u**5*X*Y3*dble(-69120) + t**4*u**4*X3*dble(-64080) + 
     -      t**4*u**4*Y3*dble(-64080) + t**5*u**3*X3*dble(-50400) + 
     -      t**3*u**5*Y3*dble(-50400) + t**3*u**5*X*dble(-48960) + 
     -      t**5*u**3*Y*dble(-48960) + t**4*u**4*X*dble(-46980) + 
     -      t**4*u**4*Y*dble(-46980) + t**3*u**5*X4*dble(-46080) + 
     -      t**5*u**3*Y4*dble(-46080) + t**4*u**4*X4*dble(-44640) + 
     -      t**4*u**4*Y4*dble(-44640) + t**3*u**5*V*dble(-40320) + 
     -      t**5*u**3*W*dble(-40320) + t**3*u**5*X3*Y*dble(-28800) + 
     -      t**5*u**3*X*Y3*dble(-28800) + t**6*u**2*Zeta2*dble(-26640) + 
     -      t**2*u**6*Zeta2*dble(-26640) + t**3*u**5*X3*dble(-23040) + 
     -      t**6*u**2*X3*Y*dble(-23040) + t**5*u**3*Y3*dble(-23040) + 
     -      t**2*u**6*X*Y3*dble(-23040) + t**4*u**4*V*dble(-20160) + 
     -      t**2*u**6*V*dble(-20160) + t**6*u**2*W*dble(-20160) + 
     -      t**4*u**4*W*dble(-20160) + t**4*u**4*X**2*dble(-20160) + 
     -      t**4*u**4*Y**2*dble(-20160) + t**4*u**4*X**2*Y**2*dble(-18720) + 
     -      t**5*u**3*X**2*dble(-15840) + t**2*u**6*X4*dble(-15840) + 
     -      t**3*u**5*Y**2*dble(-15840) + t**6*u**2*Y4*dble(-15840) + 
     -      t**2*u**6*X*dble(-14580) + t**6*u**2*Y*dble(-14580) + 
     -      t**5*u**3*X4*dble(-14400) + t**5*u**3*X**2*Y**2*dble(-14400) + 
     -      t**3*u**5*X**2*Y**2*dble(-14400) + t**3*u**5*Y4*dble(-14400) + 
     -      t**5*u**3*X*dble(-14040) + t**3*u**5*Y*dble(-14040) + 
     -      t**3*u**5*X**2*dble(-12960) + t**5*u**3*Y**2*dble(-12960) + 
     -      t**6*u**2*X3*dble(-11520) + t**2*u**6*Y3*dble(-11520) + 
     -      t**6*u**2*X*Zeta3*dble(-5760) + t**2*u**6*Y*Zeta3*dble(-5760) + 
     -      t**6*u**2*X**2*Y**2*dble(-5040) + t**2*u**6*X**2*Y**2*dble(-5040) + 
     -      t**6*u**2*Zeta3*dble(-3600) + t**2*u**6*Zeta3*dble(-3600) + 
     -      t**2*u**6*X*X3*dble(-2880) + t**2*u**6*X3*Y*dble(-2880) + 
     -      t**6*u**2*X*Y3*dble(-2880) + t**6*u**2*Y*Y3*dble(-2880) + 
     -      t**2*u**6*X**2*dble(-1800) + t**6*u**2*Y**2*dble(-1800) + 
     -      t**6*u**2*X*dble(-1440) + t**2*u**6*Y*dble(-1440) + 
     -      t**3*u**5*X**3*Y*dble(-1440) + t**5*u**3*X*Y**3*dble(-1440) + 
     -      t**5*u**3*Zeta3*dble(-1440) + t**3*u**5*Zeta3*dble(-1440) + 
     -      t**6*u**2*X**2*dble(-1080) + t**2*u**6*Y**2*dble(-1080) + 
     -      t**2*u**6*X**3*Y*dble(-720) + t**6*u**2*X*Y**3*dble(-720) + 
     -      s**2*t**3*u**2*X2*dble(-720)*
     -       (s*X**2*dble(2) + t*(Y**2 + Y*dble(-3) + Zeta2*dble(12)) + 
     -         X*dble(-2)*(t*dble(8) + u*dble(19))) + t**8*X**4*dble(180) + 
     -      u**8*Y**4*dble(180) + t**2*u**6*X**4*dble(360) + t**6*u**2*Y**4*dble(360) + 
     -      t**7*u*X**3*dble(720) + t*u**7*Y**3*dble(720) + 
     -      s**2*t**2*u**3*Y2*(-(u*X**2) + u*X*dble(3) + 
     -         dble(2)*(-(s*Y**2) + u*Zeta2*dble(-6) + Y*(u*dble(8) + t*dble(19))))*
     -       dble(720) + t**7*u*X**4*dble(1080) + t**2*u**6*X**2*Y*dble(1080) + 
     -      t**6*u**2*X*Y**2*dble(1080) + t*u**7*Y**4*dble(1080) + 
     -      t**7*u*X**2*dble(1440) + t*u**7*Y**2*dble(1440) + 
     -      t**2*u**6*X**3*dble(2160) + t**2*u**6*X3*dble(2160) + 
     -      t**6*u**2*Y**3*dble(2160) + t**6*u**2*Y3*dble(2160) + 
     -      t**3*u**5*X**4*dble(2280) + t**5*u**3*Y**4*dble(2280) + 
     -      t**6*u**2*X**4*dble(3000) + t**2*u**6*Y**4*dble(3000) + 
     -      t**6*u**2*X**3*dble(3960) + t**2*u**6*Y**3*dble(3960) + 
     -      t**8*X**2*Zeta2*dble(4320) + u**8*Y**2*Zeta2*dble(4320) + 
     -      t**4*u**4*Zeta3*dble(4320) + t**4*u**4*X**4*dble(4860) + 
     -      t**4*u**4*Y**4*dble(4860) + t**5*u**3*X**4*dble(5040) + 
     -      t**3*u**5*Y**4*dble(5040) + t**2*u**6*X**2*Zeta2*dble(5760) + 
     -      t**6*u**2*Y**2*Zeta2*dble(5760) + t**5*u**3*X*Zeta3*dble(5760) + 
     -      t**3*u**5*Y*Zeta3*dble(5760) + t**6*u**2*X*X3*dble(8640) + 
     -      t**2*u**6*Y*Y3*dble(8640) + t**7*u*X*Zeta2*dble(8640) + 
     -      t**2*u**6*X*Zeta2*dble(8640) + t**6*u**2*Y*Zeta2*dble(8640) + 
     -      t*u**7*Y*Zeta2*dble(8640) + t**6*u**2*X**2*Y*dble(9000) + 
     -      t**2*u**6*X*Y**2*dble(9000) + t**4*u**4*X**3*Y*dble(9360) + 
     -      t**4*u**4*X*Y**3*dble(9360) + t**3*u**5*X**3*dble(9840) + 
     -      t**5*u**3*Y**3*dble(9840) + t**6*u**2*X**3*Y*dble(10080) + 
     -      t**2*u**6*X*Y**3*dble(10080) + t**5*u**3*X**3*dble(11280) + 
     -      t**3*u**5*Y**3*dble(11280) + t**3*u**5*X**2*Y*dble(13680) + 
     -      t**5*u**3*X*Y**2*dble(13680) + t**5*u**3*X*Y*dble(14400) + 
     -      t**3*u**5*X*Y*dble(14400) + t**4*u**4*X**3*dble(15720) + 
     -      t**4*u**4*Y**3*dble(15720) + t**4*u**4*X*X3*dble(17280) + 
     -      t**4*u**4*Y*Y3*dble(17280) + t**6*u**2*dble(17595) + t**2*u**6*dble(17595) + 
     -      t**4*u**4*X*Y*dble(20160) + t**5*u**3*X**3*Y*dble(20160) + 
     -      t**3*u**5*X*Y**3*dble(20160) + t**5*u**3*X*X3*dble(23040) + 
     -      t**3*u**5*Y*Y3*dble(23040) + t**7*u*X**2*Zeta2*dble(25920) + 
     -      t*u**7*Y**2*Zeta2*dble(25920) + t**2*u**6*X*Zeta3*dble(25920) + 
     -      t**6*u**2*Y*Zeta3*dble(25920) + t**5*u**3*X**2*Y*dble(29520) + 
     -      t**3*u**5*X*Y**2*dble(29520) + t**6*u**2*X*Zeta2*dble(31680) + 
     -      t**2*u**6*Y*Zeta2*dble(31680) + t**4*u**4*dble(32310) + 
     -      t**4*u**4*X**2*Y*dble(33120) + t**4*u**4*X*Y**2*dble(33120) + 
     -      t**5*u**3*dble(33750) + t**3*u**5*dble(33750) + 
     -      t**6*u**2*X*Y*Zeta2*dble(37440) + t**2*u**6*X*Y*Zeta2*dble(37440) + 
     -      t**3*u**5*X**2*Zeta2*dble(41760) + t**5*u**3*Y**2*Zeta2*dble(41760) + 
     -      t**6*u**2*Zeta4*dble(45540) + t**2*u**6*Zeta4*dble(45540) + 
     -      t**4*u**4*X*Zeta3*dble(54720) + t**4*u**4*Y*Zeta3*dble(54720) + 
     -      t**3*u**5*X*Zeta2*dble(59040) + t**5*u**3*Y*Zeta2*dble(59040) + 
     -      t**3*u**5*X*Zeta3*dble(69120) + t**5*u**3*Y*Zeta3*dble(69120) + 
     -      t**6*u**2*X**2*Zeta2*dble(74880) + t**2*u**6*Y**2*Zeta2*dble(74880) + 
     -      t**5*u**3*X*Zeta2*dble(79200) + t**3*u**5*Y*Zeta2*dble(79200) + 
     -      t**4*u**4*X**2*Zeta2*dble(102240) + t**4*u**4*Y**2*Zeta2*dble(102240) + 
     -      t**4*u**4*X*Zeta2*dble(106560) + t**4*u**4*Y*Zeta2*dble(106560) + 
     -      t**5*u**3*Zeta4*dble(114120) + t**3*u**5*Zeta4*dble(114120) + 
     -      t**5*u**3*X**2*Zeta2*dble(119520) + t**3*u**5*Y**2*Zeta2*dble(119520) + 
     -      t**5*u**3*X*Y*Zeta2*dble(132480) + t**3*u**5*X*Y*Zeta2*dble(132480) + 
     -      t**4*u**4*Zeta4*dble(137160) + t**4*u**4*X*Y*Zeta2*dble(190080))*
     -    dble(dble(1)/dble(360)))/(s**2*t**2*u**2*(t**2 + u**2))
            sumQsq=0d0 ! DEBUG -- remove buggy contribution
                       cH2qq = EE*CF*nf*TF + AA*CF*sumQsq*TF + DD*CA*CF + BB*CF**2
                end if
                return
        end if
        return
        end
