      subroutine Hampgggsq(order,p1,p2,p3,p4,p5,msq0,hard1,hard2,
     & msq0qq,hard1qq,hard2qq)
c--- returns msq0, hard1, hard2 which can be used to reconstruct the
c--- complete hard function for ab -> H+c using:
c---  H = msq0 * (1 + [as/2/pi]*hard1 + [as/2/pi]^2*hard2)
      implicit none
      include 'types.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'hpls.f'
      include 'scale.f'
      include 'constants.f'
      include 'masses.f'
      include 'scet_const.f'
!     This routine calculates the amplitude for a right-handed quark
!     0--> q^+(1)+qb^1(2)+g^+(3)+l^-(4)+lb^+(5)
!     according to Eq.22 of 1309.3245v3
!     index is helicity of outgoing gluon line
      real(dp):: uu,vv,qsq,xlf,hard1,hard2
      integer:: order,p1,p2,p3,p4,p5,q1,q2,q3,iperm,region,k
      complex(dp):: al,be,ga,de,alC,beC,gaC,deC,sq0,sq1
      complex(dp):: amp(0:2,2),ampqqbgll,C0,lnrat
      real(dp):: hnppp,hnppm,hnqqg,
     & s12,s13,s23,msq0,msq1,msq2,res0,res1,res2,
     & evolve1Lgg,evolve2Lgg,evolve1Lqqb,evolve2Lqqb,
     & c1,res0qq,res1qq,res2qq,msq0qq,msq1qq,msq2qq,hard1qq,hard2qq,Ct1,Ct2
      complex(dp):: alp0,bet0,gam0,alp1,bet1,gam1,alp2,bet2,gam2
      real(dp):: alR(0:2),alI(0:2),
     & Halpha_2a1re,Halpha_2a1im,
     & beR(0:2),beI(0:2),
     & Hbeta_2a1re,Hbeta_2a1im,
     & Hbeta_4a1re,Hbeta_4a1im
      real(dp):: gaR(0:2),gaI(0:2),
     & Hgamma_2a1re,Hgamma_2a1im,
     & Hgamma_3a1re,Hgamma_3a1im,
     & Hgamma_4a1re,Hgamma_4a1im
      real(dp):: 
     & Halpha_2a0re,Halpha_2a0im,
     & Hbeta_2a0re,Hbeta_2a0im,
     & Hbeta_4a0re,Hbeta_4a0im
      real(dp):: 
     & Hgamma_2a0re,Hgamma_2a0im,
     & Hgamma_3a0re,Hgamma_3a0im,
     & Hgamma_4a0re,Hgamma_4a0im
      real(dp):: 
     & Halpha_2a2re,Halpha_2a2im,
     & Hbeta_2a2re,Hbeta_2a2im,
     & Hbeta_4a2re,Hbeta_4a2im
      real(dp):: 
     & Hgamma_2a2re,Hgamma_2a2im,
     & Hgamma_3a2re,Hgamma_3a2im,
     & Hgamma_4a2re,Hgamma_4a2im
      real(dp):: LgamHgg0,LgamHgg1,LgamHqq0,LgamHqq1,Lrat
      real(dp), parameter::
     & UGamHgg0= 6*CA,
     & UGamHgg1=
     &     - 40/3._dp*CA*TR*nf
     &     + 134/3._dp*CA**2
     &     - 2*CA**2*pisq,
     & UGamHqq0= + 2*CA + 4*CF,
     & UGamHqq1=
     &     - 40/9._dp*CA*TR*nf
     &     + 134/9._dp*CA**2
     &     - 2/3._dp*CA**2*pisq
     &     - 80/9._dp*CF*TR*nf
     &     + 268/9._dp*CF*CA
     &     - 4/3._dp*CF*CA*pisq
c--- statement functions for hn of Eqs.(50), (51)
      hnppp(s12,s23,s13,qsq,alp0,alp1)
     & =qsq**4/(two*s12*s23*s13)*real(alp0*conjg(alp1),dp)
      hnppm(s12,s23,s13,qsq,bet0,bet1)
     & =s12**3/(two*s23*s13)*real(bet0*conjg(bet1),dp)
c--- statement function for qqg hn of Eq.(54)
      hnqqg(s12,s23,s13,gam0,gam1)
     & =s23**2/(two*s12)*real(gam0*conjg(gam1),dp)
c--- statement function for one-loop evolution
      evolve1Lgg(musq)=
     &    -UGamHgg0*log(abs(s(p1,p2)/musq))**2/two
     &    -LgamHgg0*log(abs(s(p1,p2)/musq))
c--- statement function for two-loop evolution
      evolve2Lgg(musq,c1)=
     & + log(abs(s(p1,p2)/musq))**4/eight * (
     &    + UGamHgg0**2
     &    )
     & + log(abs(s(p1,p2)/musq))**3/six * (
     &    + 3*UGamHgg0*LgamHgg0
     &    + be0*UGamHgg0
     &    )
     & + log(abs(s(p1,p2)/musq))**2/two * (
     &    + LgamHgg0**2
     &    - UGamHgg1
     &    - c1*UGamHgg0
     &    + be0*LgamHgg0
     &    )
     & + log(abs(s(p1,p2)/musq)) * (
     &    - LgamHgg1
     &    - c1*LgamHgg0
     &    - be0*c1
     &    )
c--- statement function for one-loop evolution
      evolve1Lqqb(musq)=
     &    -UGamHqq0*log(abs(s(p1,p2)/musq))**2/two
     &    -LgamHqq0*log(abs(s(p1,p2)/musq))
c--- statement function for two-loop evolution
      evolve2Lqqb(musq,c1)=
     & + log(abs(s(p1,p2)/musq))**4/eight * (
     &    + UGamHqq0**2
     &    )
     & + log(abs(s(p1,p2)/musq))**3/six * (
     &    + 3*UGamHqq0*LgamHqq0
     &    + be0*UGamHqq0
     &    )
     & + log(abs(s(p1,p2)/musq))**2/two * (
     &    + LgamHqq0**2
     &    - UGamHqq1
     &    - c1*UGamHqq0
     &    + be0*LgamHqq0
     &    )
     & + log(abs(s(p1,p2)/musq)) * (
     &    - LgamHqq1
     &    - c1*LgamHqq0
     &    - be0*c1
     &    )
      
      xlf=real(nf,dp)
      
      if (1 == 2) then
c---- ! Check Becher et al. result at a single PS point
        if ((p1 == 2) .and. (p2 == 1)) then
        s(p1,p2)=1._dp
        s(p1,p3)=-0.4_dp
        s(p4,p5)=0.1_dp**2
        s(p2,p3)=s(p4,p5)-s(p1,p2)-s(p1,p3)
        s(p2,p1)=s(p1,p2)
        s(p3,p1)=s(p1,p3)
        s(p3,p2)=s(p2,p3)
        endif
        if ((p1 == 5) .and. (p2 == 1)) then
        s(p2,p3)=1._dp
        s(p1,p3)=-0.4_dp
        s(p4,p5)=0.1_dp**2
        s(p1,p2)=s(p4,p5)-s(p1,p3)-s(p2,p3)
        s(p2,p1)=s(p1,p2)
        s(p3,p1)=s(p1,p3)
        s(p3,p2)=s(p2,p3)
        endif
        if ((p1 == 2) .and. (p2 == 5)) then
        s(p1,p3)=1._dp
        s(p2,p3)=-0.4_dp
        s(p4,p5)=0.1_dp**2
        s(p1,p2)=s(p4,p5)-s(p1,p3)-s(p2,p3)
        s(p2,p1)=s(p1,p2)
        s(p3,p1)=s(p1,p3)
        s(p3,p2)=s(p2,p3)
        endif
        scale=0.6_dp
        musq=scale**2
      endif
      
c      write(6,*) 'region',region
c      write(6,*) 's(p1,p2)',s(p1,p2)
c      write(6,*) 's(p1,p3)',s(p1,p3)
c      write(6,*) 's(p2,p3)',s(p2,p3)

      qsq=s(p1,p2)+s(p1,p3)+s(p2,p3)

c--- logarithm appearing in evolution expressions
c--- note: take absolute value in order to work for all crossings
      Lrat=log(abs(s(p1,p2)**2/(s(p1,p3)*s(p2,p3))))
c--- constants for evolution components
      LgamHgg0 = - 2*CA*Lrat
      LgamHgg1 =
     &     + 256/9._dp*CA*TR*nf
     &     - 2/3._dp*CA*TR*nf*pisq
     &     + 40/9._dp*CA*Lrat*TR*nf
     &     - 692/9._dp*CA**2
     &     + 6*CA**2*zeta3
     &     + 11/6._dp*CA**2*pisq
     &     - 134/9._dp*CA**2*Lrat
     &     + 2/3._dp*CA**2*Lrat*pisq
     &     + 12*CF*TR*nf
     &     + 3*be1
     &     + 2*be1 ! Ct correction term, cf. Eq.(B.5)

c--- note: different from functions in W/Z routines due to additional
c--- beta-function factors arising from more powers of alpha-s at LO
      LgamHqq0 =
     &     - 2*CA*Lrat
     &     - 6*CF
     &     + 2*be0    
      LgamHqq1 =
     &     + 256/27._dp*CA*TR*nf
     &     - 2/9._dp*CA*TR*nf*pisq
     &     + 40/9._dp*CA*Lrat*TR*nf
     &     - 692/27._dp*CA**2
     &     + 2*CA**2*zeta3
     &     + 11/18._dp*CA**2*pisq
     &     - 134/9._dp*CA**2*Lrat
     &     + 2/3._dp*CA**2*Lrat*pisq
     &     + 368/27._dp*CF*TR*nf
     &     + 4/3._dp*CF*TR*nf*pisq
     &     - 961/27._dp*CF*CA
     &     + 52*CF*CA*zeta3
     &     - 11/3._dp*CF*CA*pisq
     &     - 3*CF**2
     &     - 48*CF**2*zeta3
     &     + 4*CF**2*pisq
     &     + 3*be1
     &     + 2*be1 ! Ct correction term, cf. Eq.(B.5)
    
      res0=zip
      res1=zip
      res2=zip
      res0qq=zip
      res1qq=zip
      res2qq=zip
      do iperm=1,4
      if     (iperm == 1) then
        q1=p1
        q2=p2
        q3=p3
      elseif (iperm == 2) then
        q1=p3
        q2=p2
        q3=p1
      elseif (iperm == 3) then
        q1=p3
        q2=p1
        q3=p2
      elseif (iperm == 4) then
        q1=p2
        q2=p1
        q3=p3
      endif

c--- determine correct region from invariants
      if (s(q1,q2) > 0) then
        region=2
      elseif (s(q1,q3) > 0) then
        region=3
      elseif (s(q2,q3) > 0) then
        region=4
      else
        write(6,*) 'region not found:'
        write(6,*) 's12,s13,s23',s(q1,q2),s(q1,q3),s(q2,q3)
        stop
      endif

      if (region == 2) then
        uu=-s(q1,q3)/s(q1,q2)
        vv=+qsq/s(q1,q2)
      elseif (region == 3) then
        uu=-s(q2,q3)/s(q1,q3)
        vv=+qsq/s(q1,q3)
      elseif (region == 4) then
        uu=-s(q1,q3)/s(q2,q3)
        vv=+qsq/s(q2,q3)
      endif
c      write(6,*) 'iperm,region',iperm,region
      
      if (order > 0) then
c--- fill arrays for 2DHPLs
        call tdhpl(uu,vv,2*order,G1,G2,G3,G4,H1,H2,H3,H4)
      else
c--- fill with dummy values
        G1=zip
        G2=zip
        G3=zip
        G4=zip
        H1=zip
        H2=zip
        H3=zip
        H4=zip
      endif

c--- Coefficients for Region 2
      if (region == 2) then
c--- region 2: relevant for gg->Hg (alpha,beta) and qa->Hg (gamma)
      alR(0)=Halpha_2a0re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf)
      alI(0)=Halpha_2a0im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf)
      beR(0)=Hbeta_2a0re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf)
      beI(0)=Hbeta_2a0im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf)
      gaR(0)=Hgamma_2a0re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf)
      gaI(0)=Hgamma_2a0im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf)
      if (order >= 1) then
      alR(1)=Halpha_2a1re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf)
      alI(1)=Halpha_2a1im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf)
      beR(1)=Hbeta_2a1re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf)
      beI(1)=Hbeta_2a1im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf)
      gaR(1)=Hgamma_2a1re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf)
      gaI(1)=Hgamma_2a1im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf)
      endif
      if (order >= 2) then
      alR(2)=Halpha_2a2re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf)
      alI(2)=Halpha_2a2im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf)
      beR(2)=Hbeta_2a2re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf)
      beI(2)=Hbeta_2a2im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf)
      gaR(2)=Hgamma_2a2re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf)
      gaI(2)=Hgamma_2a2im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf)
      endif

c--- Coefficients for Region 3
      elseif (region == 3) then
c--- region 3: relevant for qq->Hg (gamma)
      gaR(0)=Hgamma_3a0re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf)
      gaI(0)=Hgamma_3a0im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf)
      if (order >= 1) then
      gaR(1)=Hgamma_3a1re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf)
      gaI(1)=Hgamma_3a1im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf)
      endif
      if (order >= 2) then
      gaR(2)=Hgamma_3a2re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf)
      gaI(2)=Hgamma_3a2im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf)
      endif

c--- Coefficients for Region 4
      elseif (region == 4) then
c--- region 4: relevant for gg->Hg (beta) and qq->Hg (gamma)
      beR(0)=Hbeta_4a0re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf)
      beI(0)=Hbeta_4a0im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf)
      gaR(0)=Hgamma_4a0re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf)
      gaI(0)=Hgamma_4a0im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf)
      if (order >= 1) then
      beR(1)=Hbeta_4a1re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf)
      beI(1)=Hbeta_4a1im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf)
      gaR(1)=Hgamma_4a1re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf)
      gaI(1)=Hgamma_4a1im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf)
      endif
      if (order >= 2) then
      beR(2)=Hbeta_4a2re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf)
      beI(2)=Hbeta_4a2im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf)
      gaR(2)=Hgamma_4a2re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf)
      gaI(2)=Hgamma_4a2im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf)
      endif
      
      else
        write(6,*) 'Hampgggsq: region should be 2,3 or 4: ',region
        stop

      endif

      alp0=cmplx(alR(0),alI(0),dp)
      bet0=cmplx(beR(0),beI(0),dp)
      gam0=cmplx(gaR(0),gaI(0),dp)
      alp1=cmplx(alR(1),alI(1),dp)
      bet1=cmplx(beR(1),beI(1),dp)
      gam1=cmplx(gaR(1),gaI(1),dp)
      alp2=cmplx(alR(2),alI(2),dp)
      bet2=cmplx(beR(2),beI(2),dp)
      gam2=cmplx(gaR(2),gaI(2),dp)
      
      if (iperm == 1) then
      res0=res0
     & +two*hnppp(s(q1,q2),s(q2,q3),s(q1,q3),qsq,alp0,alp0)
      res1=res1
     & +two*hnppp(s(q1,q2),s(q2,q3),s(q1,q3),qsq,alp0,alp1)
      res2=res2
     & +two*hnppp(s(q1,q2),s(q2,q3),s(q1,q3),qsq,alp0,alp2)
     &     +hnppp(s(q1,q2),s(q2,q3),s(q1,q3),qsq,alp1,alp1)
      endif
      if (iperm <= 3) then
      res0=res0
     & +two*hnppm(s(q1,q2),s(q2,q3),s(q1,q3),qsq,bet0,bet0)
      res1=res1
     & +two*hnppm(s(q1,q2),s(q2,q3),s(q1,q3),qsq,bet0,bet1)
      res2=res2
     & +two*hnppm(s(q1,q2),s(q2,q3),s(q1,q3),qsq,bet0,bet2)
     &     +hnppm(s(q1,q2),s(q2,q3),s(q1,q3),qsq,bet1,bet1)
      endif
      
      if ((iperm == 1) .or. (iperm == 4)) then
c--- note interchange of q1 and q2 in this call wrt. above,
c--- per footnote 3 of 1309.3245
c--- seems to be necessary for qqb~ channel but not for qg
      if ((p1 == 2) .and. (p2 == 1)) then
      res0qq=res0qq
     & +two*hnqqg(s(q1,q2),s(q1,q3),s(q2,q3),gam0,gam0)
      res1qq=res1qq
     & +two*hnqqg(s(q1,q2),s(q1,q3),s(q2,q3),gam0,gam1)
      res2qq=res2qq
     & +two*hnqqg(s(q1,q2),s(q1,q3),s(q2,q3),gam0,gam2)
     &     +hnqqg(s(q1,q2),s(q1,q3),s(q2,q3),gam1,gam1)
      else
      res0qq=res0qq
     & +two*hnqqg(s(q1,q2),s(q2,q3),s(q1,q3),gam0,gam0)
      res1qq=res1qq
     & +two*hnqqg(s(q1,q2),s(q2,q3),s(q1,q3),gam0,gam1)
      res2qq=res2qq
     & +two*hnqqg(s(q1,q2),s(q2,q3),s(q1,q3),gam0,gam2)
     &     +hnqqg(s(q1,q2),s(q2,q3),s(q1,q3),gam1,gam1)
      endif
      endif
      
c      write(6,*) 'alp0,alp1',alp0,alp1
c      write(6,*) 'bet0,bet1',bet0,bet1
c      write(6,*) 'gam0,gam1',gam0,gam1
          
      enddo ! end of iperm
      
      msq0=res0
      msq0qq=res0qq
      if (order > 0) then
        msq1=two*res1
c--- RG running to restore scale dependence from musq=qsq
        msq1=msq1
     &   +(evolve1Lgg(musq)-evolve1Lgg(qsq))*res0
!        write(6,*) 'gg 1L/LO/fourpi',msq1/msq0/fourpi
        msq1qq=two*res1qq
c--- RG running to restore scale dependence from musq=qsq
        msq1qq=msq1qq
     &   +(evolve1Lqqb(musq)-evolve1Lqqb(qsq))*res0qq
!        write(6,*) 'qq 1L/LO/fourpi',msq1qq/msq0qq/fourpi
        if (order == 1) then
c--- if just doing NLO calculation, add correction from 
c--- Wilson coefficient |Ct|^2 here, in units of (as/4/pi)
          msq1=msq1+(5*CA-3*CF)*res0*two
          msq1qq=msq1qq+(5*CA-3*CF)*res0qq*two
        endif
      else
        msq1=zip
        msq1qq=zip
      endif
      if (order > 1) then
        msq2=two*res2
c--- RG running to restore scale dependence from musq=qsq
        c1=two*res1/res0-evolve1Lgg(qsq)
        msq2=msq2
     &   +(evolve2Lgg(musq,c1)
     &    -evolve2Lgg(qsq,c1))*res0
!        write(6,*) 'gg 2L/LO/fourpi**2',msq2/msq0/fourpi**2
        msq2qq=two*res2qq
c--- RG running to restore scale dependence from musq=qsq
        c1=two*res1qq/res0qq-evolve1Lqqb(qsq)
        msq2qq=msq2qq
     &   +(evolve2Lqqb(musq,c1)
     &    -evolve2Lqqb(qsq,c1))*res0qq
!        write(6,*) 'qq 2L/LO/fourpi**2',msq2qq/msq0qq/fourpi**2
!        pause
      else
        msq2=zip
        msq2qq=zip
      endif
      
c--- convert from coefficient of [as/4/pi]^n to coefficient of [as/2/pi]^n
      msq1=msq1/two
      msq1qq=msq1qq/two
      msq2=msq2/four
      msq2qq=msq2qq/four

c--- normalize coefficients to LO
      hard1=msq1/msq0
      hard1qq=msq1qq/msq0qq
      hard2=msq2/msq0
      hard2qq=msq2qq/msq0qq
      
c      hard1=hard1-pisq/twelve*3._dp*CA
c      hard1qq=hard1qq-pisq/twelve*(2._dp*CF+CA)
      
      return
      end


