      subroutine Zampqqbgsq(order,p1,p2,p3,p4,p5,NFV,za,zb,amps0,hard1,hard2)
c--- returns msq0, hard1, hard2 which can be used to reconstruct the
c--- complete hard function for ab -> Z+c using:
c---  H = amps0 * (1 + [as/2/pi]*hard1 + [as/2/pi]^2*hard2)
c---
c--- note that amps0, hard1 and hard2 are arrays wjere
c---    1st index is helicity of quark line
c---    2nd index is helicity of lepton line
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'hpls.f'
      include 'scale.f'
      include 'scet_const.f'
      include 'ewcharge.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      real(dp):: l(nf),r(nf),le,ln,re,rn,sin2w,qq1,l1,r1,qq2,l2,r2
      common/zcouple/l,r,qq1,l1,r1,qq2,l2,r2,le,ln,re,rn,sin2w
!     This routine calculates the amplitude for a right-handed quark
!     0--> q^+(1)+qb^1(2)+g^+(3)+l^-(4)+lb^+(5)
!     according to Eq.22 of 1309.3245v3
!     1st index is helicity quark line
!     2nd index is helicity lepton line
!     helicity of gluon is summe.e-_dpover
!
!     Order of calculation is passed in via integer "order"
!     Region is passed in via integer "region"
!
      integer:: p1,p2,p3,p4,p5,q1,q2,q3,hq,hl,k,iperm,region,order
      complex(dp):: al,be,ga,de
      complex(dp):: amps(0:2,2,2,2),ampqqbgll
      real(dp):: hn,s12,s13,s23,qsq,res0,res1,res2,msq0,msq1,msq2,
     & evolve1Lqqb,evolve2Lqqb,c1,hard1(2,2),hard2(2,2)
      real(dp):: uu,vv,amps0(2,2),amps1(2,2),amps2(2,2),xlf,NFV
      complex(dp):: alp0,bet0,gam0,alp1,bet1,gam1,alp2,bet2,gam2
      real(dp):: alR(0:2),alI(0:2),
     & Alpha_2a1re,Alpha_2a1im,
     & Alpha_3a1re,Alpha_3a1im,
     & Alpha_4a1re,Alpha_4a1im
      real(dp):: beR(0:2),beI(0:2),
     & Beta_2a1re,Beta_2a1im,
     & Beta_3a1re,Beta_3a1im,
     & Beta_4a1re,Beta_4a1im
      real(dp):: gaR(0:2),gaI(0:2),
     & Gamma_2a1re,Gamma_2a1im,
     & Gamma_3a1re,Gamma_3a1im,
     & Gamma_4a1re,Gamma_4a1im
      real(dp):: 
     & Alpha_2a0re,Alpha_2a0im,
     & Alpha_3a0re,Alpha_3a0im,
     & Alpha_4a0re,Alpha_4a0im
      real(dp):: 
     & Beta_2a0re,Beta_2a0im,
     & Beta_3a0re,Beta_3a0im,
     & Beta_4a0re,Beta_4a0im
      real(dp):: 
     & Gamma_2a0re,Gamma_2a0im,
     & Gamma_3a0re,Gamma_3a0im,
     & Gamma_4a0re,Gamma_4a0im
      real(dp):: 
     & Alpha_2a2re,Alpha_2a2im,
     & Alpha_3a2re,Alpha_3a2im,
     & Alpha_4a2re,Alpha_4a2im
      real(dp):: 
     & Beta_2a2re,Beta_2a2im,
     & Beta_3a2re,Beta_3a2im,
     & Beta_4a2re,Beta_4a2im
      real(dp):: 
     & Gamma_2a2re,Gamma_2a2im,
     & Gamma_3a2re,Gamma_3a2im,
     & Gamma_4a2re,Gamma_4a2im
      real(dp):: LgamHqq0,LgamHqq1,Lrat
      real(dp), parameter::
     & UGamHqq0= + 2*CA + 4*CF,
     & UGamHqq1=
     &     - 40/9._dp*CA*TR*nf
     &     + 134/9._dp*CA**2
     &     - 2/3._dp*CA**2*pisq
     &     - 80/9._dp*CF*TR*nf
     &     + 268/9._dp*CF*CA
     &     - 4/3._dp*CF*CA*pisq
c--- statement function for hn of Eq.(29)
      hn(s12,s23,s13,qsq,alp0,bet0,gam0,alp1,bet1,gam1)
     & =one/(s13*s23*qsq)*(
     & +real(alp0*conjg(alp1),dp)*s12*(2*s12*qsq+s13*s23)
     & +real(bet0*conjg(bet1),dp)*s23*(2*s23*qsq+s12*s13)
     & +real(gam0*conjg(gam1),dp)*s13*(2*s13*qsq+s12*s23)
     & +real((alp0*conjg(bet1)+alp1*conjg(bet0)),dp)*s12*s23*(2*qsq-s13)
     & -real((alp0*conjg(gam1)+alp1*conjg(gam0)),dp)*s12*s13*s23
     & -real((bet0*conjg(gam1)+bet1*conjg(gam0)),dp)*s13*s23*(2*qsq-s12))
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
        scale=0.6_dp
        musq=scale**2
      endif

c      write(6,*) 'region',region
c      write(6,*) 's(p1,p2)',s(p1,p2)
c      write(6,*) 's(p1,p3)',s(p1,p3)
c      write(6,*) 's(p2,p3)',s(p2,p3)

c--- logarithm appearing in evolution expressions
c--- note: take absolute value in order to work for all crossings
      Lrat=log(abs(s(p1,p2)**2/(s(p1,p3)*s(p2,p3))))
c--- constants for evolution components
      LgamHqq0 =
     &     - 2*CA*Lrat
     &     - 6*CF
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
     &     + be1

      res0=zip
      res1=zip
      res2=zip
      do iperm=1,2
      if (iperm == 1) then
        q1=p1
        q2=p2
        q3=p3
      else
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
      endif
c      write(6,*) 'iperm,region',iperm,region

      if (region == 2) then
        uu=-s(q1,q3)/s(q1,q2)
        vv=+s(p4,p5)/s(q1,q2)
      elseif (region == 3) then
        uu=-s(q2,q3)/s(q1,q3)
        vv=+s(p4,p5)/s(q1,q3)
      elseif (region == 4) then
        uu=-s(q1,q3)/s(q2,q3)
        vv=+s(p4,p5)/s(q2,q3)
      endif
            
      if (order > 0) then
c--- fill arrays for 2DHPLs
        call tdhpl(uu,vv,4,G1,G2,G3,G4,H1,H2,H3,H4)
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
        alR(0)=Alpha_2a0re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
        alI(0)=Alpha_2a0im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
        beR(0)=Beta_2a0re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
        beI(0)=Beta_2a0im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
        gaR(0)=Gamma_2a0re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
        gaI(0)=Gamma_2a0im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)

        if (order > 0) then
          alR(1)=Alpha_2a1re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
          alI(1)=Alpha_2a1im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
          beR(1)=Beta_2a1re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
          beI(1)=Beta_2a1im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
          gaR(1)=Gamma_2a1re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
          gaI(1)=Gamma_2a1im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
        endif
      
        if (order > 1) then
          alR(2)=Alpha_2a2re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
          alI(2)=Alpha_2a2im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
          beR(2)=Beta_2a2re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
          beI(2)=Beta_2a2im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
          gaR(2)=Gamma_2a2re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
          gaI(2)=Gamma_2a2im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
        endif
      
c--- Coefficients for Region 3
      elseif (region == 3) then
        alR(0)=Alpha_3a0re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
        alI(0)=Alpha_3a0im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
        beR(0)=Beta_3a0re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
        beI(0)=Beta_3a0im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
        gaR(0)=Gamma_3a0re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
        gaI(0)=Gamma_3a0im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)

        if (order > 0) then
          alR(1)=Alpha_3a1re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
          alI(1)=Alpha_3a1im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
          beR(1)=Beta_3a1re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
          beI(1)=Beta_3a1im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
          gaR(1)=Gamma_3a1re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
          gaI(1)=Gamma_3a1im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
        endif

        if (order > 1) then
          alR(2)=Alpha_3a2re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
          alI(2)=Alpha_3a2im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
          beR(2)=Beta_3a2re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
          beI(2)=Beta_3a2im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
          gaR(2)=Gamma_3a2re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
          gaI(2)=Gamma_3a2im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
        endif

c--- Coefficients for Region 4
      elseif (region == 4) then
        alR(0)=Alpha_4a0re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
        alI(0)=Alpha_4a0im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
        beR(0)=Beta_4a0re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
        beI(0)=Beta_4a0im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
        gaR(0)=Gamma_4a0re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
        gaI(0)=Gamma_4a0im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)

        if (order > 0) then
          alR(1)=Alpha_4a1re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
          alI(1)=Alpha_4a1im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
          beR(1)=Beta_4a1re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
          beI(1)=Beta_4a1im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
          gaR(1)=Gamma_4a1re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
          gaI(1)=Gamma_4a1im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
        endif

        if (order > 1) then
          alR(2)=Alpha_4a2re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
          alI(2)=Alpha_4a2im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
          beR(2)=Beta_4a2re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
          beI(2)=Beta_4a2im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
          gaR(2)=Gamma_4a2re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
          gaI(2)=Gamma_4a2im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
        endif

      else
        write(6,*) 'Zampqqbgsq: region should be 2,3 or 4: ',region
        stop

      endif

      do k=0,2
      al=cmplx(alR(k),alI(k),dp)
      be=cmplx(beR(k),beI(k),dp)
      ga=cmplx(gaR(k),gaI(k),dp)
      de=s(q1,q2)*(al-be-ga)/(2._dp*(s(p1,p2)+s(p2,p3)+s(p3,p1)))

!     1st index is order of amplitude
!     2nd index is helicity quark line
!     3rd index is helicity lepton line
!     4th index is helicity of outgoing gluon line
      if (iperm == 2) then
      amps(k,1,1,2)=-ampqqbgll(p2,p1,p3,p4,p5,al,be,ga,de,za,zb)
      amps(k,2,2,1)=-ampqqbgll(p2,p1,p3,p4,p5,al,be,ga,de,zb,za)
      amps(k,1,2,2)=-ampqqbgll(p2,p1,p3,p5,p4,al,be,ga,de,za,zb)
      amps(k,2,1,1)=-ampqqbgll(p2,p1,p3,p5,p4,al,be,ga,de,zb,za)
      endif
      
      if (iperm == 1) then
      amps(k,1,1,1)=+ampqqbgll(p1,p2,p3,p5,p4,al,be,ga,de,zb,za)
      amps(k,2,2,2)=+ampqqbgll(p1,p2,p3,p5,p4,al,be,ga,de,za,zb)
      amps(k,1,2,1)=+ampqqbgll(p1,p2,p3,p4,p5,al,be,ga,de,zb,za)
      amps(k,2,1,2)=+ampqqbgll(p1,p2,p3,p4,p5,al,be,ga,de,za,zb)
      endif
      
      enddo
 
c---- This bit of code for testing squared MEs against Becher et al.
      alp0=cmplx(alR(0),alI(0),dp)
      bet0=cmplx(beR(0),beI(0),dp)
      gam0=cmplx(gaR(0),gaI(0),dp)
      alp1=cmplx(alR(1),alI(1),dp)
      bet1=cmplx(beR(1),beI(1),dp)
      gam1=cmplx(gaR(1),gaI(1),dp)
      alp2=cmplx(alR(2),alI(2),dp)
      bet2=cmplx(beR(2),beI(2),dp)
      gam2=cmplx(gaR(2),gaI(2),dp)

      res0=res0
     & +two*hn(s(q1,q2),s(q2,q3),s(q1,q3),s(p4,p5),
     &         alp0,bet0,gam0,alp0,bet0,gam0)
     
      res1=res1
     & +two*hn(s(q1,q2),s(q2,q3),s(q1,q3),s(p4,p5),
     &        alp0,bet0,gam0,alp1,bet1,gam1)
     
      res2=res2
     & +two*hn(s(q1,q2),s(q2,q3),s(q1,q3),s(p4,p5),
     &        alp0,bet0,gam0,alp2,bet2,gam2)
     &     +hn(s(q1,q2),s(q2,q3),s(q1,q3),s(p4,p5),
     &        alp1,bet1,gam1,alp1,bet1,gam1)
     
      enddo ! end of iperm loop
      
c---- This bit of code for testing squared MEs against Becher et al.
      msq0=res0
c        write(6,*) 'LO not decayed',
c     &   msq0*esq*(L(2)**2+R(2)**2)/two*gsq*V*aveqq
      if (order > 0) then
        msq1=two*res1
c--- RG running to restore scale dependence from musq=s(p4,p5)
        msq1=msq1
     &   +(evolve1Lqqb(musq)-evolve1Lqqb(s(p4,p5)))*res0
c        write(6,*) '1L/LO/fourpi',msq1/res0/fourpi
c        write(6,*) '1L not decayed',
c     &   msq1*esq*(L(2)**2+R(2)**2)/two*gsq*V*aveqq*ason4pi
      endif
      if (order > 1) then
        msq2=two*res2
c--- RG running to restore scale dependence from musq=s(p4,p5)
        c1=two*res1/res0-evolve1Lqqb(s(p4,p5))
        msq2=msq2
     &   +(evolve2Lqqb(musq,c1)
     &    -evolve2Lqqb(s(p4,p5),c1))*res0
c        write(6,*) '2L/LO/fourpi**2',msq2/res0/fourpi**2
c        write(6,*) '2L not decayed',
c     &   msq2*esq*(L(2)**2+R(2)**2)/two*gsq*V*aveqq*ason4pi**2
      endif
      
      do hq=1,2
      do hl=1,2
        amps0(hq,hl)=abs(amps(0,hq,hl,1))**2+abs(amps(0,hq,hl,2))**2
        if (order > 0) then
          amps1(hq,hl)=two*(
     &     +real(amps(0,hq,hl,1)*conjg(amps(1,hq,hl,1)),dp)
     &     +real(amps(0,hq,hl,2)*conjg(amps(1,hq,hl,2)),dp))
c--- RG running to restore scale dependence from musq=s(p4,p5)
          amps1(hq,hl)=amps1(hq,hl)
     &     +(evolve1Lqqb(musq)-evolve1Lqqb(s(p4,p5)))*amps0(hq,hl)
        else
          amps1(hq,hl)=zip
        endif
        if (order > 1) then
          amps2(hq,hl)=two*(
     &       +real(amps(0,hq,hl,1)*conjg(amps(2,hq,hl,1)),dp)
     &       +real(amps(0,hq,hl,2)*conjg(amps(2,hq,hl,2)),dp))
     &     +real(amps(1,hq,hl,1)*conjg(amps(1,hq,hl,1)),dp)
     &     +real(amps(1,hq,hl,2)*conjg(amps(1,hq,hl,2)),dp)
c--- RG running to restore scale dependence from musq=s(p4,p5)
          c1=amps1(hq,hl)/amps0(hq,hl)-evolve1Lqqb(musq)
          amps2(hq,hl)=amps2(hq,hl)
     &     +(evolve2Lqqb(musq,c1)
     &      -evolve2Lqqb(s(p4,p5),c1))*amps0(hq,hl)
        else
          amps2(hq,hl)=zip
        endif
      enddo
      enddo
      
c--- convert from coefficient of [as/4/pi]^n to coefficient of [as/2/pi]^n
      amps1(:,:)=amps1(:,:)/two
      amps2(:,:)=amps2(:,:)/four
      msq1=msq1/two
      msq2=msq2/four

c--- normalize coefficients to LO
      hard1(:,:)=amps1(:,:)/amps0(:,:)
      hard2(:,:)=amps2(:,:)/amps0(:,:)
      
      return
      end


