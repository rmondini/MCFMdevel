      subroutine Zampgggsq(order,p1,p2,p3,p4,p5,za,zb,amps2)
c--- returns msq0, hard1, hard2 which can be used to reconstruct the
c--- complete hard function for gg -> Z+g using:
c---  H = [as/2/pi]^2 * amps2
c---
c--- note that amps2 is an array where index is helicity of lepton line
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'hpls.f'
      include 'scale.f'
!     This routine calculates the amplitude for
!     0--> g(1)+g(2)+g(3)+l^-(4)+lb^+(5)
!     according to Eq.22 of 1309.3245v3
!     1st index is helicity quark line
!     2nd index is helicity lepton line
!     helicity of gluon is summed over
!
!     Order of calculation is passed in via integer "order"
!
      integer:: p1,p2,p3,p4,p5,q1,q2,q3,hl,ha,hb,hc,j,k,iperm,region,order
      complex(dp):: al,be,ga,de,alC,beC,gaC,deC
      complex(dp):: amps(2,2,2,2),ampgggppp,ampgggpmm
      real(dp):: hnpmm,hnppp,s12,s13,s23,qsq,res0,res1,res2,
     & msq0,msq1,msq2,c1,hard1(2,2),hard2(2,2)
      real(dp):: uu,vv,amps2(2)
      complex(dp):: alpa0,alpa1,alpb0,alpb1,alpc0,alpc1
      complex(dp):: beta0,beta1,betb0,betb1,betc0,betc1
      real(dp)::
     & al1R,al1I,al2R,al2I,al3R,al3I,be1R,be1I,be2R,be2I,be3R,be3I
      real(dp)::
     & GGAlpha1_2a1re,GGAlpha1_2a1im,
     & GGAlpha2_2a1re,GGAlpha2_2a1im,
     & GGAlpha3_2a1re,GGAlpha3_2a1im
      real(dp)::
     & GGBeta1_2a1re,GGBeta1_2a1im,
     & GGBeta2_2a1re,GGBeta2_2a1im,
     & GGBeta3_2a1re,GGBeta3_2a1im
      real(dp)::
     & GGAlpha1_4a1re,GGAlpha1_4a1im,
     & GGAlpha2_4a1re,GGAlpha2_4a1im,
     & GGAlpha3_4a1re,GGAlpha3_4a1im
c--- statement function for hn of Eqs.(42), (43)
      hnpmm(s12,s23,s13,qsq,alpa0,alpb0,alpc0,alpa1,alpb1,alpc1)
     & =one/(eight*qsq*s12*s13)*(
     & +real(alpa0*conjg(alpa1),dp)*qsq*(s23*qsq+2*s12*s13)
     & +real(alpb0*conjg(alpb1),dp)*s23*(s12+s23)**2
     & +real(alpc0*conjg(alpc1),dp)*s23*(s13+s23)**2
     & -real((alpa0*conjg(alpb1)+alpa1*conjg(alpb0)),dp)
     &    *s23*(s12+s23)*qsq
     & +real((alpa0*conjg(alpc1)+alpa1*conjg(alpc0)),dp)
     &    *s23*(s13+s23)*qsq
     & -real((alpb0*conjg(alpc1)+alpb1*conjg(alpc0)),dp)
     &    *s23*(s23*qsq-s12*s13))
      hnppp(s12,s23,s13,qsq,beta0,betb0,betc0,beta1,betb1,betc1)
     & =one/(eight*qsq*s12)*(
     & +real(beta0*conjg(beta1),dp)*s13*(s12+s13)**2/s23
     & +real(betb0*conjg(betb1),dp)*s23*(s12+s23)**2/s13
     & +real(betc0*conjg(betc1),dp)*(2*s12*qsq+s13*s23)
     & -real((beta0*conjg(betb1)+beta1*conjg(betb0)),dp)
     &    *(s12*qsq-s13*s23)
     & +real((beta0*conjg(betc1)+beta1*conjg(betc0)),dp)
     &    *s13*(s12+s13)
     & +real((betb0*conjg(betc1)+betb1*conjg(betc0)),dp)
     &    *s23*(s12+s23))

c--- contribution only enters at NNLO
      if (order < 2) then
        amps2(:)=zip
        return
      endif

      if (1 == 2) then
c---- ! Check Becher et al. result at a single PS point
        s(p1,p2)=1._dp
        s(p1,p3)=-0.4_dp
        s(p4,p5)=0.1_dp**2
        s(p2,p3)=s(p4,p5)-s(p1,p2)-s(p1,p3)
        s(p2,p1)=s(p1,p2)
        s(p3,p1)=s(p1,p3)
        s(p3,p2)=s(p2,p3)
        scale=0.6_dp
        musq=scale**2
      endif

c      write(6,*) 'region',region
c      write(6,*) 's(p1,p2)',s(p1,p2)
c      write(6,*) 's(p1,p3)',s(p1,p3)
c      write(6,*) 's(p2,p3)',s(p2,p3)

      res2=zip
      amps(:,:,:,:)=czip
      do iperm=1,3
      if     (iperm == 1) then
        q1=p1
        q2=p2
        q3=p3
      elseif (iperm == 2) then
        q1=p2
        q2=p1
        q3=p3
      elseif (iperm == 3) then
        q1=p3
        q2=p2
        q3=p1
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
        write(6,*) 'Region 3 not anticipated in Zampgggsq'
        stop
        uu=-s(q2,q3)/s(q1,q3)
        vv=+s(p4,p5)/s(q1,q3)
      elseif (region == 4) then
        uu=-s(q1,q3)/s(q2,q3)
        vv=+s(p4,p5)/s(q2,q3)
      endif
      
c--- fill arrays for 2DHPLs
      call tdhpl(uu,vv,4,G1,G2,G3,G4,H1,H2,H3,H4)

c--- Coefficients for Region 2: need alpha and beta
      if (region == 2) then
        al1R=GGAlpha1_2a1re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4)
        al1I=GGAlpha1_2a1im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4)
        al2R=GGAlpha2_2a1re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4)
        al2I=GGAlpha2_2a1im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4)
        al3R=GGAlpha3_2a1re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4)
        al3I=GGAlpha3_2a1im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4)

        be1R=GGBeta1_2a1re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4)
        be1I=GGBeta1_2a1im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4)
        be2R=GGBeta2_2a1re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4)
        be2I=GGBeta2_2a1im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4)
        be3R=GGBeta3_2a1re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4)
        be3I=GGBeta3_2a1im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4)
      
c--- Coefficients for Region 4: only need alpha
      elseif (region == 4) then
        al1R=GGAlpha1_4a1re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4)
        al1I=GGAlpha1_4a1im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4)
        al2R=GGAlpha2_4a1re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4)
        al2I=GGAlpha2_4a1im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4)
        al3R=GGAlpha3_4a1re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4)
        al3I=GGAlpha3_4a1im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4)
      
      else
        write(6,*) 'Zampgggsq: region should be 2 or 4: ',region
        stop

      endif

      alpa1=cmplx(al1R,al1I,dp)
      alpb1=cmplx(al2R,al2I,dp)
      alpc1=cmplx(al3R,al3I,dp)
      beta1=cmplx(be1R,be1I,dp)
      betb1=cmplx(be2R,be2I,dp)
      betc1=cmplx(be3R,be3I,dp)

! result without Z decay
      if (iperm == 1) then
      res2=res2
     & +two*hnppp(s(q1,q2),s(q2,q3),s(q1,q3),s(p4,p5),
     &        beta1,betb1,betc1,beta1,betb1,betc1)
      endif  
      res2=res2
     & +two*hnpmm(s(q1,q2),s(q2,q3),s(q1,q3),s(p4,p5),
     &        alpa1,alpb1,alpc1,alpa1,alpb1,alpc1)
      
!     1st index is helicity of lepton line
!     2nd,3rd,4th indices are helicities of gluons
      if (iperm == 1) then
      amps(2,2,2,2)=ampgggppp(q1,q2,q3,p4,p5,beta1,betb1,betc1,za,zb)
      amps(1,2,2,2)=ampgggppp(q1,q2,q3,p5,p4,beta1,betb1,betc1,za,zb)
      amps(1,1,1,1)=ampgggppp(q1,q2,q3,p4,p5,beta1,betb1,betc1,zb,za)
      amps(2,1,1,1)=ampgggppp(q1,q2,q3,p5,p4,beta1,betb1,betc1,zb,za)

      amps(2,2,1,1)=ampgggpmm(q1,q2,q3,p4,p5,alpa1,alpb1,alpc1,za,zb)
      amps(1,2,1,1)=ampgggpmm(q1,q2,q3,p5,p4,alpa1,alpb1,alpc1,za,zb)
      amps(1,1,2,2)=ampgggpmm(q1,q2,q3,p4,p5,alpa1,alpb1,alpc1,zb,za)
      amps(2,1,2,2)=ampgggpmm(q1,q2,q3,p5,p4,alpa1,alpb1,alpc1,zb,za)
      endif

      if (iperm == 2) then
      amps(2,1,2,1)=ampgggpmm(q1,q2,q3,p4,p5,alpa1,alpb1,alpc1,za,zb)
      amps(1,1,2,1)=ampgggpmm(q1,q2,q3,p5,p4,alpa1,alpb1,alpc1,za,zb)
      amps(1,2,1,2)=ampgggpmm(q1,q2,q3,p4,p5,alpa1,alpb1,alpc1,zb,za)
      amps(2,2,1,2)=ampgggpmm(q1,q2,q3,p5,p4,alpa1,alpb1,alpc1,zb,za)
      endif
      
      if (iperm == 3) then
      amps(2,1,1,2)=ampgggpmm(q1,q2,q3,p4,p5,alpa1,alpb1,alpc1,za,zb)
      amps(1,1,1,2)=ampgggpmm(q1,q2,q3,p5,p4,alpa1,alpb1,alpc1,za,zb)
      amps(1,2,2,1)=ampgggpmm(q1,q2,q3,p4,p5,alpa1,alpb1,alpc1,zb,za)
      amps(2,2,2,1)=ampgggpmm(q1,q2,q3,p5,p4,alpa1,alpb1,alpc1,zb,za)
      endif
 
      enddo ! end of iperm loop
      
c---- This bit of code for testing squared MEs against Becher et al.
      msq2=res2
c      write(6,*) '2L/fourpi**2',msq2/fourpi**2
c      write(6,*) '2L not decayed',
c     &   (avegg*msq2*(40._dp/3._dp)*esq*gsq*ason4pi**2)
c     &  *(half*(2._dp*(L(2)+R(2))+3._dp*(L(1)+R(1))))**2
      
      do hl=1,2
        amps2(hl)=zip
        do ha=1,2
        do hb=1,2
        do hc=1,2
        amps2(hl)=amps2(hl)
     &   +real(amps(hl,ha,hb,hc)*conjg(amps(hl,ha,hb,hc)),dp)
c        write(6,*) hl,ha,hb,hc,amps(hl,ha,hb,hc)
        enddo
        enddo
        enddo
      enddo
      
c--- convert from coefficient of [as/4/pi]^n to coefficient of [as/2/pi]^n
      amps2(:)=amps2(:)/four

      return
      end


