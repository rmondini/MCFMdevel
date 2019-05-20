!#################################################################
!##### O(as) and O(as^2) soft contribution to H->bbg type of
!##### processes. Calculates soft1(-1:1) and soft2(-1:3), i.e.
!##### contributions to delta(tau1), L0(tau1), L1(tau1), and so on.
!##### Factor of alphas/(2*pi) and (alphas/(2*pi))**2 have been 
!##### extracted out.
!#################################################################
!---  RM, Nov 18
      
      subroutine soft_dec_3j(order,y12,y13,y23,soft1,soft2)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'scale.f'
      include 'nf.f'
      include 'scet_const.f'
      integer,intent(in)::order
      real(dp),intent(in)::y12,y13,y23
      real(dp),intent(out)::soft1(-1:1),soft2(-1:3)
      real(dp)::Iijm(6),Lss,Lsd,Isum12,Isum3456
      real(dp)::soft2ab(-1:3),soft2nab(-1:3)
      real(dp)::c000,c010,c020,c030,c001,c011,c021,c031,
     & c002,c012,c022,c032,c003,c013,c023,c033,
     & c100,c110,c120,c130,c101,c111,c121,c131,
     & c102,c112,c122,c132,c103,c113,c123,c133,
     & c200,c210,c220,c230,c201,c211,c221,c231,
     & c202,c212,c222,c232,c203,c213,c223,c233,
     & c300,c310,c320,c330,c301,c311,c321,c331,
     & c302,c312,c322,c332,c303,c313,c323,c333

      soft1(:)=zip
      soft2(:)=zip
      soft2ab(:)=zip
      soft2nab(:)=zip

      if (order < 1) return

      call computeIijmnobug(y12,y23,y13,Iijm)
      Lss=log(y13)+log(y23)
      Lsd=log(y13)-log(y23)
      Isum12=Iijm(1)+Iijm(2)
      Isum3456=Iijm(3)+Iijm(4)+Iijm(5)+Iijm(6)

!---- contribution to delta(tau)
      soft1(-1) = CA*(-0.25_dp*(Lsd**2+Lss**2)+Isum12-Isum3456
     & +half*zeta2+half*log(y12)**2)
     & +CF*(-two*Isum12+zeta2-log(y12)**2)
!---- contribution to L0(tau)
      soft1(0) = two*CA*(Lss-log(y12))+four*log(y12)*CF
!---- contribution to L1(tau)
      soft1(1) = -four*CA-eight*CF

      if (order < 2) return

!---- abelian contribution
      soft2ab(3)=8*(CA+2*CF)**2
      soft2ab(2)=-12*(CA+2*CF)*(CA*(log(y13)+log(y23)-log(y12))
     & +2*log(y12)*CF)

      soft2ab(1) =  + CA**2 * (
     &     - 18._dp*zeta2
     &     + 6._dp*log(y13)**2
     &     - 8._dp*log(y13)*log(y12)
     &     + 8._dp*log(y13)*log(y23)
     &     + 2._dp*log(y12)**2
     &     - 8._dp*log(y12)*log(y23)
     &     + 6._dp*log(y23)**2
     &     + 4._dp*Iijm(6)
     &     + 4._dp*Iijm(5)
     &     + 4._dp*Iijm(3)
     &     + 4._dp*Iijm(4)
     &     - 4._dp*Iijm(2)
     &     - 4._dp*Iijm(1)
     &     )
      soft2ab(1) = soft2ab(1) + CF*CA * (
     &     - 72._dp*zeta2
     &     + 4._dp*log(y13)**2
     &     + 16._dp*log(y13)*log(y12)
     &     - 16._dp*log(y12)**2
     &     + 16._dp*log(y12)*log(y23)
     &     + 4._dp*log(y23)**2
     &     + 8._dp*Iijm(6)
     &     + 8._dp*Iijm(5)
     &     + 8._dp*Iijm(3)
     &     + 8._dp*Iijm(4)
     &     )
      soft2ab(1) = soft2ab(1) + CF**2 * (
     &     - 72._dp*zeta2
     &     + 24._dp*log(y12)**2
     &     + 16._dp*Iijm(2)
     &     + 16._dp*Iijm(1)
     &     )

      soft2ab(0) =  + CA**2 * (
     &     + 16._dp*zeta3
     &     + 9._dp*log(y13)*zeta2
     &     - log(y13)**3
     &     + log(y13)**2*log(y12)
     &     - log(y13)**2*log(y23)
     &     + log(y13)*log(y12)**2
     &     - log(y13)*log(y23)**2
     &     - 2._dp*log(y13)*Iijm(6)
     &     - 2._dp*log(y13)*Iijm(5)
     &     - 2._dp*log(y13)*Iijm(3)
     &     - 2._dp*log(y13)*Iijm(4)
     &     + 2._dp*log(y13)*Iijm(2)
     &     + 2._dp*log(y13)*Iijm(1)
     &     - 9._dp*log(y12)*zeta2
     &     )
      soft2ab(0) = soft2ab(0) + CA**2 * (
     &     - log(y12)**3
     &     + log(y12)**2*log(y23)
     &     + log(y12)*log(y23)**2
     &     + 2._dp*log(y12)*Iijm(6)
     &     + 2._dp*log(y12)*Iijm(5)
     &     + 2._dp*log(y12)*Iijm(3)
     &     + 2._dp*log(y12)*Iijm(4)
     &     - 2._dp*log(y12)*Iijm(2)
     &     - 2._dp*log(y12)*Iijm(1)
     &     + 9._dp*log(y23)*zeta2
     &     - log(y23)**3
     &     - 2._dp*log(y23)*Iijm(6)
     &     - 2._dp*log(y23)*Iijm(5)
     &     - 2._dp*log(y23)*Iijm(3)
     &     )
      soft2ab(0) = soft2ab(0) + CA**2 * (
     &     - 2._dp*log(y23)*Iijm(4)
     &     + 2._dp*log(y23)*Iijm(2)
     &     + 2._dp*log(y23)*Iijm(1)
     &     )
      soft2ab(0) = soft2ab(0) + CF*CA * (
     &     + 64._dp*zeta3
     &     + 18._dp*log(y13)*zeta2
     &     - 2._dp*log(y13)**2*log(y12)
     &     - 2._dp*log(y13)*log(y12)**2
     &     - 4._dp*log(y13)*Iijm(2)
     &     - 4._dp*log(y13)*Iijm(1)
     &     + 4._dp*log(y12)**3
     &     - 2._dp*log(y12)**2*log(y23)
     &     - 2._dp*log(y12)*log(y23)**2
     &     - 4._dp*log(y12)*Iijm(6)
     &     - 4._dp*log(y12)*Iijm(5)
     &     - 4._dp*log(y12)*Iijm(3)
     &     - 4._dp*log(y12)*Iijm(4)
     &     + 8._dp*log(y12)*Iijm(2)
     &     )
      soft2ab(0) = soft2ab(0) + CF*CA * (
     &     + 8._dp*log(y12)*Iijm(1)
     &     + 18._dp*log(y23)*zeta2
     &     - 4._dp*log(y23)*Iijm(2)
     &     - 4._dp*log(y23)*Iijm(1)
     &     )
      soft2ab(0) = soft2ab(0) + CF**2 * (
     &     + 64._dp*zeta3
     &     + 36._dp*log(y12)*zeta2
     &     - 4._dp*log(y12)**3
     &     - 8._dp*log(y12)*Iijm(2)
     &     - 8._dp*log(y12)*Iijm(1)
     &     )

      soft2ab(-1) =  + CA**2 * (
     &     - 27._dp/40._dp*zeta2**2
     &     - 8._dp*log(y13)*zeta3
     &     - 9._dp/4._dp*log(y13)**2*zeta2
     &     + 1._dp/8._dp*log(y13)**4
     &     - 1._dp/4._dp*log(y13)**2*log(y12)**2
     &     + 1._dp/4._dp*log(y13)**2*log(y23)**2
     &     + 1._dp/2._dp*log(y13)**2*Iijm(6)
     &     + 1._dp/2._dp*log(y13)**2*Iijm(5)
     &     + 1._dp/2._dp*log(y13)**2*Iijm(3)
     &     + 1._dp/2._dp*log(y13)**2*Iijm(4)
     &     - 1._dp/2._dp*log(y13)**2*Iijm(2)
     &     - 1._dp/2._dp*log(y13)**2*Iijm(1)
     &     + 4._dp*log(y13)*log(y12)*zeta2
     &     - 4._dp*log(y13)*log(y23)*zeta2
     &     )
      soft2ab(-1) = soft2ab(-1) + CA**2 * (
     &     + 8._dp*log(y12)*zeta3
     &     - 7._dp/4._dp*log(y12)**2*zeta2
     &     + 1._dp/8._dp*log(y12)**4
     &     - 1._dp/4._dp*log(y12)**2*log(y23)**2
     &     - 1._dp/2._dp*log(y12)**2*Iijm(6)
     &     - 1._dp/2._dp*log(y12)**2*Iijm(5)
     &     - 1._dp/2._dp*log(y12)**2*Iijm(3)
     &     - 1._dp/2._dp*log(y12)**2*Iijm(4)
     &     + 1._dp/2._dp*log(y12)**2*Iijm(2)
     &     + 1._dp/2._dp*log(y12)**2*Iijm(1)
     &     + 4._dp*log(y12)*log(y23)*zeta2
     &     - 8._dp*log(y23)*zeta3
     &     - 9._dp/4._dp*log(y23)**2*zeta2
     &     + 1._dp/8._dp*log(y23)**4
     &     )
      soft2ab(-1) = soft2ab(-1) + CA**2 * (
     &     + 1._dp/2._dp*log(y23)**2*Iijm(6)
     &     + 1._dp/2._dp*log(y23)**2*Iijm(5)
     &     + 1._dp/2._dp*log(y23)**2*Iijm(3)
     &     + 1._dp/2._dp*log(y23)**2*Iijm(4)
     &     - 1._dp/2._dp*log(y23)**2*Iijm(2)
     &     - 1._dp/2._dp*log(y23)**2*Iijm(1)
     &     - 1._dp/2._dp*Iijm(6)*zeta2
     &     + 1._dp/2._dp*Iijm(6)**2
     &     + Iijm(6)*Iijm(5)
     &     + Iijm(6)*Iijm(3)
     &     + Iijm(6)*Iijm(4)
     &     - Iijm(6)*Iijm(2)
     &     - Iijm(6)*Iijm(1)
     &     - 1._dp/2._dp*Iijm(5)*zeta2
     &     )
      soft2ab(-1) = soft2ab(-1) + CA**2 * (
     &     + 1._dp/2._dp*Iijm(5)**2
     &     + Iijm(5)*Iijm(3)
     &     + Iijm(5)*Iijm(4)
     &     - Iijm(5)*Iijm(2)
     &     - Iijm(5)*Iijm(1)
     &     - 1._dp/2._dp*Iijm(3)*zeta2
     &     + 1._dp/2._dp*Iijm(3)**2
     &     + Iijm(3)*Iijm(4)
     &     - Iijm(3)*Iijm(2)
     &     - Iijm(3)*Iijm(1)
     &     - 1._dp/2._dp*Iijm(4)*zeta2
     &     + 1._dp/2._dp*Iijm(4)**2
     &     - Iijm(4)*Iijm(2)
     &     - Iijm(4)*Iijm(1)
     &     )
      soft2ab(-1) = soft2ab(-1) + CA**2 * (
     &     + 1._dp/2._dp*Iijm(2)*zeta2
     &     + 1._dp/2._dp*Iijm(2)**2
     &     + Iijm(2)*Iijm(1)
     &     + 1._dp/2._dp*Iijm(1)*zeta2
     &     + 1._dp/2._dp*Iijm(1)**2
     &     )
      soft2ab(-1) = soft2ab(-1) + CF*CA * (
     &     - 27._dp/10._dp*zeta2**2
     &     - 16._dp*log(y13)*zeta3
     &     - 1._dp/2._dp*log(y13)**2*zeta2
     &     + 1._dp/2._dp*log(y13)**2*log(y12)**2
     &     + log(y13)**2*Iijm(2)
     &     + log(y13)**2*Iijm(1)
     &     - 8._dp*log(y13)*log(y12)*zeta2
     &     + 8._dp*log(y12)**2*zeta2
     &     - 1._dp/2._dp*log(y12)**4
     &     + 1._dp/2._dp*log(y12)**2*log(y23)**2
     &     + log(y12)**2*Iijm(6)
     &     + log(y12)**2*Iijm(5)
     &     + log(y12)**2*Iijm(3)
     &     + log(y12)**2*Iijm(4)
     &     )
      soft2ab(-1) = soft2ab(-1) + CF*CA * (
     &     - 2._dp*log(y12)**2*Iijm(2)
     &     - 2._dp*log(y12)**2*Iijm(1)
     &     - 8._dp*log(y12)*log(y23)*zeta2
     &     - 16._dp*log(y23)*zeta3
     &     - 1._dp/2._dp*log(y23)**2*zeta2
     &     + log(y23)**2*Iijm(2)
     &     + log(y23)**2*Iijm(1)
     &     - Iijm(6)*zeta2
     &     + 2._dp*Iijm(6)*Iijm(2)
     &     + 2._dp*Iijm(6)*Iijm(1)
     &     - Iijm(5)*zeta2
     &     + 2._dp*Iijm(5)*Iijm(2)
     &     + 2._dp*Iijm(5)*Iijm(1)
     &     - Iijm(3)*zeta2
     &     )
      soft2ab(-1) = soft2ab(-1) + CF*CA * (
     &     + 2._dp*Iijm(3)*Iijm(2)
     &     + 2._dp*Iijm(3)*Iijm(1)
     &     - Iijm(4)*zeta2
     &     + 2._dp*Iijm(4)*Iijm(2)
     &     + 2._dp*Iijm(4)*Iijm(1)
     &     - 2._dp*Iijm(2)**2
     &     - 4._dp*Iijm(2)*Iijm(1)
     &     - 2._dp*Iijm(1)**2
     &     )
      soft2ab(-1) = soft2ab(-1) + CF**2 * (
     &     - 27._dp/10._dp*zeta2**2
     &     - 32._dp*log(y12)*zeta3
     &     - 9._dp*log(y12)**2*zeta2
     &     + 1._dp/2._dp*log(y12)**4
     &     + 2._dp*log(y12)**2*Iijm(2)
     &     + 2._dp*log(y12)**2*Iijm(1)
     &     - 2._dp*Iijm(2)*zeta2
     &     + 2._dp*Iijm(2)**2
     &     + 4._dp*Iijm(2)*Iijm(1)
     &     - 2._dp*Iijm(1)*zeta2
     &     + 2._dp*Iijm(1)**2
     &     )

!---- non-abelian contribution
      soft2nab(3)=zip
      soft2nab(2)=2*be0*(CA+2*CF)
      soft2nab(1)=+CA*(-half*Ga1
     & -two*log(y13)*be0+two*log(y12)*be0-two*log(y23)*be0)
     & +CF*(-Ga1-four*log(y12)*be0)
      soft2nab(0)=CA*(
     &     - quarter*gams1
     &     - half*zeta2*be0
     &     + quarter*log(y13)*Ga1
     &     + half*log(y13)**2*be0
     &     - quarter*log(y12)*Ga1
     &     - half*log(y12)**2*be0
     &     + quarter*log(y23)*Ga1
     &     + half*log(y23)**2*be0
     &     + Iijm(5)*be0
     &     + Iijm(6)*be0
     &     + Iijm(3)*be0
     &     + Iijm(4)*be0
     &     - Iijm(1)*be0
     &     - Iijm(2)*be0
     &     )
      soft2nab(0) = soft2nab(0) + CF * (
     &     - half*gams1
     &     - zeta2*be0
     &     + half*log(y12)*Ga1
     &     + log(y12)**2*be0
     &     + two*Iijm(1)*be0
     &     + two*Iijm(2)*be0
     &     )

!---- CEMW fit 1711.09984
      c000=40.36899351167569_dp
      c001=25.62789055081602_dp
      c002=-8.736690397854694_dp
      c003=-2.114404899583132_dp
      c010=26.13736916163353_dp
      c011=2.689865691050357_dp
      c012=-0.5515149013106969_dp
      c013=0.5002962443905451_dp
      c020=-8.592225292354909_dp
      c021=-2.371774585802533_dp
      c022=1.706294721953703_dp
      c023=0.2791885713486318_dp
      c030=-2.108712689473365_dp
      c031=0.04554347592404617_dp
      c032=0.5960886850291142_dp
      c033=-0.03789129411205423_dp
      c100=14.67940638634606_dp
      c101=-1.786642833784238_dp
      c102=-0.8376186408480450_dp
      c103=-0.7671773626451408_dp
      c110=-1.815173614514120_dp
      c111=-7.383702376905395_dp
      c112=-5.365482713615537_dp
      c113=-0.8654462906420545_dp
      c120=0.2341817689090819_dp
      c121=-9.647937031219083_dp
      c122=-5.753210482964973_dp
      c123=-0.9288308952872950_dp
      c130=-0.4315847674561632_dp
      c131=-1.990762803587165_dp
      c132=-0.6488570120772139_dp
      c133=-0.1212740605786565_dp
      c200=-2.441158778463018_dp
      c201=-0.8473608891348289_dp
      c202=-2.235994141288788_dp
      c203=-0.8551720435319591_dp
      c210=0.3306300397987830_dp
      c211=-6.081268122109750_dp
      c212=-8.082000043182264_dp
      c213=-2.039532189769125_dp
      c220=2.430077237125532_dp
      c221=-4.096188898362771_dp
      c222=-5.362694355097635_dp
      c223=-1.259591065936183_dp
      c230=0.4283812380848217_dp
      c231=-0.6952642925422763_dp
      c232=-0.6236594120772266_dp
      c233=-0.1391711344196455_dp
      c300=0.04761734323240323_dp
      c301=-0.2223776441027224_dp
      c302=-0.6667735807177795_dp
      c303=-0.3698205811417911_dp
      c310=0.1603960831737783_dp
      c311=-0.9843435058271673_dp
      c312=-1.570616193506059_dp
      c313=-0.5449941237467155_dp
      c320=0.8608741841650420_dp
      c321=0.2488227601645385_dp
      c322=-0.5170169998270092_dp
      c323=-0.2045839157158000_dp
      c330=0.01511441241649222_dp
      c331=-0.05222533897401384_dp
      c332=-0.06296920559206375_dp
      c333=-0.01947806237908754_dp

      soft2nab(-1)= c000 + c100*log(y12) + c200*log(y12)**2 + 
     &  c300*log(y12)**3 + c010*log(y13) + 
     &  c110*log(y12)*log(y13) + c210*log(y12)**2*log(y13) + 
     &  c310*log(y12)**3*log(y13) + c020*log(y13)**2 + 
     &  c120*log(y12)*log(y13)**2 + c220*log(y12)**2*log(y13)**2 + 
     &  c320*log(y12)**3*log(y13)**2 + c030*log(y13)**3 + 
     &  c130*log(y12)*log(y13)**3 + c230*log(y12)**2*log(y13)**3 + 
     &  c330*log(y12)**3*log(y13)**3 + c001*log(y23) + 
     &  c101*log(y12)*log(y23) + c201*log(y12)**2*log(y23) + 
     &  c301*log(y12)**3*log(y23) + c011*log(y13)*log(y23) + 
     &  c111*log(y12)*log(y13)*log(y23) + 
     &  c211*log(y12)**2*log(y13)*log(y23) + 
     &  c311*log(y12)**3*log(y13)*log(y23) + 
     &  c021*log(y13)**2*log(y23) + 
     &  c121*log(y12)*log(y13)**2*log(y23) + 
     &  c221*log(y12)**2*log(y13)**2*log(y23) + 
     &  c321*log(y12)**3*log(y13)**2*log(y23) + 
     &  c031*log(y13)**3*log(y23) + 
     &  c131*log(y12)*log(y13)**3*log(y23) + 
     &  c231*log(y12)**2*log(y13)**3*log(y23) + 
     &  c331*log(y12)**3*log(y13)**3*log(y23) + c002*log(y23)**2 + 
     &  c102*log(y12)*log(y23)**2 + c202*log(y12)**2*log(y23)**2 + 
     &  c302*log(y12)**3*log(y23)**2 + c012*log(y13)*log(y23)**2 + 
     &  c112*log(y12)*log(y13)*log(y23)**2 + 
     &  c212*log(y12)**2*log(y13)*log(y23)**2 + 
     &  c312*log(y12)**3*log(y13)*log(y23)**2 + 
     &  c022*log(y13)**2*log(y23)**2 + 
     &  c122*log(y12)*log(y13)**2*log(y23)**2 + 
     &  c222*log(y12)**2*log(y13)**2*log(y23)**2 + 
     &  c322*log(y12)**3*log(y13)**2*log(y23)**2 + 
     &  c032*log(y13)**3*log(y23)**2 + 
     &  c132*log(y12)*log(y13)**3*log(y23)**2 + 
     &  c232*log(y12)**2*log(y13)**3*log(y23)**2 + 
     &  c332*log(y12)**3*log(y13)**3*log(y23)**2 + 
     &  c003*log(y23)**3 + c103*log(y12)*log(y23)**3 + 
     &  c203*log(y12)**2*log(y23)**3 + 
     &  c303*log(y12)**3*log(y23)**3 + c013*log(y13)*log(y23)**3 + 
     &  c113*log(y12)*log(y13)*log(y23)**3 + 
     &  c213*log(y12)**2*log(y13)*log(y23)**3 + 
     &  c313*log(y12)**3*log(y13)*log(y23)**3 + 
     &  c023*log(y13)**2*log(y23)**3 + 
     &  c123*log(y12)*log(y13)**2*log(y23)**3 + 
     &  c223*log(y12)**2*log(y13)**2*log(y23)**3 + 
     &  c323*log(y12)**3*log(y13)**2*log(y23)**3 + 
     &  c033*log(y13)**3*log(y23)**3 + 
     &  c133*log(y12)*log(y13)**3*log(y23)**3 + 
     &  c233*log(y12)**2*log(y13)**3*log(y23)**3 + 
     &  c333*log(y12)**3*log(y13)**3*log(y23)**3 

!---- contribution to delta(tau)
      soft2(-1) = soft2ab(-1)+soft2nab(-1)
!---- contribution to L0(tau)
      soft2(0) = soft2ab(0)+soft2nab(0)
!---- contribution to L1(tau)
      soft2(1) = soft2ab(1)+soft2nab(1)
!---- contribution to L2(tau)
      soft2(2) = soft2ab(2)+soft2nab(2)
!---- contribution to L3(tau)
      soft2(3) = soft2ab(3)+soft2nab(3)

!      soft2(-1) = soft2nab(-1) ! debug
!      soft2(0) = soft2nab(0) ! debug
!      soft2(1) = soft2nab(1) ! debug
!      soft2(2) = soft2nab(2) ! debug
!      soft2(3) = soft2nab(3) ! debug

      return
      end

