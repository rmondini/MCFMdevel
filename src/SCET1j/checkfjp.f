      subroutine checkfjp(p)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'nf.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'zcouple.f'
      include 'scale.f'
      real(dp):: p(mxpart,4),aemmz,dot
      common/em/aemmz
      
c--- point from FJP ! interchange 1<->2 and 3<->4
      p(1,4)=150.00000000000000_dp
      p(1,1)=0.0000000000000000_dp
      p(1,2)=0.0000000000000000_dp
      p(1,3)=150.00000000000000_dp
      p(2,4)=150.00000000000000_dp
      p(2,1)=0.0000000000000000_dp
      p(2,2)=0.0000000000000000_dp
      p(2,3)=-150.00000000000000_dp
      p(5,4)=-136.14136934373334_dp
      p(5,1)=58.950942178836499_dp
      p(5,2)=-102.10602700780002_dp
      p(5,3)=-68.070684671866687_dp
      p(3,4)=24.412283644899802_dp
      p(3,1)=17.262091309558095_dp
      p(3,2)=0.0000000000000000_dp
      p(3,3)=17.262091309558098_dp
      p(4,4)=139.44634701136687_dp
      p(4,1)=41.688850869278404_dp
      p(4,2)=-102.10602700780002_dp
      p(4,3)=-85.332775981424788_dp
c---- crossing to MCFM notation
      p(1,:)=-p(1,:)
      p(2,:)=-p(2,:)
      p(5,:)=-p(5,:)
c---- small change in alpha-s
      as=0.1180024244524787_dp
      gsq=fourpi*as
      ason2pi=as/twopi
      ason4pi=as/fourpi
      musq=2._dp*dot(p,3,4)
      scale=sqrt(musq)

c--- EW parameters
      xw=0.2226459_dp
      aemmz=1._dp/132.3384323_dp
      gwsq=fourpi*aemmz/xw
      esq=gwsq*xw
      gw=sqrt(gwsq)
      call couplz(xw)
      q1=qe
      l1=le
      r1=re
      
      return
      end
      
