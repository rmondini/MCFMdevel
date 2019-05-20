      subroutine lumxmsq_hb(p,xx,z1,z2,QB,order,xmsq)
      implicit none
      include 'types.f'
c---- H+j production
c---- yb**2 term only (i.e. the crossing of H->bbg decay)
c     b(-p1) + g(-p2) --> H(p3+p4) + b(p5)
c     b~(-p1) + g(-p2) --> H(p3+p4) + b~(p5)
c     g(-p1) + b(-p2) --> H(p3+p4) + b(p5)
c     g(-p1) + b~(-p2) --> H(p3+p4) + b~(p5)
c
c     b(-p1) + b~(-p2) --> H(p3+p4) + g(p5)
c     b~(-p1) + b(-p2) --> H(p3+p4) + g(p5)
c---- RM, May 2019
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'sprods_com.f'
      include 'scale.f'
      include 'facscale.f'
      include 'ewcharge.f'
      include 'hdecaymode.f'
      include 'scet_const.f'
      include 'hbbparams.f'
      include 'taucut.f'
      integer:: j,k,ih1,ih2,order
      real(dp),intent(out):: xmsq
      real(dp) :: xvarqg,yvarqg,zvarqg,uvarqg,vvarqg,
     & xvargq,yvargq,zvargq,uvargq,vvargq,
     & xvarqa,yvarqa,zvarqa,uvarqa,vvarqa,
     & y12,y15,y25,s34
      real(dp):: p(mxpart,4),xx(2),hardqg(2),hardgq(2),hardqa(2),
     & Iijm(6),
     & soft1qg(-1:1),soft2qg(-1:3),soft2qg_nab(-1:3),
     & soft1gq(-1:1),soft2gq(-1:3),soft2gq_nab(-1:3),
     & soft1qa(-1:1),soft2qa(-1:3),soft2qa_nab(-1:3),
     & beama0(-5:5),beamb0(-5:5),tauc,getdynamictau,
     & beama1(-5:5,-1:1),beamb1(-5:5,-1:1),
     & beama2(-5:5,-1:3),beamb2(-5:5,-1:3),
     & jet1q(-1:1),jet2q(-1:3),jet1g(-1:1),jet2g(-1:3),
     & z1,z2,QB(2),bitqg,bitgq,lobitqg,lobitgq,
     & bitqa,lobitqa,assemblejet
      logical inclqa
      common/density/ih1,ih2

!---------------------
!==== whether to include bb~->Hg and b~b->Hg
      inclqa=.true.
!---------------------

C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C~~~~~~~~~~~~~~~~~ the Z1jet process has mb_eff=0 ~~~~~~~~~~~~~~~~~~~
      mb_eff=2.95_dp
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
c---- debug
!      write(*,*) '~~~the beginning~~~'
!      order=2
c---- debug

      call dotem(5,p,s)

      s34=(p(3,4)+p(4,4))**2
     &   -(p(3,1)+p(4,1))**2-(p(3,2)+p(4,2))**2-(p(3,3)+p(4,3))**2

c---- debug
!      write(*,*) 'p(1)=',p(1,:)
!      write(*,*) 'p(2)=',p(2,:)
!      write(*,*) 'p(3)=',p(3,:)
!      write(*,*) 'p(4)=',p(4,:)
!      write(*,*) 'p(5)=',p(5,:)
!      write(*,*) 's(1,1)=',s(1,1)
!      write(*,*) 's(2,2)=',s(2,2)
!      write(*,*) 's(3,3)=',s(3,3)
!      write(*,*) 's(4,4)=',s(4,4)
!      write(*,*) 's(5,5)=',s(5,5)
!      write(*,*) 's(3,4)=',s(3,4)
!      write(*,*) 's34=',s34
!      write(*,*) 'xx=',xx
!      write(*,*) 'z1=',z1
!      write(*,*) 'z2=',z2
!      write(*,*) 'QB=',QB
!      write(*,*) 'order=',order
!      write(*,*) 's(1,2)=',s(1,2)
!      write(*,*) 's(1,5)=',s(1,5)
c---- debug

c---- beam functions
      if (order >= 0) then
      call fdist_select(ih1,1,xx(1),facscale,beama0)
      call fdist_select(ih2,2,xx(2),facscale,beamb0)
      endif
      if (order >= 1) then
      call xbeam1bis(ih1,1,z1,xx(1),QB(1),beama1)
      call xbeam1bis(ih2,2,z2,xx(2),QB(2),beamb1)
      endif
      if (order >= 2) then
      call xbeam2bis(ih1,1,z1,xx(1),QB(1),beama2)
      call xbeam2bis(ih2,2,z2,xx(2),QB(2),beamb2)
      endif

c---- debug
!      write(*,*) 'beama1=',beama1
!      write(*,*) 'beama2=',beama2
c---- debug

c---- jet functions
      call jetq(order,two*p(5,4),jet1q,jet2q)
      call jetg(order,two*p(5,4),jet1g,jet2g)

c---- debug
!      write(*,*) 'jet1q=',jet1q
!      write(*,*) 'jet2q=',jet2q
!      write(*,*) 'jet1g=',jet1g
!      write(*,*) 'jet2g=',jet2g
c---- debug

c---- soft functions
      y12=s(1,2)/p(1,4)/p(2,4)/four
      y15=s(1,5)/p(1,4)/p(5,4)/four
      y25=s(2,5)/p(2,4)/p(5,4)/four

c---- debug
!      write(*,*) 'y12=',y12
!      write(*,*) 'y15=',y15
!      write(*,*) 'y25=',y25
!      y15=0.80_dp/two
!      y25=one-y15
c---- debug

c---- option 1
c      call soft_qgq(order,y12,y15,y25,soft1qg,soft2qg)
c      call soft_qgq(order,y12,y25,y15,soft1gq,soft2gq)
c      not implemented for soft_qag
c---- option 2
      call computeIijmnobug(y12,y25,y15,Iijm)
      call soft_ab_qgq(order,y12,y25,y15,Iijm,1,2,3,4,5,6,soft1qg,soft2qg)
      call soft_ab_qgq(order,y12,y15,y25,Iijm,2,1,6,5,4,3,soft1gq,soft2gq)
      call soft_ab_qag(order,y12,y25,y15,Iijm,1,2,3,4,5,6,soft1qa,soft2qa)
      if (order > 1) then
         call soft_nab_qgq(order,y12,y25,y15,Iijm,1,2,3,4,5,6,soft2qg_nab)
         call soft_nab_qgq(order,y12,y15,y25,Iijm,2,1,6,5,4,3,soft2gq_nab)
         call soft_nab_qag(order,y12,y25,y15,Iijm,1,2,3,4,5,6,soft2qa_nab)
         soft2qg(:)=soft2qg(:)+soft2qg_nab(:)
         soft2gq(:)=soft2gq(:)+soft2gq_nab(:)
         soft2qa(:)=soft2qa(:)+soft2qa_nab(:)
      endif

c---- debug
!      write(*,*) 'soft1qg=',soft1qg
!      write(*,*) 'soft2qg=',soft2qg
!      write(*,*) 'soft1gq=',soft1gq
!      write(*,*) 'soft2gq=',soft2gq
!      write(*,*) 'soft1qa=',soft1qa
!      write(*,*) 'soft2qa=',soft2qa
c---- debug

c---- hard function
      xvarqg=s(1,5)/s34
      yvarqg=s(2,5)/s34
      zvarqg=s(1,2)/s34
      uvarqg=-yvarqg/zvarqg
      vvarqg=one/zvarqg

      xvargq=s(2,5)/s34
      yvargq=s(1,5)/s34
      zvargq=s(1,2)/s34
      uvargq=-yvargq/zvargq
      vvargq=one/zvargq

      xvarqa=s(1,2)/s34
      yvarqa=s(1,5)/s34
      zvarqa=s(2,5)/s34
      uvarqa=-yvarqa/xvarqa
      vvarqa=one/xvarqa

c---- debug
!      uvarqa=0.0800000641679_dp
!      vvarqa=0.359999861089_dp
!      musq=15000_dp
c---- debug

c---- debug
!      write(*,*) 'hmass=',hmass
!      write(*,*) 'musq',musq
!      write(*,*) 'xvarqg=',xvarqg
!      write(*,*) 'yvarqg=',yvarqg
!      write(*,*) 'zvarqg=',zvarqg
!      write(*,*) 'uvarqg=',uvarqg
!      write(*,*) 'vvarqg=',vvarqg
!      write(*,*)
!      write(*,*) 'xvargq=',xvargq
!      write(*,*) 'yvargq=',yvargq
!      write(*,*) 'zvargq=',zvargq
!      write(*,*) 'uvargq=',uvargq
!      write(*,*) 'vvargq=',vvargq
!      write(*,*)
!      write(*,*) 'xvarqa=',xvarqa
!      write(*,*) 'yvarqa=',yvarqa
!      write(*,*) 'zvarqa=',zvarqa
!      write(*,*) 'uvarqa=',uvarqa
!      write(*,*) 'vvarqa=',vvarqa
c---- debug

      call hard_hbbg_crossed(order,2,uvarqg,vvarqg,musq,hardqg)
      call hard_hbbg_crossed(order,2,uvargq,vvargq,musq,hardgq)

      if(inclqa) then
         call hard_hbbg_crossed(order,1,uvarqa,vvarqa,musq,hardqa)
      endif

c---- debug
!      write(*,*) 'hardqg=',hardqg
!      write(*,*) 'hardgq=',hardgq
!      write(*,*) 'hardqa=',hardqa
c---- debug

!==== LO bg->Hb msq
!==== -1 is from crossing from Hbbg
!==== avg over spin and color
      lobitqg=32._dp*pi**2*gwsq*mb_eff**2/(four*wmass**2)*(xvarqg**2+one)/(yvarqg*zvarqg)
      lobitqg=-ason2pi*lobitqg*xn*cf*aveqg
!==== LO gb->Hb msq
!==== -1 is from crossing from Hbbg
!==== avg over spin and color
      lobitgq=32._dp*pi**2*gwsq*mb_eff**2/(four*wmass**2)*(xvargq**2+one)/(yvargq*zvargq)
      lobitgq=-ason2pi*lobitgq*xn*cf*aveqg
!==== LO bb~->Hb msq
!==== avg over spin and color
      lobitqa=32._dp*pi**2*gwsq*mb_eff**2/(four*wmass**2)*(xvarqa**2+one)/(yvarqa*zvarqa)
      lobitqa=ason2pi*lobitqa*xn*cf*aveqq

c---- debug
!      write(*,*) 'gwsq=',gwsq
!      write(*,*) 'mb_eff=',mb_eff
!      write(*,*) 'wmass=',wmass
!      write(*,*) 'ason2pi=',ason2pi
!      write(*,*) 'lobitqg=',lobitqg
!      write(*,*) 'lobitgq=',lobitgq
!      write(*,*) 'xn=',xn
!      write(*,*) 'cf=',cf
!      write(*,*) 'aveqg=',aveqg
c---- debug

c---- assemble below-cut piece

      if (dynamictau) then
        tauc=getdynamictau(p)
      else
        tauc=taucut
      endif

c---- debug
!      write(*,*) 'dynamictau=',dynamictau
!      write(*,*) 'tauc=',taucut
c---- debug

!==== bg
      bitqg=assemblejet(order,tauc,
     &     beama0(5),beamb0(0),beama1(5,:),beamb1(0,:),
     &     beama2(5,:),beamb2(0,:),soft1qg,soft2qg,jet1q,jet2q,hardqg)
!==== b~g
      bitqg=bitqg+assemblejet(order,tauc,
     &     beama0(-5),beamb0(0),beama1(-5,:),beamb1(0,:),
     &     beama2(-5,:),beamb2(0,:),soft1qg,soft2qg,jet1q,jet2q,hardqg)
!==== gb
      bitgq=assemblejet(order,tauc,
     &     beama0(0),beamb0(5),beama1(0,:),beamb1(5,:),
     &     beama2(0,:),beamb2(5,:),soft1gq,soft2gq,jet1q,jet2q,hardgq)
!==== gb~
      bitgq=bitgq+assemblejet(order,tauc,
     &     beama0(0),beamb0(-5),beama1(0,:),beamb1(-5,:),
     &     beama2(0,:),beamb2(-5,:),soft1gq,soft2gq,jet1q,jet2q,hardgq)

      bitqg=bitqg*lobitqg
      bitgq=bitgq*lobitgq
      xmsq=bitqg+bitgq

      if(inclqa) then
!==== bb~
         bitqa=assemblejet(order,tauc,
     &     beama0(5),beamb0(-5),beama1(5,:),beamb1(-5,:),
     &     beama2(5,:),beamb2(-5,:),soft1qa,soft2qa,jet1g,jet2g,hardqa)
!==== b~b
         bitqa=bitqa+assemblejet(order,tauc,
     &     beama0(-5),beamb0(5),beama1(-5,:),beamb1(5,:),
     &     beama2(-5,:),beamb2(5,:),soft1qa,soft2qa,jet1g,jet2g,hardqa)

         bitqa=bitqa*lobitqa
         xmsq=xmsq+bitqa
      endif

c---- debug
!      write(*,*) 'bitqg=',lobitqg
!      write(*,*) 'bitgq=',lobitgq
!      write(*,*) 'bitqa=',lobitqa
!      write(*,*) 'xmsq=',xmsq
!      write(*,*) '~~~the end~~~'
!      stop
c---- debug

      return
      end



