      subroutine qqb_dirgam_vbis(p,msq,order)
      implicit none
      include 'types.f'
!---  CW and JC May 16,
!---- hard functions for ga + jet 
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'zcouple.f'
      include 'ewcharge.f'
      include 'sprods_com.f'
      include 'nflav.f'
      integer:: j,k,order
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),fac,Qsum,msqgg
      real(dp):: qqbg0,qqbg1,qqbg2,qqb,qqb_nf
      real(dp):: qgq0,qgq1,qgq2,qg,qg_nf
      real(dp):: gqq0,gqq1,gqq2,gq,gq_nf
      real(dp):: qbq,gqb,qbg,qbq_nf,gqb_nf,qbg_nf
      real(dp):: Fqqb(3),Fgq(3),Fqg(3)
      integer, parameter:: inite1x1=1,inite0x2=2,initeNf=3
      
      Qsum=zip
      do j=1,nf
         Qsum=Qsum+Q(j)
      enddo
      
      msq(:,:)=zip
      msqgg=zip
      Fqqb(:)=zip
      Fgq(:)=zip
      Fqg(:)=zip
      
      fac=gsq*esq

! MadLoop phase space point for checking gg -> a g
!      p(1,4)=-5.0000000000000000e+02
!      p(1,1)=0.0000000000000000e+00
!      p(1,2)=0.0000000000000000e+00
!      p(1,3)=-5.0000000000000000e+02
!      p(2,4)=-5.0000000000000000e+02
!      p(2,1)=0.0000000000000000e+00 
!      p(2,2)=0.0000000000000000e+00
!      p(2,3)=+5.0000000000000000e+02
!      p(3,4)= 5.0000000000000000e+02
!      p(3,1)=1.1092428444383300e+02
!      p(3,2)=4.4483078948812101e+02
!      p(3,3)=-1.9955292993087900e+02
!      p(4,4)= 5.0000000000000000e+02
!      p(4,1)=-1.1092428444383300e+02
!      p(4,2)=-4.4483078948812101e+02
!      p(4,3)=1.9955292993087900e+02
      
c--- compute matrix elements

!===== qqb channel
      call gamamp_qqb_gag(order,p,1,2,3,4,qqbg0,qqbg1,qqbg2)
!===== q g channel 
      call gamamp_qg_gaq(order,p,1,2,3,4,qgq0,qgq1,qgq2)
!===== g q channel 
      call gamamp_qg_gaq(order,p,2,1,3,4,gqq0,gqq1,gqq2)

!=====call to finite functions if order ==2
      if(order==2) then
         call AGTY_dirgam_assemble_schan(p,1,2,3,Fqqb)
         call AGTY_dirgam_assemble_uchan(p,1,2,3,Fgq)
         call AGTY_dirgam_assemble_uchan(p,2,1,3,Fqg)
         call fill_gg_ga(p,1,2,4,3,msqgg)
      endif

      qqb_nf=zip
      gq_nf=zip
      qg_nf=zip
      if(order.eq.0) then
c--- apply overall factor
         qqb= aveqq*fac*(qqbg0)
         qg= -aveqg*fac*(qgq0)
         gq= -aveqg*fac*(gqq0)
      elseif(order.eq.1) then 
c--- apply overall factor and dress with ason2pi
         qqb= aveqq*fac*(ason2pi*qqbg1)
         qg= -aveqg*fac*(ason2pi*qgq1)
         gq= -aveqg*fac*(ason2pi*gqq1)
      elseif(order.eq.2) then
c--- apply overall factor and dress with ason2pi^2
         qqb= aveqq*fac*ason2pi**2*(qqbg2+Fqqb(inite0x2)+Fqqb(inite1x1))
         qg= -aveqg*fac*ason2pi**2*(qgq2+Fqg(inite0x2)+Fqg(inite1x1))
         gq= -aveqg*fac*ason2pi**2*(gqq2+Fgq(inite0x2)+Fgq(inite1x1))
         qqb_nf= aveqq*fac*ason2pi**2*Fqqb(initeNf)
         qg_nf= -aveqg*fac*ason2pi**2*Fqg(initeNf)
         gq_nf= -aveqg*fac*ason2pi**2*Fgq(initeNf)
      endif

! obtain remaining channels (qbq, gqb and qbg) by symmetry
      qbq=qqb
      gqb=gq
      qbg=qg
      qbq_nf=qqb_nf
      gqb_nf=gq_nf
      qbg_nf=qg_nf
      
      do j=-nflav,nflav
      do k=-nflav,nflav

         if( (j .ne. 0) .and. (k .ne. 0) .and. (j .ne. -k)) goto 20

         if((j.eq.0).and.(k.eq.0).and.(order==2)) then
            msq(j,k)=avegg*msqgg
         endif
         if((j.eq.-k).and.(j.ne.0).and.(j.gt.0)) then
            msq(j,k)=qqb*Q(j)**2
!-------additional nf-loop if order ==2
            if(order==2) then
               msq(j,k)=msq(j,k)+Q(j)*Qsum*qqb_nf
            endif
         elseif((j.eq.-k).and.(j.ne.0).and.(j.lt.0)) then
            msq(j,k)=qbq*Q(k)**2
!-------additional nf-loop if order ==2
            if(order==2) then
               msq(j,k)=msq(j,k)+Q(k)*Qsum*qbq_nf
            endif
         elseif((j.gt.0).and.(k.eq.0)) then
            msq(j,k)=qg*Q(j)**2
!-------additional nf-loop if order ==2
            if(order==2) then
               msq(j,k)=msq(j,k)+Q(j)*Qsum*qg_nf
            endif
         elseif((j.lt.0).and.(k.eq.0)) then
            msq(j,k)=qbg*Q(-j)**2
!-------additional nf-loop if order ==2
            if(order==2) then
               msq(j,k)=msq(j,k)+Q(-j)*Qsum*qbg_nf
            endif
         elseif((j.eq.0).and.(k.gt.0)) then
            msq(j,k)=gq*Q(k)**2
!-------additional nf-loop if order ==2
            if(order==2) then
               msq(j,k)=msq(j,k)+Q(k)*Qsum*gq_nf
            endif
         elseif((j.eq.0).and.(k.lt.0)) then
            msq(j,k)=gqb*Q(-k)**2
!-------additional nf-loop if order ==2
            if(order==2) then
               msq(j,k)=msq(j,k)+Q(-k)*Qsum*gqb_nf
            endif
         endif
            
 20      continue
      enddo
      enddo
      
!      if (order == 2) then
!        write(6,*) 'msq(0,0)',msq(0,0)
!        stop
!      endif
      
      return
      end
 

