      subroutine qqb_gam2j_gs(p,msq)
      implicit none
      include 'types.f'
************************************************************************
*     Author: J.M. Campbell                                            *
*     June, 2015.                                                      *
************************************************************************
c---Matrix element SUBTRACTION squared averaged over initial colors and spins
c     q(-p1)+qbar(-p2) -->  gamma(p3) + f(p4) + f(p4) + g(p5)
c
c     where the partons are either q(p4) and qbar(p4) [Qflag = .true.]
c                               or g4p4) and g(p4)    [Gflag = .true.]
c
c---This routine does not function if both Qflag and Gflag=.true.
c
c--- all momenta are outgoing
       
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'ptilde.f'
      include 'qqgg.f'
      include 'flags.f'
      include 'lc.f'
      include 'pp.f'
      include 'msq_cs.f'

c--- np6,np12,... = n+6,n+12,...
c -- nd counts the dipoles
      integer:: j,k,n,np6,np12,np18,np21
c--- slightly obtuse notation, fn=-nf, to simplify declaration lines
      real(dp):: p(mxpart,4),msq(maxd,fn:nf,fn:nf)
      real(dp):: 
     & msq16_2(fn:nf,fn:nf),msq26_1(fn:nf,fn:nf),
     & msq14_2(fn:nf,fn:nf),msq24_1(fn:nf,fn:nf),
     & msq15_2(fn:nf,fn:nf),msq25_1(fn:nf,fn:nf),
     & msq16_4(fn:nf,fn:nf),msq16_5(fn:nf,fn:nf),
     & msq14_5(fn:nf,fn:nf),msq54_6(fn:nf,fn:nf),
     & msq26_4(fn:nf,fn:nf),msq26_5(fn:nf,fn:nf),
     & msq46_5(fn:nf,fn:nf),msq56_4(fn:nf,fn:nf),
     & msq46_5v(fn:nf,fn:nf),msq25_1v(fn:nf,fn:nf),
     & msq14_2v(fn:nf,fn:nf),msq15_2v(fn:nf,fn:nf),
     & msq16_2v(fn:nf,fn:nf),msq24_5(fn:nf,fn:nf),
     & dummy(fn:nf,fn:nf),dummyv(fn:nf,fn:nf),
     & sub16_2(4),sub26_1(4),sub14_2(4),sub24_1(4),
     & sub15_2(4),sub25_1(4),sub24_5(4),
     & sub16_4(4),sub46_1(4),sub26_4(4),sub46_2(4),
     & sub16_5(4),sub56_1(4),sub26_5(4),sub56_2(4),
     & sub14_5(4),sub54_1(4),sub54_2(4),sub54_6(4),
     & sub46_5(4),sub56_4(4),dsubv,dsub(4),
     & sub46_5v,sub46_1v,sub46_2v,sub56_1v,sub56_2v,sub16_2v,sub54_6v,
     & sub25_1v,sub24_1v,sub26_1v,sub14_2v,sub15_2v,sub54_1v,sub54_2v,
     & sub56_4v
      real(dp)::
     & m16_2(0:2,fn:nf,fn:nf),m26_1(0:2,fn:nf,fn:nf),
     & m15_2(0:2,fn:nf,fn:nf),m25_1(0:2,fn:nf,fn:nf),
     & m14_2(0:2,fn:nf,fn:nf),m24_1(0:2,fn:nf,fn:nf),
     & m46_5(0:2,fn:nf,fn:nf),m56_4(0:2,fn:nf,fn:nf),
     & m54_6(0:2,fn:nf,fn:nf),
     & m16_4(0:2,fn:nf,fn:nf),
     & m16_5(0:2,fn:nf,fn:nf),
     & m26_4(0:2,fn:nf,fn:nf),
     & m26_5(0:2,fn:nf,fn:nf),
     & m14_5(0:2,fn:nf,fn:nf),
     & m24_5(0:2,fn:nf,fn:nf),
     & m16_2x(0:2,ppmax),
     & m26_1x(0:2,ppmax),
     & m14_2x(0:2,ppmax),
     & m15_2x(0:2,ppmax),
     & m24_1x(0:2,ppmax),
     & m25_1x(0:2,ppmax),
     & m46_5x(0:2,ppmax),
     & m54_6x(0:2,ppmax),
     & m56_4x(0:2,ppmax),
     & m16_4x(0:2,ppmax),
     & m16_5x(0:2,ppmax),
     & m26_4x(0:2,ppmax),
     & m26_5x(0:2,ppmax),
     & m14_5x(0:2,ppmax),
     & m24_5x(0:2,ppmax)
           
      real(dp):: 
     & msq1a_b(6,0:2,fn:nf,fn:nf),msqba_1(6,0:2,fn:nf,fn:nf),
     & msqab_c(6,0:2,fn:nf,fn:nf),
     & msqbc_2(6,0:2,fn:nf,fn:nf),msq2c_b(6,0:2,fn:nf,fn:nf),
     & sub1a_b(6,4),subba_1(6,4),subab_c(6,4),
     & subbc_2(6,4),sub2c_b(6,4),
     & msq1a_bv(6,0:2,fn:nf,fn:nf),msqba_1v(6,0:2,fn:nf,fn:nf),
     & msqab_cv(6,0:2,fn:nf,fn:nf),
     & msqbc_2v(6,0:2,fn:nf,fn:nf),msq2c_bv(6,0:2,fn:nf,fn:nf),
     & sub1a_bv(6),subba_1v(6),subab_cv(6),
     & subbc_2v(6),sub2c_bv(6),
     & msq1b_2(6,0:2,fn:nf,fn:nf),msq2b_1(6,0:2,fn:nf,fn:nf),
     & sub1b_2(6,4),sub2b_1(6,4),
     & msq1b_2v(6,0:2,fn:nf,fn:nf),msq2b_1v(6,0:2,fn:nf,fn:nf),
     & sub1b_2v(6),sub2b_1v(6)
      real(dp):: 
     & m46_1g(0:2,fn:nf,fn:nf),m46_1vg(0:2,fn:nf,fn:nf),
     & m46_1vx(ppmax),
     & m25_1g(0:2,fn:nf,fn:nf),m25_1vg(0:2,fn:nf,fn:nf),
     & m25_1vx(ppmax),
     & m24_1g(0:2,fn:nf,fn:nf),m24_1vg(0:2,fn:nf,fn:nf),
     & m24_1vx(ppmax),
     & m26_1g(0:2,fn:nf,fn:nf),m26_1vg(0:2,fn:nf,fn:nf),
     & m26_1vx(ppmax),
     & m56_1g(0:2,fn:nf,fn:nf),m56_1vg(0:2,fn:nf,fn:nf),
     & m56_1vx(ppmax),
     & m56_2g(0:2,fn:nf,fn:nf),m56_2vg(0:2,fn:nf,fn:nf),
     & m56_2vx(ppmax),
     & m46_2g(0:2,fn:nf,fn:nf),m46_2vg(0:2,fn:nf,fn:nf),
     & m46_2vx(ppmax),
     & m14_2g(0:2,fn:nf,fn:nf),m14_2vg(0:2,fn:nf,fn:nf),
     & m14_2vx(ppmax),
     & m15_2g(0:2,fn:nf,fn:nf),m15_2vg(0:2,fn:nf,fn:nf),
     & m16_2g(0:2,fn:nf,fn:nf),m16_2vg(0:2,fn:nf,fn:nf),
     & m15_2vx(ppmax),
     & m16_2vx(ppmax),
     & m54_1g(0:2,fn:nf,fn:nf),m54_1vg(0:2,fn:nf,fn:nf),
     & m54_1vx(ppmax),
     & m54_2g(0:2,fn:nf,fn:nf),m54_2vg(0:2,fn:nf,fn:nf),
     & m54_2vx(ppmax),
     & m14_5g(0:2,fn:nf,fn:nf),m24_5g(0:2,fn:nf,fn:nf),
     & m56_4vx(ppmax),m56_4vg(0:2,fn:nf,fn:nf),
     & m56_4g(0:2,fn:nf,fn:nf),
     & m54_6vx(ppmax),m54_6vg(0:2,fn:nf,fn:nf),
     & m54_6g(0:2,fn:nf,fn:nf),
     & m46_5vx(ppmax),m46_5vg(0:2,fn:nf,fn:nf),
     & m46_5g(0:2,fn:nf,fn:nf)
      
      real(dp):: mqq(0:2,fn:nf,fn:nf),
     & msqx(0:2,ppmax),mg(0:2,-nf:nf,-nf:nf),
     & mvg(0:2,-nf:nf,-nf:nf),mvxg(ppmax)
    
      external qqb_gam2j,qqb_gam2j_gvec
      external qqb_gam2jx_new,qqb_gam2j_gvecx_new,donothing_gvecx_new

      integer,parameter:: a(6)=(/4,4,6,5,5,6/),b(6)=(/5,6,4,6,4,5/),
     & c(6)=(/6,5,5,4,6,4/),
     & jj(-nf:nf)=(/-1,-2,-1,-2,-1,0,1,2,1,2,1/),
     & kk(-nf:nf)=(/-1,-2,-1,-2,-1,0,1,2,1,2,1/),
     & pntr(4:6,4:6)=reshape((/0,2,2,1,0,2,1,1,0/),(/3,3/))

      if (Qflag .and. Gflag) then
        write(6,*) 'Both Qflag and Gflag cannot be true'
        write(6,*) 'They are set in file options.DAT'
        write(6,*) 'Failed in qqb_gam2j_gs.f'
        stop
      endif
      
      ndmax=24

c-- initialize the matrix elements to zero
      msq(:,:,:)=zip

      if (Gflag) then

c---arguments of dips:
c---    1  dipole number
c---    2  momentum
c---    3  emitter
c---    4  emitted
c---    4  spectator
c---    5  interference-free subtraction, equiv. to AP kernel for qq,qg
c---    6  correlation piece of subtraction, relevant only for gg,gq
c---    8  lowest order matrix elements at rescaled momentum, msq(j,k)
c---    9  lowest order matrix elements at rescaled momentum
c---        with emitter contracted with appropriate vector, msqv(j,k)
c---   10  appropriate lowest order calculating routine for msq
c---   11  appropriate lowest order calculating routine for msqv

c--- final-final
      do n=1,6
      call dips(n,p,a(n),b(n),c(n),dsub,dsubv,dummy,dummyv,
     & qqb_gam2j,qqb_gam2j_gvec)
      call storedip(msqab_c,msqab_cv,dsub,dsubv,subab_c,subab_cv,n)
      enddo

c--- initial-final/final-initial
      do n=1,6
      np6=n+6
      np12=n+12
      call dips(np6,p,1,a(n),b(n),dsub,dsubv,dummy,dummyv,
     & qqb_gam2j,qqb_gam2j_gvec)
      call storedip(msq1a_b,msq1a_bv,dsub,dsubv,sub1a_b,sub1a_bv,n)
      call dips(np6,p,b(n),a(n),1,dsub,dsubv,dummy,dummyv,
     & qqb_gam2j,qqb_gam2j_gvec)
      call storedip(msqba_1,msqba_1v,dsub,dsubv,subba_1,subba_1v,n)
      call dips(np12,p,2,c(n),b(n),dsub,dsubv,dummy,dummyv,
     & qqb_gam2j,qqb_gam2j_gvec)
      call storedip(msq2c_b,msq2c_bv,dsub,dsubv,sub2c_b,sub2c_bv,n)
      call dips(np12,p,b(n),c(n),2,dsub,dsubv,dummy,dummyv,
     & qqb_gam2j,qqb_gam2j_gvec)
      call storedip(msqbc_2,msqbc_2v,dsub,dsubv,subbc_2,subbc_2v,n)
      enddo

c--- initial-initial
      do n=1,3
      np18=n+18
      np21=n+21
      call dips(np18,p,1,b(n),2,dsub,dsubv,dummy,dummyv,
     & qqb_gam2j,qqb_gam2j_gvec)
      call storedip(msq1b_2,msq1b_2v,dsub,dsubv,sub1b_2,sub1b_2v,n)
      call dips(np21,p,2,b(n),1,dsub,dsubv,dummy,dummyv,
     & qqb_gam2j,qqb_gam2j_gvec)
      call storedip(msq2b_1,msq2b_1v,dsub,dsubv,sub2b_1,sub2b_1v,n)
      enddo

c-- fill the matrix elements
      do j=-nf,nf
      do k=-nf,nf
      
c--- QUARK-ANTIQUARK contributions
      if    ((j>0).and.(k<0)) then
c------ leading colour
          if ((colourchoice == 1) .or. (colourchoice == 0)) then
          do n=1,6   
          msq(n,j,k)   =subab_c(n,gg)/2._dp
     &                  *(msqab_c(n,pntr(a(n),c(n)),j,k)
     &                   +msqab_c(n,pntr(c(n),a(n)),j,k))*xn/3._dp
     &                 +subab_cv(n)/2._dp
     &                  *(msqab_cv(n,pntr(a(n),c(n)),j,k)
     &                   +msqab_cv(n,pntr(c(n),a(n)),j,k))*xn/3._dp
          msq(6+n,j,k) =(sub1a_b(n,qq)+subba_1(n,gg)/2._dp)
     &                  *msq1a_b(n,pntr(b(n),c(n)),j,k)*xn/3._dp
     &                 +subba_1v(n)/2._dp
     &                  *msqba_1v(n,pntr(b(n),c(n)),j,k)*xn/3._dp
          msq(12+n,j,k)=(subbc_2(n,gg)/2._dp+sub2c_b(n,qq))
     &                  *msq2c_b(n,pntr(a(n),b(n)),j,k)*xn/3._dp
     &                 +subbc_2v(n)/2._dp
     &                  *msqbc_2v(n,pntr(a(n),b(n)),j,k)*xn/3._dp
          enddo
          endif
c------ sub-leading colour
          if ((colourchoice == 2) .or. (colourchoice == 0)) then
          do n=1,6 
          msq(6+n,j,k)=msq(6+n,j,k)
     &                 +(sub1a_b(n,qq)+subba_1(n,gg)/2._dp)
     &                  *msq1a_b(n,0,j,k)*xn/3._dp
     &                 +subba_1v(n)/2._dp
     &                  *msqba_1v(n,0,j,k)*xn/3._dp  
          msq(12+n,j,k)=msq(12+n,j,k)
     &                 +(subbc_2(n,gg)/2._dp+sub2c_b(n,qq))
     &                  *msq2c_b(n,0,j,k)*xn/3._dp
     &                 +subbc_2v(n)/2._dp
     &                  *msqbc_2v(n,0,j,k)*xn/3._dp
          enddo
          do n=1,3
          msq(18+n,j,k)=msq(18+n,j,k)
     &                   -sub1b_2(n,qq)
     &                  *(msq1b_2(n,1,j,k)+msq1b_2(n,2,j,k))/xn/3._dp
          msq(21+n,j,k)=msq(21+n,j,k)
     &                   -sub2b_1(n,qq)
     &                  *(msq2b_1(n,1,j,k)+msq2b_1(n,2,j,k))/xn/3._dp
         enddo
         endif
c------ sub-sub-leading colour
         if ((colourchoice == 3) .or. (colourchoice == 0)) then
          do n=1,3
          msq(18+n,j,k)=msq(18+n,j,k)
     &                  -sub1b_2(n,qq)
     &                  *msq1b_2(n,0,j,k)*(xn+1._dp/xn)/3._dp
          msq(21+n,j,k)=msq(21+n,j,k)
     &                  -sub2b_1(n,qq)
     &                  *msq2b_1(n,0,j,k)*(xn+1._dp/xn)/3._dp
         enddo         
         endif
c--- ANTIQUARK-QUARK contributions
      elseif((j<0).and.(k>0)) then
c------ leading colour
          if ((colourchoice == 1) .or. (colourchoice == 0)) then
          do n=1,6   
          msq(n,j,k)   =subab_c(n,gg)/2._dp
     &                  *(msqab_c(n,pntr(a(n),c(n)),j,k)
     &                   +msqab_c(n,pntr(c(n),a(n)),j,k))*xn/3._dp
     &                 +subab_cv(n)/2._dp
     &                  *(msqab_cv(n,pntr(a(n),c(n)),j,k)
     &                   +msqab_cv(n,pntr(c(n),a(n)),j,k))*xn/3._dp
          msq(6+n,j,k) =(sub1a_b(n,qq)+subba_1(n,gg)/2._dp)
     &                  *msq1a_b(n,pntr(c(n),b(n)),j,k)*xn/3._dp
     &                 +subba_1v(n)/2._dp
     &                  *msqba_1v(n,pntr(c(n),b(n)),j,k)*xn/3._dp
          msq(12+n,j,k)=(subbc_2(n,gg)/2._dp+sub2c_b(n,qq))
     &                  *msq2c_b(n,pntr(b(n),a(n)),j,k)*xn/3._dp
     &                 +subbc_2v(n)/2._dp
     &                  *msqbc_2v(n,pntr(b(n),a(n)),j,k)*xn/3._dp
          enddo
          endif
c------ sub-leading colour
          if ((colourchoice == 2) .or. (colourchoice == 0)) then
          do n=1,6 
          msq(6+n,j,k)=msq(6+n,j,k)
     &                 +(sub1a_b(n,qq)+subba_1(n,gg)/2._dp)
     &                  *msq1a_b(n,0,j,k)*xn/3._dp
     &                 +subba_1v(n)/2._dp
     &                  *msqba_1v(n,0,j,k)*xn/3._dp  
          msq(12+n,j,k)=msq(12+n,j,k)
     &                 +(subbc_2(n,gg)/2._dp+sub2c_b(n,qq))
     &                  *msq2c_b(n,0,j,k)*xn/3._dp
     &                 +subbc_2v(n)/2._dp
     &                  *msqbc_2v(n,0,j,k)*xn/3._dp
          enddo
          do n=1,3
          msq(18+n,j,k)=msq(18+n,j,k)
     &                   -sub1b_2(n,qq)
     &                  *(msq1b_2(n,1,j,k)+msq1b_2(n,2,j,k))/xn/3._dp
          msq(21+n,j,k)=msq(21+n,j,k)
     &                   -sub2b_1(n,qq)
     &                  *(msq2b_1(n,1,j,k)+msq2b_1(n,2,j,k))/xn/3._dp
         enddo
         endif
c------ sub-sub-leading colour
         if ((colourchoice == 3) .or. (colourchoice == 0)) then
          do n=1,3
          msq(18+n,j,k)=msq(18+n,j,k)
     &                  -sub1b_2(n,qq)
     &                  *msq1b_2(n,0,j,k)*(xn+1._dp/xn)/3._dp
          msq(21+n,j,k)=msq(21+n,j,k)
     &                  -sub2b_1(n,qq)
     &                  *msq2b_1(n,0,j,k)*(xn+1._dp/xn)/3._dp
         enddo         
         endif
c--- GLUON-GLUON contributions
      elseif ((j==0) .and. (k==0)) then
c------ leading colour
          if ((colourchoice == 1) .or. (colourchoice == 0)) then
c--- choose n=2 which is (4,6,5)
          do n=2,2
          msq(18+n,j,k)=sub1b_2(n,gg)
     &                  *(msq1b_2(n,1,j,k)+msq1b_2(n,2,j,k))*xn
     &                 +sub1b_2v(n)
     &                  *(msq1b_2v(n,1,j,k)+msq1b_2v(n,2,j,k))*xn
          msq(21+n,j,k)=sub2b_1(n,gg)
     &                  *(msq2b_1(n,1,j,k)+msq2b_1(n,2,j,k))*xn
     &                 +sub2b_1v(n)
     &                  *(msq2b_1v(n,1,j,k)+msq2b_1v(n,2,j,k))*xn
          enddo
c--- choose n=3,5 which is (6,4,5) and (6,5,4)
          do n=3,6,3
          msq(6+n,j,k)=(sub1a_b(n,gg)+subba_1(n,qq))
     &                  *msq1a_b(n,pntr(b(n),c(n)),j,k)*xn
     &                 +sub1a_bv(n)
     &                  *msq1a_bv(n,pntr(b(n),c(n)),j,k)*xn
          enddo
c--- choose n=1,4 which is (4,5,6) and (5,4,6)
          do n=1,5,4
          msq(12+n,j,k)=(sub2c_b(n,gg)+subbc_2(n,qq))
     &                   *msq2c_b(n,pntr(a(n),b(n)),j,k)*xn
     &                  +sub2c_bv(n)
     &                   *msq2c_bv(n,pntr(a(n),b(n)),j,k)*xn
         enddo
c--- additional (qg) collinear contributions
c--- choose n=3 which is (6,4,5)
          do n=3,3
          msq(18+n,j,k)=msq(18+n,j,k)+sub1b_2(n,qg)*
     &      (msq1b_2(n,1,+1,k)+msq1b_2(n,1,+2,k)+msq1b_2(n,1,+3,k)
     &      +msq1b_2(n,1,+4,k)+msq1b_2(n,1,+5,k)
     &      +msq1b_2(n,2,+1,k)+msq1b_2(n,2,+2,k)+msq1b_2(n,2,+3,k)
     &      +msq1b_2(n,2,+4,k)+msq1b_2(n,2,+5,k))*xn*(avegg/aveqg)
          msq(21+n,j,k)=msq(21+n,j,k)+sub2b_1(n,qg)*
     &      (msq2b_1(n,1,j,+1)+msq2b_1(n,1,j,+2)+msq2b_1(n,1,j,+3)
     &      +msq2b_1(n,1,j,+4)+msq2b_1(n,1,j,+5)
     &      +msq2b_1(n,2,j,+1)+msq2b_1(n,2,j,+2)+msq2b_1(n,2,j,+3)
     &      +msq2b_1(n,2,j,+4)+msq2b_1(n,2,j,+5))*xn*(avegg/aveqg)
          enddo
c--- choose n=1 which is (4,5,6)
          do n=1,1
          msq(18+n,j,k)=msq(18+n,j,k)+sub1b_2(n,qg)*
     &      (msq1b_2(n,1,-1,k)+msq1b_2(n,1,-2,k)+msq1b_2(n,1,-3,k)
     &      +msq1b_2(n,1,-4,k)+msq1b_2(n,1,-5,k)
     &      +msq1b_2(n,2,-1,k)+msq1b_2(n,2,-2,k)+msq1b_2(n,2,-3,k)
     &      +msq1b_2(n,2,-4,k)+msq1b_2(n,2,-5,k))*xn*(avegg/aveqg)
          msq(21+n,j,k)=msq(21+n,j,k)+sub2b_1(n,qg)*
     &      (msq2b_1(n,1,j,-1)+msq2b_1(n,1,j,-2)+msq2b_1(n,1,j,-3)
     &      +msq2b_1(n,1,j,-4)+msq2b_1(n,1,j,-5)
     &      +msq2b_1(n,2,j,-1)+msq2b_1(n,2,j,-2)+msq2b_1(n,2,j,-3)
     &      +msq2b_1(n,2,j,-4)+msq2b_1(n,2,j,-5))*xn*(avegg/aveqg)
          enddo
          endif
c------ sub-leading colour
          if ((colourchoice == 2) .or. (colourchoice == 0)) then
          do n=2,4,2
c-- (a,b,c) = (4,6,5) and (5,6,4)
          msq(n,j,k)   =msq(n,j,k)
     &                   -subab_c(n,qq)
     &                  *(msqab_c(n,1,j,k)+msqab_c(n,2,j,k))/xn
          enddo
          do n=3,6,3
c-- (a,b,c) = (6,4,5) and (6,5,4)
          msq(6+n,j,k)=msq(6+n,j,k)
     &                 +(sub1a_b(n,gg)+subba_1(n,qq))
     &                  *msq1a_b(n,0,j,k)*xn
     &                 +sub1a_bv(n)
     &                  *msq1a_bv(n,0,j,k)*xn
          enddo
c--- choose n=1,4 which is (4,5,6) and (5,4,6)
          do n=1,5,4
          msq(12+n,j,k)=msq(12+n,j,k)
     &                  +(sub2c_b(n,gg)+subbc_2(n,qq))
     &                   *msq2c_b(n,0,j,k)*xn
     &                  +sub2c_bv(n)
     &                   *msq2c_bv(n,0,j,k)*xn
         enddo
c--- additional (qg) collinear contributions
c--- choose n=3 which is (6,4,5)
          do n=3,3
          msq(18+n,j,k)=msq(18+n,j,k)+sub1b_2(n,qg)*(avegg/aveqg)*
     &    (-(msq1b_2(n,1,+1,k)+msq1b_2(n,1,+2,k)+msq1b_2(n,1,+3,k)
     &      +msq1b_2(n,1,+4,k)+msq1b_2(n,1,+5,k)
     &      +msq1b_2(n,2,+1,k)+msq1b_2(n,2,+2,k)+msq1b_2(n,2,+3,k)
     &      +msq1b_2(n,2,+4,k)+msq1b_2(n,2,+5,k))/xn
     &     +(msq1b_2(n,0,+1,k)+msq1b_2(n,0,+2,k)+msq1b_2(n,0,+3,k)
     &      +msq1b_2(n,0,+4,k)+msq1b_2(n,0,+5,k))*2._dp*xn)
          msq(21+n,j,k)=msq(21+n,j,k)+sub2b_1(n,qg)*(avegg/aveqg)*
     &    (-(msq2b_1(n,1,j,+1)+msq2b_1(n,1,j,+2)+msq2b_1(n,1,j,+3)
     &      +msq2b_1(n,1,j,+4)+msq2b_1(n,1,j,+5)
     &      +msq2b_1(n,2,j,+1)+msq2b_1(n,2,j,+2)+msq2b_1(n,2,j,+3)
     &      +msq2b_1(n,2,j,+4)+msq2b_1(n,2,j,+5))/xn
     &     +(msq2b_1(n,0,j,+1)+msq2b_1(n,0,j,+2)+msq2b_1(n,0,j,+3)
     &      +msq2b_1(n,0,j,+4)+msq2b_1(n,0,j,+5))*2._dp*xn)
          enddo
c--- choose n=1 which is (4,5,6)
          do n=1,1
          msq(18+n,j,k)=msq(18+n,j,k)+sub1b_2(n,qg)*(avegg/aveqg)*
     &    (-(msq1b_2(n,1,-1,k)+msq1b_2(n,1,-2,k)+msq1b_2(n,1,-3,k)
     &      +msq1b_2(n,1,-4,k)+msq1b_2(n,1,-5,k)
     &      +msq1b_2(n,2,-1,k)+msq1b_2(n,2,-2,k)+msq1b_2(n,2,-3,k)
     &      +msq1b_2(n,2,-4,k)+msq1b_2(n,2,-5,k))/xn
     &     +(msq1b_2(n,0,-1,k)+msq1b_2(n,0,-2,k)+msq1b_2(n,0,-3,k)
     &      +msq1b_2(n,0,-4,k)+msq1b_2(n,0,-5,k))*2._dp*xn)
          msq(21+n,j,k)=msq(21+n,j,k)+sub2b_1(n,qg)*(avegg/aveqg)*
     &    (-(msq2b_1(n,1,j,-1)+msq2b_1(n,1,j,-2)+msq2b_1(n,1,j,-3)
     &      +msq2b_1(n,1,j,-4)+msq2b_1(n,1,j,-5)
     &      +msq2b_1(n,2,j,-1)+msq2b_1(n,2,j,-2)+msq2b_1(n,2,j,-3)
     &      +msq2b_1(n,2,j,-4)+msq2b_1(n,2,j,-5))/xn
     &     +(msq2b_1(n,0,j,-1)+msq2b_1(n,0,j,-2)+msq2b_1(n,0,j,-3)
     &      +msq2b_1(n,0,j,-4)+msq2b_1(n,0,j,-5))*2._dp*xn)
          enddo
          endif
c------ sub-sub-leading colour
          if ((colourchoice == 3) .or. (colourchoice == 0)) then
          do n=2,4,2
c-- (a,b,c) = (4,6,5) and (5,6,4)
          msq(n,j,k)   =msq(n,j,k)
     &                  -subab_c(n,qq)
     &                  *msqab_c(n,0,j,k)*(xn+1._dp/xn)
          enddo
c--- additional (qg) collinear contributions
c--- choose n=3 which is (6,4,5)
          do n=3,3
          msq(18+n,j,k)=msq(18+n,j,k)+sub1b_2(n,qg)*(avegg/aveqg)*
     &      (msq1b_2(n,0,+1,k)+msq1b_2(n,0,+2,k)+msq1b_2(n,0,+3,k)
     &      +msq1b_2(n,0,+4,k)+msq1b_2(n,0,+5,k))*(-xn-1._dp/xn)
          msq(21+n,j,k)=msq(21+n,j,k)+sub2b_1(n,qg)*(avegg/aveqg)*
     &      (msq2b_1(n,0,j,+1)+msq2b_1(n,0,j,+2)+msq2b_1(n,0,j,+3)
     &      +msq2b_1(n,0,j,+4)+msq2b_1(n,0,j,+5))*(-xn-1._dp/xn)
          enddo
c--- choose n=1 which is (4,5,6)
          do n=1,1
          msq(18+n,j,k)=msq(18+n,j,k)+sub1b_2(n,qg)*(avegg/aveqg)*
     &      (msq1b_2(n,0,-1,k)+msq1b_2(n,0,-2,k)+msq1b_2(n,0,-3,k)
     &      +msq1b_2(n,0,-4,k)+msq1b_2(n,0,-5,k))*(-xn-1._dp/xn)
          msq(21+n,j,k)=msq(21+n,j,k)+sub2b_1(n,qg)*(avegg/aveqg)*
     &      (msq2b_1(n,0,j,-1)+msq2b_1(n,0,j,-2)+msq2b_1(n,0,j,-3)
     &      +msq2b_1(n,0,j,-4)+msq2b_1(n,0,j,-5))*(-xn-1._dp/xn)
          enddo
          endif
c--- QUARK-GLUON contributions
      elseif ((j>0) .and. (k==0)) then
c------ leading colour
          if ((colourchoice == 1) .or. (colourchoice == 0)) then
          do n=1,2  
c-- (a,b,c) = (4,5,6) and (4,6,5)
          msq(n,j,k)   =subab_c(n,qq)
     &                  *msqab_c(n,pntr(a(n),c(n)),j,k)*xn/2._dp
          msq(12+n,j,k)=(subbc_2(n,gg)/2._dp+sub2c_b(n,gg))
     &                  *(msq2c_b(n,pntr(a(n),b(n)),j,k)
     &                   +msq2c_b(n,pntr(b(n),a(n)),j,k))*xn/2._dp
     &                 +subbc_2v(n)/2._dp
     &                  *(msqbc_2v(n,pntr(a(n),b(n)),j,k)
     &                   +msqbc_2v(n,pntr(b(n),a(n)),j,k))*xn/2._dp
     &                 +sub2c_bv(n)
     &                  *(msq2c_bv(n,pntr(a(n),b(n)),j,k)
     &                   +msq2c_bv(n,pntr(b(n),a(n)),j,k))*xn/2._dp
          msq(18+n,j,k)=sub1b_2(n,qq)
     &                  *msq1b_2(n,pntr(a(n),c(n)),j,k)*xn/2._dp
          msq(21+n,j,k)=sub2b_1(n,gg)
     &                  *msq2b_1(n,pntr(a(n),c(n)),j,k)*xn/2._dp
     &                 +sub2b_1v(n)
     &                  *msq2b_1v(n,pntr(a(n),c(n)),j,k)*xn/2._dp
          enddo
          do n=4,6,2
c-- (a,b,c) = (5,6,4) and (6,5,4)
          msq(n,j,k)   =subab_c(n,gg)/2._dp
     &                  *msqab_c(n,pntr(c(n),a(n)),j,k)*xn/2._dp
     &                 +subab_cv(n)/2._dp
     &                  *msqab_cv(n,pntr(c(n),a(n)),j,k)*xn/2._dp
          msq(6+n,j,k) =(sub1a_b(n,qq)+subba_1(n,gg)/2._dp)
     &                  *msq1a_b(n,pntr(b(n),c(n)),j,k)*xn/2._dp
     &                 +subba_1v(n)/2._dp
     &                  *msqba_1v(n,pntr(b(n),c(n)),j,k)*xn/2._dp
          enddo
          do n=3,5,2
c-- (a,b,c) = (6,4,5) and (5,4,6)
          msq(12+n,j,k)=(subbc_2(n,qq)+sub2c_b(n,gg))
     &                  *msq2c_b(n,pntr(a(n),b(n)),j,k)*xn/2._dp
     &                 +sub2c_bv(n)
     &                  *msq2c_bv(n,pntr(a(n),b(n)),j,k)*xn/2._dp
          enddo      
c--- additional (qg) collinear contributions
c--- choose n=3 which is (6,4,5)
          do n=3,3
          msq(21+n,j,k)=msq(21+n,j,k)+sub2b_1(n,qg)*
     &      (msq2b_1(n,1,j,-1)+msq2b_1(n,1,j,-2)+msq2b_1(n,1,j,-3)
     &      +msq2b_1(n,1,j,-4)+msq2b_1(n,1,j,-5)
     &      +msq2b_1(n,2,j,-1)+msq2b_1(n,2,j,-2)+msq2b_1(n,2,j,-3)
     &      +msq2b_1(n,2,j,-4)+msq2b_1(n,2,j,-5))*xn*(aveqg/aveqq)
          enddo
          endif
c------ sub-leading colour
          if ((colourchoice == 2) .or. (colourchoice == 0)) then
          do n=4,6,2
c-- (a,b,c) = (5,6,4) and (6,5,4)
          msq(6+n,j,k) =msq(6+n,j,k)
     &                 +(sub1a_b(n,qq)+subba_1(n,gg)/2._dp)
     &                  *msq1a_b(n,0,j,k)*xn/2._dp
     &                 +subba_1v(n)/2._dp
     &                  *msqba_1v(n,0,j,k)*xn/2._dp
          enddo
          do n=3,5,2
c-- (a,b,c) = (6,4,5) and (5,4,6)
          msq(6+n,j,k)=msq(6+n,j,k)
     &                 -(subba_1(n,qq)+sub1a_b(n,qq))
     &                  *(msq1a_b(n,1,j,k)+msq1a_b(n,2,j,k))/xn/2._dp
          msq(12+n,j,k)=msq(12+n,j,k)
     &                 +(subbc_2(n,qq)+sub2c_b(n,gg))
     &                  *msq2c_b(n,0,j,k)*xn/2._dp
     &                 +sub2c_bv(n)
     &                  *msq2c_bv(n,0,j,k)*xn/2._dp
          enddo      
          do n=1,2  
c-- (a,b,c) = (4,5,6) and (4,6,5)
          msq(n,j,k)   =msq(n,j,k)
     &                 +subab_c(n,qq)
     &                  *msqab_c(n,0,j,k)*xn/2._dp
          msq(18+n,j,k)=msq(18+n,j,k)
     &                 +sub1b_2(n,qq)
     &                  *msq1b_2(n,0,j,k)*xn/2._dp
          msq(21+n,j,k)=msq(21+n,j,k)
     &                 +sub2b_1(n,gg)
     &                  *msq2b_1(n,0,j,k)*xn/2._dp
     &                 +sub2b_1v(n)
     &                  *msq2b_1v(n,0,j,k)*xn/2._dp
          enddo    
c-- [n=4] (a,b,c) = (5,6,4)
          msq(4,j,k)  =msq(4,j,k)
     &                 +subab_c(4,gg)
     &                  *msqab_c(4,0,j,k)*xn/2._dp
     &                 +subab_cv(4)
     &                  *msqab_cv(4,0,j,k)*xn/2._dp
c--- additional (qg) collinear contributions
c--- choose n=3 which is (6,4,5)
          do n=3,3
          msq(21+n,j,k)=msq(21+n,j,k)+sub2b_1(n,qg)*(aveqg/aveqq)*
     &    (-(msq2b_1(n,1,j,-1)+msq2b_1(n,1,j,-2)+msq2b_1(n,1,j,-3)
     &      +msq2b_1(n,1,j,-4)+msq2b_1(n,1,j,-5)
     &      +msq2b_1(n,2,j,-1)+msq2b_1(n,2,j,-2)+msq2b_1(n,2,j,-3)
     &      +msq2b_1(n,2,j,-4)+msq2b_1(n,2,j,-5))/xn
     &     +(msq2b_1(n,0,j,-1)+msq2b_1(n,0,j,-2)+msq2b_1(n,0,j,-3)
     &      +msq2b_1(n,0,j,-4)+msq2b_1(n,0,j,-5))*2._dp*xn)
          enddo          
          endif
c------ sub-sub-leading colour
          if ((colourchoice == 3) .or. (colourchoice == 0)) then
          do n=3,5,2
c-- (a,b,c) = (6,4,5) and (5,4,6)
          msq(6+n,j,k)=msq(6+n,j,k)
     &                 -(subba_1(n,qq)+sub1a_b(n,qq))
     &                  *msq1a_b(n,0,j,k)*(xn+1._dp/xn)/2._dp
          enddo
c--- additional (qg) collinear contributions
c--- choose n=3 which is (6,4,5)
          do n=3,3
          msq(21+n,j,k)=msq(21+n,j,k)+sub2b_1(n,qg)*(aveqg/aveqq)*
     &      (msq2b_1(n,0,j,-1)+msq2b_1(n,0,j,-2)+msq2b_1(n,0,j,-3)
     &      +msq2b_1(n,0,j,-4)+msq2b_1(n,0,j,-5))*(-xn-1._dp/xn)
          enddo          
          endif
c--- GLUON-QUARK contributions
      elseif ((j==0) .and. (k>0)) then
c------ leading colour
          if ((colourchoice == 1) .or. (colourchoice == 0)) then
          do n=1,2  
c-- (a,b,c) = (4,5,6) and (4,6,5)
          msq(n,j,k)   =subab_c(n,qq)
     &                  *msqab_c(n,pntr(a(n),c(n)),j,k)*xn/2._dp
          msq(12+n,j,k) =(sub2c_b(n,qq)+subbc_2(n,gg)/2._dp)
     &                    *msq2c_b(n,pntr(c(n),a(n)),j,k)*xn/2._dp
     &                   +subbc_2v(n)/2._dp
     &                    *msqbc_2v(n,pntr(c(n),a(n)),j,k)*xn/2._dp
          msq(18+n,j,k)=sub1b_2(n,gg)
     &                  *msq1b_2(n,pntr(a(n),c(n)),j,k)*xn/2._dp
     &                 +sub1b_2v(n)
     &                  *msq1b_2v(n,pntr(a(n),c(n)),j,k)*xn/2._dp
          msq(21+n,j,k)=sub2b_1(n,qq)
     &                  *msq2b_1(n,pntr(a(n),c(n)),j,k)*xn/2._dp
          enddo
          do n=4,6,2
c-- (a,b,c) = (5,6,4) and (6,5,4)
          msq(n,j,k)   =subab_c(n,gg)/2._dp
     &                  *msqab_c(n,pntr(c(n),a(n)),j,k)*xn/2._dp
     &                 +subab_cv(n)/2._dp
     &                  *msqab_cv(n,pntr(c(n),a(n)),j,k)*xn/2._dp
          msq(6+n,j,k)=(subba_1(n,gg)/2._dp+sub1a_b(n,gg))
     &                 *(msq1a_b(n,pntr(c(n),b(n)),j,k)
     &                  +msq1a_b(n,pntr(b(n),c(n)),j,k))*xn/2._dp
     &                +subba_1v(n)/2._dp
     &                 *(msqba_1v(n,pntr(c(n),b(n)),j,k)
     &                  +msqba_1v(n,pntr(b(n),c(n)),j,k))*xn/2._dp
     &                +sub1a_bv(n)
     &                 *(msq1a_bv(n,pntr(c(n),b(n)),j,k)
     &                  +msq1a_bv(n,pntr(b(n),c(n)),j,k))*xn/2._dp
          enddo
          do n=3,5,2
c-- (a,b,c) = (5,4,6) and (6,4,5)
          msq(6+n,j,k)=(subba_1(n,qq)+sub1a_b(n,gg))
     &                  *msq1a_b(n,pntr(c(n),b(n)),j,k)*xn/2._dp
     &                 +sub1a_bv(n)
     &                  *msq1a_bv(n,pntr(c(n),b(n)),j,k)*xn/2._dp
          enddo              
c--- additional (qg) collinear contributions
c--- choose n=3 which is (6,4,5)
          do n=3,3
          msq(18+n,j,k)=msq(18+n,j,k)+sub1b_2(n,qg)*
     &      (msq1b_2(n,1,-1,k)+msq1b_2(n,1,-2,k)+msq1b_2(n,1,-3,k)
     &      +msq1b_2(n,1,-4,k)+msq1b_2(n,1,-5,k)
     &      +msq1b_2(n,2,-1,k)+msq1b_2(n,2,-2,k)+msq1b_2(n,2,-3,k)
     &      +msq1b_2(n,2,-4,k)+msq1b_2(n,2,-5,k))*xn*(aveqg/aveqq)
          enddo
          endif
c------ sub-leading colour
          if ((colourchoice == 2) .or. (colourchoice == 0)) then
          do n=1,2
c-- (a,b,c) = (4,6,5) and (4,5,6)
          msq(12+n,j,k) =msq(12+n,j,k)
     &                  +(sub2c_b(n,qq)+subbc_2(n,gg)/2._dp)
     &                   *msq2c_b(n,0,j,k)*xn/2._dp
     &                  +subbc_2v(n)/2._dp
     &                   *msqbc_2v(n,0,j,k)*xn/2._dp
          enddo
          do n=3,5,2
c-- (a,b,c) = (6,4,5) and (5,4,6)
          msq(12+n,j,k)=msq(12+n,j,k)
     &                 -(subbc_2(n,qq)+sub2c_b(n,qq))
     &                  *(msq2c_b(n,1,j,k)+msq2c_b(n,2,j,k))/xn/2._dp
          msq(6+n,j,k) =msq(6+n,j,k)
     &                 +(subba_1(n,qq)+sub1a_b(n,gg))
     &                  *msq1a_b(n,0,j,k)*xn/2._dp
     &                 +sub1a_bv(n)
     &                  *msq1a_bv(n,0,j,k)*xn/2._dp
          enddo      
          do n=1,2  
c-- (a,b,c) = (4,5,6) and (4,6,5)
          msq(n,j,k)   =msq(n,j,k)
     &                 +subab_c(n,qq)
     &                  *msqab_c(n,0,j,k)*xn/2._dp
          msq(18+n,j,k)=msq(18+n,j,k)
     &                 +sub1b_2(n,gg)
     &                  *msq1b_2(n,0,j,k)*xn/2._dp
     &                 +sub1b_2v(n)
     &                  *msq1b_2v(n,0,j,k)*xn/2._dp
          msq(21+n,j,k)=msq(21+n,j,k)
     &                 +sub2b_1(n,qq)
     &                  *msq2b_1(n,0,j,k)*xn/2._dp
          enddo    
c-- [n=4] (a,b,c) = (5,6,4)
          msq(4,j,k)  =msq(4,j,k)
     &                 +subab_c(4,gg)
     &                  *msqab_c(4,0,j,k)*xn/2._dp                
     &                 +subab_cv(4)
     &                  *msqab_cv(4,0,j,k)*xn/2._dp                
c--- choose n=3 which is (6,4,5)
          do n=3,3
          msq(18+n,j,k)=msq(18+n,j,k)+sub1b_2(n,qg)*(aveqg/aveqq)*
     &    (-(msq1b_2(n,1,-1,k)+msq1b_2(n,1,-2,k)+msq1b_2(n,1,-3,k)
     &      +msq1b_2(n,1,-4,k)+msq1b_2(n,1,-5,k)
     &      +msq1b_2(n,2,-1,k)+msq1b_2(n,2,-2,k)+msq1b_2(n,2,-3,k)
     &      +msq1b_2(n,2,-4,k)+msq1b_2(n,2,-5,k))/xn
     &     +(msq1b_2(n,0,-1,k)+msq1b_2(n,0,-2,k)+msq1b_2(n,0,-3,k)
     &      +msq1b_2(n,0,-4,k)+msq1b_2(n,0,-5,k))*2._dp*xn)
          enddo
          endif
c------ sub-sub-leading colour
          if ((colourchoice == 3) .or. (colourchoice == 0)) then
          do n=3,5,2
c-- (a,b,c) = (6,4,5) and (5,4,6)
          msq(12+n,j,k)=msq(12+n,j,k)
     &                 -(subbc_2(n,qq)+sub2c_b(n,qq))
     &                  *msq2c_b(n,0,j,k)*(xn+1._dp/xn)/2._dp
          enddo
c--- choose n=3 which is (6,4,5)
          do n=3,3
          msq(18+n,j,k)=msq(18+n,j,k)+sub1b_2(n,qg)*(aveqg/aveqq)*
     &      (msq1b_2(n,0,-1,k)+msq1b_2(n,0,-2,k)+msq1b_2(n,0,-3,k)
     &      +msq1b_2(n,0,-4,k)+msq1b_2(n,0,-5,k))*(-xn-1._dp/xn)
          enddo
          endif
c--- GLUON-ANTIQUARK contributions
      elseif ((j==0) .and. (k<0)) then
c------ leading colour
          if ((colourchoice == 1) .or. (colourchoice == 0)) then
          do n=1,2  
c-- (a,b,c) = (4,5,6) and (4,6,5)
          msq(n,j,k)   =subab_c(n,qq)
     &                  *msqab_c(n,pntr(c(n),a(n)),j,k)*xn/2._dp
          msq(12+n,j,k) =(sub2c_b(n,qq)+subbc_2(n,gg)/2._dp)
     &                    *msq2c_b(n,pntr(a(n),c(n)),j,k)*xn/2._dp
     &                   +subbc_2v(n)/2._dp
     &                    *msqbc_2v(n,pntr(a(n),c(n)),j,k)*xn/2._dp
          msq(18+n,j,k)=sub1b_2(n,gg)
     &                  *msq1b_2(n,pntr(c(n),a(n)),j,k)*xn/2._dp
     &                 +sub1b_2v(n)
     &                  *msq1b_2v(n,pntr(c(n),a(n)),j,k)*xn/2._dp
          msq(21+n,j,k)=sub2b_1(n,qq)
     &                  *msq2b_1(n,pntr(c(n),a(n)),j,k)*xn/2._dp
          enddo
          do n=4,6,2
c-- (a,b,c) = (5,6,4) and (6,5,4)
          msq(n,j,k)   =subab_c(n,gg)/2._dp
     &                  *msqab_c(n,pntr(a(n),c(n)),j,k)*xn/2._dp
     &                 +subab_cv(n)/2._dp
     &                  *msqab_cv(n,pntr(a(n),c(n)),j,k)*xn/2._dp
          msq(6+n,j,k)=(subba_1(n,gg)/2._dp+sub1a_b(n,gg))
     &                 *(msq1a_b(n,pntr(b(n),c(n)),j,k)
     &                  +msq1a_b(n,pntr(c(n),b(n)),j,k))*xn/2._dp
     &                +subba_1v(n)/2._dp
     &                 *(msqba_1v(n,pntr(b(n),c(n)),j,k)
     &                  +msqba_1v(n,pntr(c(n),b(n)),j,k))*xn/2._dp
     &                +sub1a_bv(n)
     &                 *(msq1a_bv(n,pntr(b(n),c(n)),j,k)
     &                  +msq1a_bv(n,pntr(c(n),b(n)),j,k))*xn/2._dp
          enddo
          do n=3,5,2
c-- (a,b,c) = (5,4,6) and (6,4,5)
          msq(6+n,j,k)=(subba_1(n,qq)+sub1a_b(n,gg))
     &                  *msq1a_b(n,pntr(b(n),c(n)),j,k)*xn/2._dp
     &                 +sub1a_bv(n)
     &                  *msq1a_bv(n,pntr(b(n),c(n)),j,k)*xn/2._dp
          enddo              
c--- choose n=1 which is (6,4,5)
          do n=3,3
          msq(18+n,j,k)=msq(18+n,j,k)+sub1b_2(n,qg)*
     &      (msq1b_2(n,1,+1,k)+msq1b_2(n,1,+2,k)+msq1b_2(n,1,+3,k)
     &      +msq1b_2(n,1,+4,k)+msq1b_2(n,1,+5,k)
     &      +msq1b_2(n,2,+1,k)+msq1b_2(n,2,+2,k)+msq1b_2(n,2,+3,k)
     &      +msq1b_2(n,2,+4,k)+msq1b_2(n,2,+5,k))*xn*(aveqg/aveqq)
          enddo
          endif
c------ sub-leading colour
          if ((colourchoice == 2) .or. (colourchoice == 0)) then
          do n=1,2
c-- (a,b,c) = (4,6,5) and (4,5,6)
          msq(12+n,j,k) =msq(12+n,j,k)
     &                  +(sub2c_b(n,qq)+subbc_2(n,gg)/2._dp)
     &                   *msq2c_b(n,0,j,k)*xn/2._dp
     &                  +subbc_2v(n)/2._dp
     &                   *msqbc_2v(n,0,j,k)*xn/2._dp
          enddo
          do n=3,5,2
c-- (a,b,c) = (6,4,5) and (5,4,6)
          msq(12+n,j,k)=msq(12+n,j,k)
     &                 -(subbc_2(n,qq)+sub2c_b(n,qq))
     &                  *(msq2c_b(n,1,j,k)+msq2c_b(n,2,j,k))/xn/2._dp
          msq(6+n,j,k) =msq(6+n,j,k)
     &                 +(subba_1(n,qq)+sub1a_b(n,gg))
     &                  *msq1a_b(n,0,j,k)*xn/2._dp
     &                 +sub1a_bv(n)
     &                  *msq1a_bv(n,0,j,k)*xn/2._dp
          enddo      
          do n=1,2  
c-- (a,b,c) = (4,5,6) and (4,6,5)
          msq(n,j,k)   =msq(n,j,k)
     &                 +subab_c(n,qq)
     &                  *msqab_c(n,0,j,k)*xn/2._dp
          msq(18+n,j,k)=msq(18+n,j,k)
     &                 +sub1b_2(n,gg)
     &                  *msq1b_2(n,0,j,k)*xn/2._dp
     &                 +sub1b_2v(n)
     &                  *msq1b_2v(n,0,j,k)*xn/2._dp
          msq(21+n,j,k)=msq(21+n,j,k)
     &                 +sub2b_1(n,qq)
     &                  *msq2b_1(n,0,j,k)*xn/2._dp
          enddo    
c-- [n=4] (a,b,c) = (5,6,4)
          msq(4,j,k)  =msq(4,j,k)
     &                 +subab_c(4,gg)
     &                  *msqab_c(4,0,j,k)*xn/2._dp                
     &                 +subab_cv(4)
     &                  *msqab_cv(4,0,j,k)*xn/2._dp                
c--- choose n=1 which is (6,4,5)
          do n=3,3
          msq(18+n,j,k)=msq(18+n,j,k)+sub1b_2(n,qg)*(aveqg/aveqq)*
     &    (-(msq1b_2(n,1,+1,k)+msq1b_2(n,1,+2,k)+msq1b_2(n,1,+3,k)
     &      +msq1b_2(n,1,+4,k)+msq1b_2(n,1,+5,k)
     &      +msq1b_2(n,2,+1,k)+msq1b_2(n,2,+2,k)+msq1b_2(n,2,+3,k)
     &      +msq1b_2(n,2,+4,k)+msq1b_2(n,2,+5,k))/xn
     &     +(msq1b_2(n,0,+1,k)+msq1b_2(n,0,+2,k)+msq1b_2(n,0,+3,k)
     &      +msq1b_2(n,0,+4,k)+msq1b_2(n,0,+5,k))*2._dp*xn)
          enddo
          endif
c------ sub-sub-leading colour
          if ((colourchoice == 3) .or. (colourchoice == 0)) then
          do n=3,5,2
c-- (a,b,c) = (6,4,5) and (5,4,6)
          msq(12+n,j,k)=msq(12+n,j,k)
     &                 -(subbc_2(n,qq)+sub2c_b(n,qq))
     &                  *msq2c_b(n,0,j,k)*(xn+1._dp/xn)/2._dp
          enddo
c--- choose n=3 which is (6,4,5)
          do n=3,3
          msq(18+n,j,k)=msq(18+n,j,k)+sub1b_2(n,qg)*(aveqg/aveqq)*
     &      (msq1b_2(n,0,+1,k)+msq1b_2(n,0,+2,k)+msq1b_2(n,0,+3,k)
     &      +msq1b_2(n,0,+4,k)+msq1b_2(n,0,+5,k))*(-xn-1._dp/xn)
          enddo
          endif
c--- ANTIQUARK-GLUON contributions
      elseif ((j<0) .and. (k==0)) then
c------ leading colour
          if ((colourchoice == 1) .or. (colourchoice == 0)) then
          do n=1,2  
c-- (a,b,c) = (4,5,6) and (4,6,5)
          msq(n,j,k)   =subab_c(n,qq)
     &                  *msqab_c(n,pntr(c(n),a(n)),j,k)*xn/2._dp
          msq(12+n,j,k)=(subbc_2(n,gg)/2._dp+sub2c_b(n,gg))
     &                  *(msq2c_b(n,pntr(b(n),a(n)),j,k)
     &                   +msq2c_b(n,pntr(a(n),b(n)),j,k))*xn/2._dp
     &                 +subbc_2v(n)/2._dp
     &                  *(msqbc_2v(n,pntr(b(n),a(n)),j,k)
     &                   +msqbc_2v(n,pntr(a(n),b(n)),j,k))*xn/2._dp
     &                 +sub2c_bv(n)
     &                  *(msq2c_bv(n,pntr(b(n),a(n)),j,k)
     &                   +msq2c_bv(n,pntr(a(n),b(n)),j,k))*xn/2._dp
          msq(18+n,j,k)=sub1b_2(n,qq)
     &                  *msq1b_2(n,pntr(c(n),a(n)),j,k)*xn/2._dp
          msq(21+n,j,k)=sub2b_1(n,gg)
     &                  *msq2b_1(n,pntr(c(n),a(n)),j,k)*xn/2._dp
     &                 +sub2b_1v(n)
     &                  *msq2b_1v(n,pntr(c(n),a(n)),j,k)*xn/2._dp
          enddo
          do n=4,6,2
c-- (a,b,c) = (5,6,4) and (6,5,4)
          msq(n,j,k)   =subab_c(n,gg)/2._dp
     &                  *msqab_c(n,pntr(a(n),c(n)),j,k)*xn/2._dp
     &                 +subab_cv(n)/2._dp
     &                  *msqab_cv(n,pntr(a(n),c(n)),j,k)*xn/2._dp
          msq(6+n,j,k) =(sub1a_b(n,qq)+subba_1(n,gg)/2._dp)
     &                  *msq1a_b(n,pntr(c(n),b(n)),j,k)*xn/2._dp
     &                 +subba_1v(n)/2._dp
     &                  *msqba_1v(n,pntr(c(n),b(n)),j,k)*xn/2._dp
          enddo
          do n=3,5,2
c-- (a,b,c) = (6,4,5) and (5,4,6)
          msq(12+n,j,k)=(subbc_2(n,qq)+sub2c_b(n,gg))
     &                  *msq2c_b(n,pntr(b(n),a(n)),j,k)*xn/2._dp
     &                 +sub2c_bv(n)
     &                  *msq2c_bv(n,pntr(b(n),a(n)),j,k)*xn/2._dp
          enddo      
c--- additional (qg) collinear contributions
c--- choose n=1 which is (6,4,5)
          do n=3,3
          msq(21+n,j,k)=msq(21+n,j,k)+sub2b_1(n,qg)*
     &      (msq2b_1(n,1,j,+1)+msq2b_1(n,1,j,+2)+msq2b_1(n,1,j,+3)
     &      +msq2b_1(n,1,j,+4)+msq2b_1(n,1,j,+5)
     &      +msq2b_1(n,2,j,+1)+msq2b_1(n,2,j,+2)+msq2b_1(n,2,j,+3)
     &      +msq2b_1(n,2,j,+4)+msq2b_1(n,2,j,+5))*xn*(aveqg/aveqq)
          enddo
          endif
c------ sub-leading colour
          if ((colourchoice == 2) .or. (colourchoice == 0)) then
          do n=4,6,2
c-- (a,b,c) = (5,6,4) and (6,5,4)
          msq(6+n,j,k) =msq(6+n,j,k)
     &                 +(sub1a_b(n,qq)+subba_1(n,gg)/2._dp)
     &                  *msq1a_b(n,0,j,k)*xn/2._dp
     &                 +subba_1v(n)/2._dp
     &                  *msqba_1v(n,0,j,k)*xn/2._dp
          enddo
          do n=3,5,2
c-- (a,b,c) = (6,4,5) and (5,4,6)
          msq(6+n,j,k)=msq(6+n,j,k)
     &                 -(subba_1(n,qq)+sub1a_b(n,qq))
     &                  *(msq1a_b(n,1,j,k)+msq1a_b(n,2,j,k))/xn/2._dp
          msq(12+n,j,k)=msq(12+n,j,k)
     &                 +(subbc_2(n,qq)+sub2c_b(n,gg))
     &                  *msq2c_b(n,0,j,k)*xn/2._dp
     &                 +sub2c_bv(n)
     &                  *msq2c_bv(n,0,j,k)*xn/2._dp
          enddo      
          do n=1,2  
c-- (a,b,c) = (4,5,6) and (4,6,5)
          msq(n,j,k)   =msq(n,j,k)
     &                 +subab_c(n,qq)
     &                  *msqab_c(n,0,j,k)*xn/2._dp
          msq(18+n,j,k)=msq(18+n,j,k)
     &                 +sub1b_2(n,qq)
     &                  *msq1b_2(n,0,j,k)*xn/2._dp
          msq(21+n,j,k)=msq(21+n,j,k)
     &                 +sub2b_1(n,gg)
     &                  *msq2b_1(n,0,j,k)*xn/2._dp
     &                 +sub2b_1v(n)
     &                  *msq2b_1v(n,0,j,k)*xn/2._dp
          enddo    
c-- [n=4] (a,b,c) = (5,6,4)
          msq(4,j,k)  =msq(4,j,k)
     &                 +subab_c(4,gg)
     &                  *msqab_c(4,0,j,k)*xn/2._dp
     &                 +subab_cv(4)
     &                  *msqab_cv(4,0,j,k)*xn/2._dp
c--- choose n=1 which is (6,4,5)
          do n=3,3
          msq(21+n,j,k)=msq(21+n,j,k)+sub2b_1(n,qg)*(aveqg/aveqq)*
     &    (-(msq2b_1(n,1,j,+1)+msq2b_1(n,1,j,+2)+msq2b_1(n,1,j,+3)
     &      +msq2b_1(n,1,j,+4)+msq2b_1(n,1,j,+5)
     &      +msq2b_1(n,2,j,+1)+msq2b_1(n,2,j,+2)+msq2b_1(n,2,j,+3)
     &      +msq2b_1(n,2,j,+4)+msq2b_1(n,2,j,+5))/xn
     &     +(msq2b_1(n,0,j,+1)+msq2b_1(n,0,j,+2)+msq2b_1(n,0,j,+3)
     &      +msq2b_1(n,0,j,+4)+msq2b_1(n,0,j,+5))*2._dp*xn)
          enddo
          endif
c------ sub-sub-leading colour
          if ((colourchoice == 3) .or. (colourchoice == 0)) then
          do n=3,5,2
c-- (a,b,c) = (6,4,5) and (5,4,6)
          msq(6+n,j,k)=msq(6+n,j,k)
     &                 -(subba_1(n,qq)+sub1a_b(n,qq))
     &                  *msq1a_b(n,0,j,k)*(xn+1._dp/xn)/2._dp
          enddo
c--- additional (qg) collinear contributions
c--- choose n=3 which is (6,4,5)
          do n=3,3
          msq(21+n,j,k)=msq(21+n,j,k)+sub2b_1(n,qg)*(aveqg/aveqq)*
     &      (msq2b_1(n,0,j,+1)+msq2b_1(n,0,j,+2)+msq2b_1(n,0,j,+3)
     &      +msq2b_1(n,0,j,+4)+msq2b_1(n,0,j,+5))*(-xn-1._dp/xn)
          enddo          
          endif
      endif

      enddo
      enddo

      endif

      
      if (Qflag) then

c---arguments of dipsx_new:
c---    1  dipole number
c---    2  momentum
c---    3  emitter
c---    4  emitted
c---    4  spectator
c---    5  interference-free subtraction, equiv. to AP kernel for qq,qg
c---    6  correlation piece of subtraction, relevant only for gg,gq
c---    8  lowest order matrix elements at rescaled momentum, msq(j,k)
c---    9  lowest order matrix elements at rescaled momentum
c---        with emitter contracted with appropriate vector, msqv(j,k)
c---   10  appropriate lowest order calculating routine for msq
c---   11  appropriate lowest order calculating routine for msqv
c---   12  4-quark contribution to lowest order matrix elements squared
c---   13  lowest order matrix elements with 4 indices, msqx(j,k,l,m)
c---         Sum_{l,m} msqx(j,k,l,m) = msq(j,k)
c---   14  2-quark contribution to lowest order matrix elements squared,
c---        separated by colours
c---   14  2-quark contribution to lowest order matrix elements squared,
c---        separated by colours, contracted with appropriate vector
c---   15  lowest order matrix elements with 4 indices and
c----        contracted with appropriate vector, msqvx(j,k,l,m)
c---         Sum_{l,m} msqvx(j,k,l,m) = msqv(j,k)

c--- calculate all the dipoles
c--- the dipole number relates the matrix elements to the transformed
c--- momenta, used in realint to perform cuts and clustering

c--- final-final
      call dipsx_new(2,p,4,6,5,sub46_5,sub46_5v,msq46_5,msq46_5v,
     & qqb_gam2jx_new,qqb_gam2j_gvecx_new,m46_5,m46_5x,m46_5g,m46_5vg,m46_5vx)
      call dipsx_new(4,p,5,6,4,sub56_4,sub56_4v,msq56_4,dummyv,
     & qqb_gam2jx_new,qqb_gam2j_gvecx_new,m56_4,m56_4x,m56_4g,m56_4vg,m56_4vx)
c--- dipole added for consistency with Gflag piece
      call dipsx_new(5,p,5,4,6,sub54_6,sub54_6v,msq54_6,dummyv,
     & qqb_gam2jx_new,qqb_gam2j_gvecx_new,m54_6,m54_6x,m54_6g,m54_6vg,m54_6vx)

c--- now the basic initial final and final initial
c--- second call for dipole 9,12 etc  only supplies new values for
c--- sub..
      call dipsx_new(7,p,1,4,5,sub14_5,dsubv,msq14_5,dummyv,
     & qqb_gam2jx_new,donothing_gvecx_new,m14_5,m14_5x,m14_5g,mvg,mvxg)
      call dipsx_new(7,p,5,4,1,sub54_1,sub54_1v,dummy,dummyv,
     & qqb_gam2jx_new,qqb_gam2j_gvecx_new,mqq,msqx,m54_1g,m54_1vg,m54_1vx)

      call dipsx_new(9,p,1,6,4,sub16_4,dsubv,msq16_4,dummyv,
     & qqb_gam2jx_new,donothing_gvecx_new,m16_4,m16_4x,mg,mvg,mvxg)
      call dipsx_new(9,p,4,6,1,sub46_1,sub46_1v,dummy,dummyv,
     & qqb_gam2jx_new,qqb_gam2j_gvecx_new,mqq,msqx,m46_1g,m46_1vg,m46_1vx)

      call dipsx_new(12,p,1,6,5,sub16_5,dsubv,msq16_5,dummyv,
     & qqb_gam2jx_new,donothing_gvecx_new,m16_5,m16_5x,mg,mvg,mvxg)
      call dipsx_new(12,p,5,6,1,sub56_1,sub56_1v,dummy,dummyv,
     & qqb_gam2jx_new,qqb_gam2j_gvecx_new,mqq,msqx,m56_1g,m56_1vg,m56_1vx)

      call dipsx_new(13,p,2,6,5,sub26_5,dsubv,msq26_5,dummyv,
     & qqb_gam2jx_new,donothing_gvecx_new,m26_5,m26_5x,mg,mvg,mvxg)
      call dipsx_new(13,p,5,6,2,sub56_2,sub56_2v,dummy,dummyv,
     & qqb_gam2jx_new,qqb_gam2j_gvecx_new,mqq,msqx,m56_2g,m56_2vg,m56_2vx)
 
      call dipsx_new(17,p,2,6,4,sub26_4,dsubv,msq26_4,dummyv,
     & qqb_gam2jx_new,donothing_gvecx_new,m26_4,m26_4x,mg,mvg,mvxg)
      call dipsx_new(17,p,4,6,2,sub46_2,sub46_2v,dummy,dummyv,
     & qqb_gam2jx_new,qqb_gam2j_gvecx_new,mqq,msqx,m46_2g,m46_2vg,m46_2vx)

      call dipsx_new(18,p,2,4,5,sub24_5,dsubv,msq24_5,dummyv,
     & qqb_gam2jx_new,donothing_gvecx_new,m24_5,m24_5x,m24_5g,mvg,mvxg)
      call dipsx_new(18,p,5,4,2,sub54_2,sub54_2v,dummy,dummyv,
     & qqb_gam2jx_new,qqb_gam2j_gvecx_new,mqq,msqx,m54_2g,m54_2vg,m54_2vx)
      
c--- calculate all the initial-initial dipoles
      call dipsx_new(19,p,1,5,2,sub15_2,sub15_2v,msq15_2,msq15_2v,
     & qqb_gam2jx_new,qqb_gam2j_gvecx_new,m15_2,m15_2x,m15_2g,m15_2vg,m15_2vx)
      call dipsx_new(20,p,1,6,2,sub16_2,sub16_2v,msq16_2,msq16_2v,
     & qqb_gam2jx_new,qqb_gam2j_gvecx_new,m16_2,m16_2x,m16_2g,m16_2vg,m16_2vx)
      call dipsx_new(21,p,1,4,2,sub14_2,sub14_2v,msq14_2,msq14_2v,
     & qqb_gam2jx_new,qqb_gam2j_gvecx_new,m14_2,m14_2x,m14_2g,m14_2vg,m14_2vx)
      call dipsx_new(22,p,2,5,1,sub25_1,sub25_1v,msq25_1,msq25_1v,
     & qqb_gam2jx_new,qqb_gam2j_gvecx_new,m25_1,m25_1x,m25_1g,m25_1vg,m25_1vx)
      call dipsx_new(23,p,2,6,1,sub26_1,sub26_1v,msq26_1,dummyv,
     & qqb_gam2jx_new,qqb_gam2j_gvecx_new,m26_1,m26_1x,m26_1g,m26_1vg,m26_1vx)
      call dipsx_new(24,p,2,4,1,sub24_1,sub24_1v,msq24_1,dummyv,
     & qqb_gam2jx_new,qqb_gam2j_gvecx_new,m24_1,m24_1x,m24_1g,m24_1vg,m24_1vx)

c--- fill the dipole contributions
      do j=-nf,nf
      do k=-nf,nf
      
      if ((j > 0) .and. (k>0)) then
c---QQ
      msq(2,j,k)=msq(2,j,k)
     & +sub46_5(qq)
     & *((xn+1._dp/xn)*m46_5(0,j,k)+two*(m46_5(1,j,k)+m46_5(2,j,k))/xn)
      msq(4,j,k)=msq(4,j,k)
     & +sub56_4(qq)
     & *((xn+1._dp/xn)*m56_4(0,j,k)+two*(m56_4(1,j,k)+m56_4(2,j,k))/xn)
      msq(9,j,k)=msq(9,j,k)
     & +sub16_4(qq)
     & *((xn-two/xn)*m16_4(2,j,k)-m16_4(0,j,k)/xn-m16_4(1,j,k)/xn)
     & +sub46_1(qq)
     & *((xn-two/xn)*m16_4(2,j,k)-m16_4(0,j,k)/xn-m16_4(1,j,k)/xn)
      msq(12,j,k)=msq(12,j,k)
     & +sub16_5(qq)
     & *((xn-two/xn)*m16_5(1,j,k)-m16_5(0,j,k)/xn-m16_5(2,j,k)/xn)
     & +sub56_1(qq)
     & *((xn-two/xn)*m16_5(1,j,k)-m16_5(0,j,k)/xn-m16_5(2,j,k)/xn)
      msq(13,j,k)=msq(13,j,k)
     & +sub26_5(qq)
     & *((xn-two/xn)*m26_5(2,j,k)-m26_5(0,j,k)/xn-m26_5(1,j,k)/xn)
     & +sub56_2(qq)
     & *((xn-two/xn)*m26_5(2,j,k)-m26_5(0,j,k)/xn-m26_5(1,j,k)/xn)
      msq(17,j,k)=msq(17,j,k)
     & +sub26_4(qq)
     & *((xn-two/xn)*m26_4(1,j,k)-m26_4(0,j,k)/xn-m26_4(2,j,k)/xn)
     & +sub46_2(qq)
     & *((xn-two/xn)*m26_4(1,j,k)-m26_4(0,j,k)/xn-m26_4(2,j,k)/xn)
      msq(20,j,k)=msq(20,j,k)
     & +sub16_2(qq)
     & *((xn+1._dp/xn)*m16_2(0,j,k)+two*(m16_2(1,j,k)+m16_2(2,j,k))/xn)
      msq(23,j,k)=msq(23,j,k)
     & +sub26_1(qq)
     & *((xn+1._dp/xn)*m26_1(0,j,k)+two*(m26_1(1,j,k)+m26_1(2,j,k))/xn)

      if (j==k) then
      msq(22,j,k)=msq(22,j,k)+0.5_dp*sub25_1(gq)*(aveqq/aveqg)
     &    *(m25_1x(0,pp(jj(j),0,jj(j),0))+m25_1x(1,pp(jj(j),0,jj(j),0))+m25_1x(2,pp(jj(j),0,jj(j),0)))
     &    +0.5_dp*sub25_1v*(aveqq/aveqg)*m25_1vx(pp(jj(j),0,jj(j),0))
      msq(24,j,k)=msq(24,j,k)+0.5_dp*sub24_1(gq)*(aveqq/aveqg)
     &    *(m24_1x(0,pp(jj(j),0,jj(j),0))+m24_1x(1,pp(jj(j),0,jj(j),0))+m24_1x(2,pp(jj(j),0,jj(j),0)))
     &    +0.5_dp*sub24_1v*(aveqq/aveqg)*(+m24_1vx(pp(jj(j),0,jj(j),0)))
      msq(21,j,k)=msq(21,j,k)+0.5_dp*sub14_2(gq)*(aveqq/aveqg)
     &    *(m14_2x(0,pp(0,kk(k),kk(k),0))+m14_2x(1,pp(0,kk(k),kk(k),0))+m14_2x(2,pp(0,kk(k),kk(k),0)))
     &    +0.5_dp*sub14_2v*(aveqq/aveqg)*m14_2vx(pp(0,kk(k),kk(k),0))
      msq(19,j,k)=msq(19,j,k)+0.5_dp*sub15_2(gq)*(aveqq/aveqg)
     &    *(m15_2x(0,pp(0,kk(k),kk(k),0))+m15_2x(1,pp(0,kk(k),kk(k),0))+m15_2x(2,pp(0,kk(k),kk(k),0)))
     &    +0.5_dp*sub15_2v*(aveqq/aveqg)*m15_2vx(pp(0,kk(k),kk(k),0))

      elseif (j .ne. k) then
      msq(22,j,k)=msq(22,j,k)+sub25_1(gq)*(aveqq/aveqg)
     &    *(m25_1x(0,pp(jj(j),0,jj(j),0))+m25_1x(1,pp(jj(j),0,jj(j),0))+m25_1x(2,pp(jj(j),0,jj(j),0)))
     &    +sub25_1v*(aveqq/aveqg)*m25_1vx(pp(jj(j),0,jj(j),0))
      msq(21,j,k)=msq(21,j,k)+sub14_2(gq)*(aveqq/aveqg)
     &    *(m14_2x(0,pp(0,kk(k),kk(k),0))+m14_2x(1,pp(0,kk(k),kk(k),0))+m14_2x(2,pp(0,kk(k),kk(k),0)))
     &    +sub14_2v*(aveqq/aveqg)*m14_2vx(pp(0,kk(k),kk(k),0))
      endif

      elseif ((j < 0) .and. (k<0)) then
c---QbarQbar
      msq(2,j,k)=msq(2,j,k)
     & +sub46_5(qq)
     & *((xn+1._dp/xn)*m46_5(0,j,k)+two*(m46_5(1,j,k)+m46_5(2,j,k))/xn)
      msq(4,j,k)=msq(4,j,k)
     & +sub56_4(qq)
     & *((xn+1._dp/xn)*m56_4(0,j,k)+two*(m56_4(1,j,k)+m56_4(2,j,k))/xn)
      msq(9,j,k)=msq(9,j,k)
     & +sub16_4(qq)
     & *((xn-two/xn)*m16_4(2,j,k)-m16_4(0,j,k)/xn-m16_4(1,j,k)/xn)
     & +sub46_1(qq)
     & *((xn-two/xn)*m16_4(2,j,k)-m16_4(0,j,k)/xn-m16_4(1,j,k)/xn)
      msq(12,j,k)=msq(12,j,k)
     & +sub16_5(qq)
     & *((xn-two/xn)*m16_5(1,j,k)-m16_5(0,j,k)/xn-m16_5(2,j,k)/xn)
     & +sub56_1(qq)
     & *((xn-two/xn)*m16_5(1,j,k)-m16_5(0,j,k)/xn-m16_5(2,j,k)/xn)
      msq(13,j,k)=msq(13,j,k)
     & +sub26_5(qq)
     & *((xn-two/xn)*m26_5(2,j,k)-m26_5(0,j,k)/xn-m26_5(1,j,k)/xn)
     & +sub56_2(qq)
     & *((xn-two/xn)*m26_5(2,j,k)-m26_5(0,j,k)/xn-m26_5(1,j,k)/xn)
      msq(17,j,k)=msq(17,j,k)
     & +sub26_4(qq)
     & *((xn-two/xn)*m26_4(1,j,k)-m26_4(0,j,k)/xn-m26_4(2,j,k)/xn)
     & +sub46_2(qq)
     & *((xn-two/xn)*m26_4(1,j,k)-m26_4(0,j,k)/xn-m26_4(2,j,k)/xn)
      msq(20,j,k)=msq(20,j,k)
     & +sub16_2(qq)
     & *((xn+1._dp/xn)*m16_2(0,j,k)+two*(m16_2(1,j,k)+m16_2(2,j,k))/xn)
      msq(23,j,k)=msq(23,j,k)
     & +sub26_1(qq)
     & *((xn+1._dp/xn)*m26_1(0,j,k)+two*(m26_1(1,j,k)+m26_1(2,j,k))/xn)

      if (j==k) then
      msq(22,j,k)=msq(22,j,k)+0.5_dp*sub25_1(gq)*(aveqq/aveqg)
     &    *(m25_1x(0,pp(jj(j),0,jj(j),0))+m25_1x(1,pp(jj(j),0,jj(j),0))+m25_1x(2,pp(jj(j),0,jj(j),0)))
     &    +0.5_dp*sub25_1v*(aveqq/aveqg)*m25_1vx(pp(jj(j),0,jj(j),0))
      msq(24,j,k)=msq(24,j,k)+0.5_dp*sub24_1(gq)*(aveqq/aveqg)
     &    *(m24_1x(0,pp(jj(j),0,jj(j),0))+m24_1x(1,pp(jj(j),0,jj(j),0))+m24_1x(2,pp(jj(j),0,jj(j),0)))
     &    +0.5_dp*sub24_1v*(aveqq/aveqg)*(+m24_1vx(pp(jj(j),0,jj(j),0)))
      msq(21,j,k)=msq(21,j,k)+0.5_dp*sub14_2(gq)*(aveqq/aveqg)
     &    *(m14_2x(0,pp(0,kk(k),kk(k),0))+m14_2x(1,pp(0,kk(k),kk(k),0))+m14_2x(2,pp(0,kk(k),kk(k),0)))
     &    +0.5_dp*sub14_2v*(aveqq/aveqg)*m14_2vx(pp(0,kk(k),kk(k),0))
      msq(19,j,k)=msq(19,j,k)+0.5_dp*sub15_2(gq)*(aveqq/aveqg)
     &    *(m15_2x(0,pp(0,kk(k),kk(k),0))+m15_2x(1,pp(0,kk(k),kk(k),0))+m15_2x(2,pp(0,kk(k),kk(k),0)))
     &    +0.5_dp*sub15_2v*(aveqq/aveqg)*m15_2vx(pp(0,kk(k),kk(k),0))

      elseif (j .ne. k) then
      msq(22,j,k)=msq(22,j,k)+sub25_1(gq)*(aveqq/aveqg)
     &    *(m25_1x(0,pp(jj(j),0,jj(j),0))+m25_1x(1,pp(jj(j),0,jj(j),0))+m25_1x(2,pp(jj(j),0,jj(j),0)))
     &    +sub25_1v*(aveqq/aveqg)*m25_1vx(pp(jj(j),0,jj(j),0))
      msq(21,j,k)=msq(21,j,k)+sub14_2(gq)*(aveqq/aveqg)
     &    *(m14_2x(0,pp(0,kk(k),kk(k),0))+m14_2x(1,pp(0,kk(k),kk(k),0))+m14_2x(2,pp(0,kk(k),kk(k),0)))
     &    +sub14_2v*(aveqq/aveqg)*m14_2vx(pp(0,kk(k),kk(k),0))
      endif


      elseif ((j > 0) .and. (k<0)) then
c---QQbar
      msq(2,j,k)=msq(2,j,k)
     & +sub46_5(qq)
     & *((xn-two/xn)*m46_5(2,j,k)-m46_5(0,j,k)/xn-m46_5(1,j,k)/xn)
      msq(4,j,k)=msq(4,j,k)
     & +sub56_4(qq)
     & *((xn-two/xn)*m56_4(2,j,k)-m56_4(0,j,k)/xn-m56_4(1,j,k)/xn)
      msq(9,j,k)=msq(9,j,k)
     & +sub16_4(qq)
     & *((xn-two/xn)*m16_4(1,j,k)-m16_4(0,j,k)/xn-m16_4(2,j,k)/xn)
     & +sub46_1(qq)
     & *((xn-two/xn)*m16_4(1,j,k)-m16_4(0,j,k)/xn-m16_4(2,j,k)/xn)
      msq(12,j,k)=msq(12,j,k)
     & +sub16_5(qq)
     & *((xn+1._dp/xn)*m16_5(0,j,k)+two*(m16_5(1,j,k)+m16_5(2,j,k))/xn)
     & +sub56_1(qq)
     & *((xn+1._dp/xn)*m16_5(0,j,k)+two*(m16_5(1,j,k)+m16_5(2,j,k))/xn)
      msq(13,j,k)=msq(13,j,k)
     & +sub26_5(qq)
     & *((xn-two/xn)*m26_5(1,j,k)-m26_5(0,j,k)/xn-m26_5(2,j,k)/xn)
     & +sub56_2(qq)
     & *((xn-two/xn)*m26_5(1,j,k)-m26_5(0,j,k)/xn-m26_5(2,j,k)/xn)
      msq(17,j,k)=msq(17,j,k)
     & +sub26_4(qq)
     & *((xn+1._dp/xn)*m26_4(0,j,k)+two*(m26_4(1,j,k)+m26_4(2,j,k))/xn)
     & +sub46_2(qq)
     & *((xn+1._dp/xn)*m26_4(0,j,k)+two*(m26_4(1,j,k)+m26_4(2,j,k))/xn)
      msq(20,j,k)=msq(20,j,k)
     & +sub16_2(qq)
     & *((xn-two/xn)*m16_2(2,j,k)-m16_2(0,j,k)/xn-m16_2(1,j,k)/xn)
      msq(23,j,k)=msq(23,j,k)
     & +sub26_1(qq)
     & *((xn-two/xn)*m26_1(2,j,k)-m26_1(0,j,k)/xn-m26_1(1,j,k)/xn)


      if (j .ne. -k) then
      msq(22,j,k)=msq(22,j,k)+sub25_1(gq)*(aveqq/aveqg)
     &    *(m25_1x(0,pp(jj(j),0,jj(j),0))+m25_1x(1,pp(jj(j),0,jj(j),0))+m25_1x(2,pp(jj(j),0,jj(j),0)))
     &    +sub25_1v*(aveqq/aveqg)*m25_1vx(pp(jj(j),0,jj(j),0))
      msq(21,j,k)=msq(21,j,k)+sub14_2(gq)*(aveqq/aveqg)
     &    *(m14_2x(0,pp(0,kk(k),kk(k),0))+m14_2x(1,pp(0,kk(k),kk(k),0))+m14_2x(2,pp(0,kk(k),kk(k),0)))
     &    +sub14_2v*(aveqq/aveqg)*m14_2vx(pp(0,kk(k),kk(k),0))
      elseif (j==-k) then
c--- combination of final-initial and final-final for consistency with Gflag
      msq( 5,j,k)=msq( 5,j,k)+real(nf,dp)
     & *(sub54_6(gq)*(m54_6g(0,j,k)+m54_6g(1,j,k)+m54_6g(2,j,k))
     &  -sub54_6v*(m54_6vg(0,j,k)+m54_6vg(1,j,k)+m54_6vg(2,j,k)))
c--- 8/6/09: symmetrized 54 between _1 and _2
      msq( 7,j,k)=msq( 7,j,k)+half*real(nf,dp)
     & *(sub54_1(gq)*(m14_5g(0,j,k)+m14_5g(1,j,k)+m14_5g(2,j,k))
     &  -sub54_1v*(m54_1vg(0,j,k)+m54_1vg(1,j,k)+m54_1vg(2,j,k)))
      msq(18,j,k)=msq(18,j,k)+half*real(nf,dp)
     & *(sub54_2(gq)*(m24_5g(0,j,k)+m24_5g(1,j,k)+m24_5g(2,j,k))
     &  -sub54_2v*(m54_2vg(0,j,k)+m54_2vg(1,j,k)+m54_2vg(2,j,k)))

      msq(22,j,k)=msq(22,j,k)+sub25_1(gq)*(aveqq/aveqg)
     &    *(m25_1x(0,pp(jj(j),0,jj(j),0))+m25_1x(1,pp(jj(j),0,jj(j),0))+m25_1x(2,pp(jj(j),0,jj(j),0)))
     &    +sub25_1v*(aveqq/aveqg)*m25_1vx(pp(jj(j),0,jj(j),0))
      msq(21,j,k)=msq(21,j,k)+sub14_2(gq)*(aveqq/aveqg)
     &    *(m14_2x(0,pp(0,kk(k),kk(k),0))+m14_2x(1,pp(0,kk(k),kk(k),0))+m14_2x(2,pp(0,kk(k),kk(k),0)))
     &    +sub14_2v*(aveqq/aveqg)*m14_2vx(pp(0,kk(k),kk(k),0))

      endif

      elseif ((j < 0) .and. (k>0)) then
c---QbarQ
      msq(2,j,k)=msq(2,j,k)
     & +sub46_5(qq)
     & *((xn-two/xn)*m46_5(2,j,k)-m46_5(0,j,k)/xn-m46_5(1,j,k)/xn)
      msq(4,j,k)=msq(4,j,k)
     & +sub56_4(qq)
     & *((xn-two/xn)*m56_4(2,j,k)-m56_4(0,j,k)/xn-m56_4(1,j,k)/xn)
      msq(9,j,k)=msq(9,j,k)
     & +sub16_4(qq)
     & *((xn+1._dp/xn)*m16_4(0,j,k)+two*(m16_4(1,j,k)+m16_4(2,j,k))/xn)
     & +sub46_1(qq)
     & *((xn+1._dp/xn)*m16_4(0,j,k)+two*(m16_4(1,j,k)+m16_4(2,j,k))/xn)
      msq(12,j,k)=msq(12,j,k)
     & +sub16_5(qq)
     & *((xn-two/xn)*m16_5(1,j,k)-m16_5(0,j,k)/xn-m16_5(2,j,k)/xn)
     & +sub56_1(qq)
     & *((xn-two/xn)*m16_5(1,j,k)-m16_5(0,j,k)/xn-m16_5(2,j,k)/xn)
      msq(13,j,k)=msq(13,j,k)
     & +sub26_5(qq)
     & *((xn+1._dp/xn)*m26_5(0,j,k)+two*(m26_5(1,j,k)+m26_5(2,j,k))/xn)
     & +sub56_2(qq)
     & *((xn+1._dp/xn)*m26_5(0,j,k)+two*(m26_5(1,j,k)+m26_5(2,j,k))/xn)
      msq(17,j,k)=msq(17,j,k)
     & +sub26_4(qq)
     & *((xn-two/xn)*m26_4(1,j,k)-m26_4(0,j,k)/xn-m26_4(2,j,k)/xn)
     & +sub46_2(qq)
     & *((xn-two/xn)*m26_4(1,j,k)-m26_4(0,j,k)/xn-m26_4(2,j,k)/xn)
      msq(20,j,k)=msq(20,j,k)
     & +sub16_2(qq)
     & *((xn-two/xn)*m16_2(2,j,k)-m16_2(0,j,k)/xn-m16_2(1,j,k)/xn)
      msq(23,j,k)=msq(23,j,k)
     & +sub26_1(qq)
     & *((xn-two/xn)*m26_1(2,j,k)-m26_1(0,j,k)/xn-m26_1(1,j,k)/xn)

      if (-j .ne. k) then
      msq(24,j,k)=msq(24,j,k)+sub24_1(gq)*(aveqq/aveqg)
     &    *(m24_1x(0,pp(jj(j),0,jj(j),0))+m24_1x(1,pp(jj(j),0,jj(j),0))+m24_1x(2,pp(jj(j),0,jj(j),0)))
     &    +sub24_1v*(aveqq/aveqg)*m24_1vx(pp(jj(j),0,jj(j),0))
      msq(19,j,k)=msq(19,j,k)+sub15_2(gq)*(aveqq/aveqg)
     &    *(m15_2x(0,pp(0,kk(k),kk(k),0))+m15_2x(1,pp(0,kk(k),kk(k),0))+m15_2x(2,pp(0,kk(k),kk(k),0)))
     &    +sub15_2v*(aveqq/aveqg)*m15_2vx(pp(0,kk(k),kk(k),0))
      elseif (-j==k) then
c--- combination of final-initial and final-final for consistency with Gflag
      msq( 5,j,k)=msq( 5,j,k)+real(nf,dp)
     & *(sub54_6(gq)*(m54_6g(0,j,k)+m54_6g(1,j,k)+m54_6g(2,j,k))
     &  -sub54_6v*(m54_6vg(0,j,k)+m54_6vg(1,j,k)+m54_6vg(2,j,k)))
c--- 8/5/09: symmetrized 54 between _1 and _2
      msq( 7,j,k)=msq( 7,j,k)+half*real(nf,dp)
     & *(sub54_1(gq)*(m14_5g(0,j,k)+m14_5g(1,j,k)+m14_5g(2,j,k))
     &  -sub54_1v*(m54_1vg(0,j,k)+m54_1vg(1,j,k)+m54_1vg(2,j,k)))
      msq(18,j,k)=msq(18,j,k)+half*real(nf,dp)
     & *(sub54_2(gq)*(m24_5g(0,j,k)+m24_5g(1,j,k)+m24_5g(2,j,k))
     &  -sub54_2v*(m54_2vg(0,j,k)+m54_2vg(1,j,k)+m54_2vg(2,j,k)))

      msq(24,j,k)=msq(24,j,k)+sub24_1(gq)*(aveqq/aveqg)
     &    *(m24_1x(0,pp(jj(j),0,jj(j),0))+m24_1x(1,pp(jj(j),0,jj(j),0))+m24_1x(2,pp(jj(j),0,jj(j),0)))
     &    +sub24_1v*(aveqq/aveqg)*m24_1vx(pp(jj(j),0,jj(j),0))
      msq(19,j,k)=msq(19,j,k)+sub15_2(gq)*(aveqq/aveqg)
     &    *(m15_2x(0,pp(0,kk(k),kk(k),0))+m15_2x(1,pp(0,kk(k),kk(k),0))+m15_2x(2,pp(0,kk(k),kk(k),0)))
     &    +sub15_2v*(aveqq/aveqg)*m15_2vx(pp(0,kk(k),kk(k),0))
      endif

      elseif (j==0) then
c---------G-Q and G-Qbar
        if    (k > 0) then
          msq( 7,j,k)=msq( 7,j,k)+(real(nf-1,dp)+0.5_dp)*half
     &    *(sub54_1(gq)
     &    *(m14_5x(1,pp(0,kk(k),0,kk(k)))+m14_5x(2,pp(0,kk(k),0,kk(k))))
     &    -sub54_1v
     &    *(m54_1vg(1,0,kk(k))+m54_1vg(2,0,kk(k))))
          msq(18,j,k)=msq(18,j,k)+(real(nf-1,dp)+0.5_dp)*half
     &    *(sub54_2(gq)
     &    *(m24_5x(0,pp(0,kk(k),0,kk(k)))+m24_5x(2,pp(0,kk(k),0,kk(k))))
     &    -sub54_2v
     &    *(m54_2vg(0,0,kk(k))+m54_2vg(2,0,kk(k))))
          msq( 5,j,k)=msq( 5,j,k)+(real(nf-1,dp)+0.5_dp)*half
     &    *(sub54_6(gq)
     &    *(m54_6x(0,pp(0,kk(k),0,kk(k)))+m54_6x(1,pp(0,kk(k),0,kk(k))))
     &    -sub54_6v
     &    *(m54_6vg(0,0,kk(k))+m54_6vg(1,0,kk(k))))

          msq(12,j,k)=msq(12,j,k)+0.5_dp*half*sub56_1(gq)
     &    *(m16_5x(1,pp(0,kk(k),kk(k),0))+m16_5x(2,pp(0,kk(k),kk(k),0)))
     &                           -0.5_dp*half*sub56_1v
     &    *(m56_1vg(1,0,kk(k))+m56_1vg(2,0,kk(k)))
          msq(13,j,k)=msq(13,j,k)+0.5_dp*half*sub56_2(gq)
     &    *(m26_5x(0,pp(0,kk(k),kk(k),0))+m26_5x(2,pp(0,kk(k),kk(k),0)))
     &                           -0.5_dp*half*sub56_2v
     &    *(m56_2vg(0,0,kk(k))+m56_2vg(2,0,kk(k)))
          msq( 4,j,k)=msq( 4,j,k)+0.5_dp*half*sub56_4(gq)
     &    *(m56_4x(0,pp(0,kk(k),kk(k),0))+m56_4x(1,pp(0,kk(k),kk(k),0)))
     &                           -0.5_dp*half*sub56_4v
     &    *(m56_4vg(0,0,kk(k))+m56_4vg(1,0,kk(k)))

          if (kk(k) == 1) then
          msq(19,j,k)=msq(19,j,k)+(xn-one/xn)*sub15_2(qg)*aveqg/aveqq*(
     &     (m15_2x(0,pp(1,1,1,1))+m15_2x(1,pp(1,1,1,1))+m15_2x(2,pp(1,1,1,1)))
     &    +2._dp*(m15_2x(0,pp(3,1,3,1))+m15_2x(1,pp(3,1,3,1))+m15_2x(2,pp(3,1,3,1)))
     &    +2._dp*(m15_2x(0,pp(2,1,2,1))+m15_2x(1,pp(2,1,2,1))+m15_2x(2,pp(2,1,2,1))))
 
          msq(20,j,k)=msq(20,j,k)+(xn-one/xn)*sub16_2(qg)*aveqg/aveqq
     &   *(0.5_dp
     &   *(m16_2x(0,pp(-1,1,1,-1))+m16_2x(1,pp(-1,1,1,-1))+m16_2x(2,pp(-1,1,1,-1)))
     &   +2._dp
     &   *(m16_2x(0,pp(-1,1,2,-2))+m16_2x(1,pp(-1,1,2,-2))+m16_2x(2,pp(-1,1,2,-2)))
     &   +2._dp
     &   *(m16_2x(0,pp(-1,1,3,-3))+m16_2x(1,pp(-1,1,3,-3))+m16_2x(2,pp(-1,1,3,-3))))

          msq(21,j,k)=msq(21,j,k)+(xn-one/xn)*sub14_2(qg)*aveqg/aveqq
     &   *(0.5_dp
     &   *(m14_2x(0,pp(-1,1,-1,1))+m14_2x(1,pp(-1,1,-1,1))+m14_2x(2,pp(-1,1,-1,1)))
     &   +2._dp
     &   *(m14_2x(0,pp(-2,1,-2,1))+m14_2x(1,pp(-2,1,-2,1))+m14_2x(2,pp(-2,1,-2,1)))
     &   +2._dp
     &   *(m14_2x(0,pp(-3,1,-3,1))+m14_2x(1,pp(-3,1,-3,1))+m14_2x(2,pp(-3,1,-3,1))))

          msq(23,j,k)=msq(23,j,k)+(aveqg/avegg)*(
     &    +2.5_dp*(sub26_1(gq)
     &    *(m26_1x(0,pp(0,0,1,-1))+m26_1x(1,pp(0,0,1,-1))+m26_1x(2,pp(0,0,1,-1)))
     &    +sub26_1v*m26_1vx(pp(0,0,1,-1)))
     &    +2._dp*(sub26_1(gq)
     &    *(m26_1x(0,pp(0,0,2,-2))+m26_1x(1,pp(0,0,2,-2))+m26_1x(2,pp(0,0,2,-2)))
     &    +sub26_1v*m26_1vx(pp(0,0,2,-2))))

          msq(24,j,k)=msq(24,j,k)+0.5_dp*(aveqg/avegg)*(sub24_1(gq)
     &    *(m24_1x(0,pp(0,0,-1,1))+m24_1x(1,pp(0,0,-1,1))+m24_1x(2,pp(0,0,-1,1)))
     &                     +sub24_1v*m24_1vx(pp(0,0,-1,1)))

         endif
         if (kk(k) == 2) then
          msq(23,j,k)=msq(23,j,k)+(aveqg/avegg)*(
     &    +3._dp*(sub26_1(gq)
     &    *(m26_1x(0,pp(0,0,1,-1))+m26_1x(1,pp(0,0,1,-1))+m26_1x(2,pp(0,0,1,-1)))
     &    +sub26_1v*m26_1vx(pp(0,0,1,-1)))
     &    +1.5_dp*(sub26_1(gq)
     &    *(m26_1x(0,pp(0,0,2,-2))+m26_1x(1,pp(0,0,2,-2))+m26_1x(2,pp(0,0,2,-2)))
     &    +sub26_1v*m26_1vx(pp(0,0,2,-2))))

          msq(24,j,k)=msq(24,j,k)+0.5_dp*(aveqg/avegg)*(sub24_1(gq)
     &    *(m24_1x(0,pp(0,0,-2,2))+m24_1x(1,pp(0,0,-2,2))+m24_1x(2,pp(0,0,-2,2)))
     &                     +sub24_1v*m24_1vx(pp(0,0,-2,2)))

          msq(19,j,k)=msq(19,j,k)+(xn-one/xn)*sub15_2(qg)*aveqg/aveqq*(
     &     3._dp*(m15_2x(0,pp(1,2,1,2))+m15_2x(1,pp(1,2,1,2))+m15_2x(2,pp(1,2,1,2)))
     &    +(m15_2x(0,pp(4,2,4,2))+m15_2x(1,pp(4,2,4,2))+m15_2x(2,pp(4,2,4,2)))
     &    +(m15_2x(0,pp(2,2,2,2))+m15_2x(1,pp(2,2,2,2))+m15_2x(2,pp(2,2,2,2))))
 
          msq(20,j,k)=msq(20,j,k)+(xn-one/xn)*sub16_2(qg)*aveqg/aveqq
     &   *(0.5_dp
     &   *(m16_2x(0,pp(-2,2,2,-2))+m16_2x(1,pp(-2,2,2,-2))+m16_2x(2,pp(-2,2,2,-2)))
     &   +3._dp
     &   *(m16_2x(0,pp(-2,2,1,-1))+m16_2x(1,pp(-2,2,1,-1))+m16_2x(2,pp(-2,2,1,-1)))
     &   +1._dp
     &   *(m16_2x(0,pp(-2,2,4,-4))+m16_2x(1,pp(-2,2,4,-4))+m16_2x(2,pp(-2,2,4,-4))))

          msq(21,j,k)=msq(21,j,k)+(xn-one/xn)*sub14_2(qg)*aveqg/aveqq
     &   *(3._dp
     &   *(m14_2x(0,pp(-1,2,-1,2))+m14_2x(1,pp(-1,2,-1,2))+m14_2x(2,pp(-1,2,-1,2)))
     &   +0.5_dp
     &   *(m14_2x(0,pp(-2,2,-2,2))+m14_2x(1,pp(-2,2,-2,2))+m14_2x(2,pp(-2,2,-2,2)))
     &   +1._dp
     &   *(m14_2x(0,pp(-4,2,-4,2))+m14_2x(1,pp(-4,2,-4,2))+m14_2x(2,pp(-4,2,-4,2))))
         endif
        elseif (k < 0) then
          msq( 7,j,k)=msq( 7,j,k)+(real(nf-1,dp)+0.5_dp)*half
     &    *(sub54_1(gq)
     &    *(m14_5x(1,pp(0,kk(k),0,kk(k)))+m14_5x(2,pp(0,kk(k),0,kk(k))))
     &    -sub54_1v
     &    *(m54_1vg(1,0,kk(k))+m54_1vg(2,0,kk(k))))
          msq(18,j,k)=msq(18,j,k)+(real(nf-1,dp)+0.5_dp)*half
     &    *(sub54_2(gq)
     &    *(m24_5x(0,pp(0,kk(k),0,kk(k)))+m24_5x(1,pp(0,kk(k),0,kk(k))))
     &    -sub54_2v
     &    *(m54_2vg(0,0,kk(k))+m54_2vg(1,0,kk(k))))
          msq( 5,j,k)=msq( 5,j,k)+(real(nf-1,dp)+0.5_dp)*half
     &    *(sub54_6(gq)
     &    *(m54_6x(0,pp(0,kk(k),0,kk(k)))+m54_6x(2,pp(0,kk(k),0,kk(k))))
     &    -sub54_6v
     &    *(m54_6vg(0,0,kk(k))+m54_6vg(2,0,kk(k))))

          msq( 9,j,k)=msq( 9,j,k)+0.5_dp*half*sub46_1(gq)
     &    *(m16_4x(1,pp(0,kk(k),0,kk(k)))+m16_4x(2,pp(0,kk(k),0,kk(k))))
     &                           -0.5_dp*half*sub46_1v
     &    *(m46_1vg(1,0,kk(k))+m46_1vg(2,0,kk(k)))
          msq(17,j,k)=msq(17,j,k)+0.5_dp*half*sub46_2(gq)
     &    *(m26_4x(0,pp(0,kk(k),0,kk(k)))+m26_4x(1,pp(0,kk(k),0,kk(k))))
     &                           -0.5_dp*half*sub46_2v
     &    *(m46_2vg(0,0,kk(k))+m46_2vg(1,0,kk(k)))
          msq( 2,j,k)=msq( 2,j,k)+0.5_dp*half*sub46_5(gq)
     &    *(m46_5x(0,pp(0,kk(k),0,kk(k)))+m46_5x(2,pp(0,kk(k),0,kk(k))))
     &                           -0.5_dp*half*sub46_5v
     &    *(m46_5vg(0,0,kk(k))+m46_5vg(2,0,kk(k)))

          if (kk(k) == -1) then
          msq(22,j,k)=msq(22,j,k)+0.5_dp*(aveqg/avegg)*(sub25_1(gq)
     &    *(m25_1x(0,pp(0,0,1,-1))+m25_1x(1,pp(0,0,1,-1))+m25_1x(2,pp(0,0,1,-1)))
     &                     +sub25_1v*m25_1vx(pp(0,0,1,-1)))

          msq(23,j,k)=msq(23,j,k)+(aveqg/avegg)*(
     &    +2.5_dp*(sub26_1(gq)
     &    *(m26_1x(0,pp(0,0,1,-1))+m26_1x(1,pp(0,0,1,-1))+m26_1x(2,pp(0,0,1,-1)))
     &    +sub26_1v*m26_1vx(pp(0,0,1,-1)))
     &    +2._dp*(sub26_1(gq)
     &    *(m26_1x(0,pp(0,0,2,-2))+m26_1x(1,pp(0,0,2,-2))+m26_1x(2,pp(0,0,2,-2)))
     &    +sub26_1v*m26_1vx(pp(0,0,2,-2))))

          msq(19,j,k)=msq(19,j,k)+(xn-one/xn)*sub15_2(qg)*aveqg/aveqq
     &   *(0.5_dp
     &   *(m15_2x(0,pp(1,-1,1,-1))+m15_2x(1,pp(1,-1,1,-1))+m15_2x(2,pp(1,-1,1,-1)))
     &   +2._dp
     &   *(m15_2x(0,pp(2,-1,2,-1))+m15_2x(1,pp(2,-1,2,-1))+m15_2x(2,pp(2,-1,2,-1)))
     &   +2._dp
     &   *(m15_2x(0,pp(3,-1,3,-1))+m15_2x(1,pp(3,-1,3,-1))+m15_2x(2,pp(3,-1,3,-1))))

          msq(20,j,k)=msq(20,j,k)+(xn-one/xn)*sub16_2(qg)*aveqg/aveqq
     &   *(0.5_dp
     &   *(m16_2x(0,pp(1,-1,1,-1))+m16_2x(1,pp(1,-1,1,-1))+m16_2x(2,pp(1,-1,1,-1)))
     &   +2._dp
     &   *(m16_2x(0,pp(1,-1,2,-2))+m16_2x(1,pp(1,-1,2,-2))+m16_2x(2,pp(1,-1,2,-2)))
     &   +2._dp
     &   *(m16_2x(0,pp(1,-1,3,-3))+m16_2x(1,pp(1,-1,3,-3))+m16_2x(2,pp(1,-1,3,-3))))

          msq(21,j,k)=msq(21,j,k)+(xn-one/xn)*sub14_2(qg)*aveqg/aveqq*(
     &     (m14_2x(0,pp(-1,-1,-1,-1))+m14_2x(1,pp(-1,-1,-1,-1))
     &     +m14_2x(2,pp(-1,-1,-1,-1)))
     &    +2._dp*(m14_2x(0,pp(-3,-1,-3,-1))+m14_2x(1,pp(-3,-1,-3,-1))
     &     +m14_2x(2,pp(-3,-1,-3,-1)))
     &    +2._dp*(m14_2x(0,pp(-2,-1,-2,-1))+m14_2x(1,pp(-2,-1,-2,-1))
     &     +m14_2x(2,pp(-2,-1,-2,-1))))

          endif
          if (kk(k) == -2) then
          msq(22,j,k)=msq(22,j,k)+0.5_dp*(aveqg/avegg)*(sub25_1(gq)
     &    *(m25_1x(0,pp(0,0,2,-2))+m25_1x(1,pp(0,0,2,-2))+m25_1x(2,pp(0,0,2,-2)))
     &                     +sub25_1v*m25_1vx(pp(0,0,2,-2)))

          msq(23,j,k)=msq(23,j,k)+(aveqg/avegg)*(
     &    +3._dp*(sub26_1(gq)
     &    *(m26_1x(0,pp(0,0,1,-1))+m26_1x(1,pp(0,0,1,-1))+m26_1x(2,pp(0,0,1,-1)))
     &    +sub26_1v*m26_1vx(pp(0,0,1,-1)))
     &    +1.5_dp*(sub26_1(gq)
     &    *(m26_1x(0,pp(0,0,2,-2))+m26_1x(1,pp(0,0,2,-2))+m26_1x(2,pp(0,0,2,-2)))
     &    +sub26_1v*m26_1vx(pp(0,0,2,-2))))

          msq(19,j,k)=msq(19,j,k)+(xn-one/xn)*sub15_2(qg)*aveqg/aveqq
     &   *(3._dp
     &   *(m15_2x(0,pp(1,-2,1,-2))+m15_2x(1,pp(1,-2,1,-2))+m15_2x(2,pp(1,-2,1,-2)))
     &   +0.5_dp
     &    *(m15_2x(0,pp(2,-2,2,-2))+m15_2x(1,pp(2,-2,2,-2))+m15_2x(2,pp(2,-2,2,-2)))
     &   +1._dp
     &   *(m15_2x(0,pp(4,-2,4,-2))+m15_2x(1,pp(4,-2,4,-2))+m15_2x(2,pp(4,-2,4,-2))))

          msq(20,j,k)=msq(20,j,k)+(xn-one/xn)*sub16_2(qg)*aveqg/aveqq
     &   *(0.5_dp
     &   *(m16_2x(0,pp(2,-2,2,-2))+m16_2x(1,pp(2,-2,2,-2))+m16_2x(2,pp(2,-2,2,-2)))
     &   +3._dp
     &   *(m16_2x(0,pp(2,-2,1,-1))+m16_2x(1,pp(2,-2,1,-1))+m16_2x(2,pp(2,-2,1,-1)))
     &   +1._dp
     &   *(m16_2x(0,pp(2,-2,4,-4))+m16_2x(1,pp(2,-2,4,-4))+m16_2x(2,pp(2,-2,4,-4))))

          msq(21,j,k)=msq(21,j,k)+(xn-one/xn)*sub14_2(qg)*aveqg/aveqq*(
     &     3._dp*(m14_2x(0,pp(-1,-2,-1,-2))+m14_2x(1,pp(-1,-2,-1,-2))
     &     +m14_2x(2,pp(-1,-2,-1,-2)))
     &    +(m14_2x(0,pp(-4,-2,-4,-2))+m14_2x(1,pp(-4,-2,-4,-2))
     &     +m14_2x(2,pp(-4,-2,-4,-2)))
     &    +(m14_2x(0,pp(-2,-2,-2,-2))+m14_2x(1,pp(-2,-2,-2,-2))
     &     +m14_2x(2,pp(-2,-2,-2,-2))))

          endif
        endif
      elseif (k==0) then
c---------Q-G and Qbar-G
        if     (j > 0) then
          msq( 7,j,k)=msq( 7,j,k)+(real(nf-1,dp)+0.5_dp)*half
     &    *(sub54_1(gq)
     &    *(m14_5x(0,pp(jj(j),0,0,jj(j)))+m14_5x(2,pp(jj(j),0,0,jj(j))))
     &    -sub54_1v
     &    *(m54_1vg(0,jj(j),0)+m54_1vg(2,jj(j),0)))
          msq(18,j,k)=msq(18,j,k)+(real(nf-1,dp)+0.5_dp)*half
     &    *(sub54_2(gq)
     &    *(m24_5x(1,pp(jj(j),0,0,jj(j)))+m24_5x(2,pp(jj(j),0,0,jj(j))))
     &    -sub54_2v
     &    *(m54_2vg(1,jj(j),0)+m54_2vg(2,jj(j),0)))
          msq( 5,j,k)=msq( 5,j,k)+(real(nf-1,dp)+0.5_dp)*half
     &    *(sub54_6(gq)
     &    *(m54_6x(0,pp(jj(j),0,0,jj(j)))+m54_6x(1,pp(jj(j),0,0,jj(j))))
     &    -sub54_6v
     &    *(m54_6vg(0,jj(j),0)+m54_6vg(1,jj(j),0)))

          msq(12,j,k)=msq(12,j,k)+0.5_dp*half*sub56_1(gq)
     &    *(m16_5x(0,pp(jj(j),0,jj(j),0))+m16_5x(2,pp(jj(j),0,jj(j),0)))
     &                           -0.5_dp*half*sub56_1v
     &    *(m56_1vg(0,jj(j),0)+m56_1vg(2,jj(j),0))
          msq(13,j,k)=msq(13,j,k)+0.5_dp*half*sub56_2(gq)
     &    *(m26_5x(1,pp(jj(j),0,jj(j),0))+m26_5x(2,pp(jj(j),0,jj(j),0)))
     &                           -0.5_dp*half*sub56_2v
     &    *(m56_2vg(1,jj(j),0)+m56_2vg(2,jj(j),0))
          msq( 4,j,k)=msq( 4,j,k)+0.5_dp*half*sub56_4(gq)
     &    *(m56_4x(0,pp(jj(j),0,jj(j),0))+m56_4x(1,pp(jj(j),0,jj(j),0)))
     &                           -0.5_dp*half*sub56_4v
     &    *(m56_4vg(0,jj(j),0)+m56_4vg(1,jj(j),0))

          if (jj(j) == 1) then
          msq(20,j,k)=msq(20,j,k)+(aveqg/avegg)*(
     &    +2.5_dp*(sub16_2(gq)
     &    *(m16_2x(0,pp(0,0,1,-1))+m16_2x(1,pp(0,0,1,-1))+m16_2x(2,pp(0,0,1,-1)))
     &    +sub16_2v*m16_2vx(pp(0,0,1,-1)))
     &    +2._dp*(sub16_2(gq)
     &    *(m16_2x(0,pp(0,0,2,-2))+m16_2x(1,pp(0,0,2,-2))+m16_2x(2,pp(0,0,2,-2)))
     &    +sub16_2v*m16_2vx(pp(0,0,2,-2))))

          msq(21,j,k)=msq(21,j,k)+0.5_dp*(aveqg/avegg)*(sub14_2(gq)
     &    *(m14_2x(0,pp(0,0,-1,1))+m14_2x(1,pp(0,0,-1,1))+m14_2x(2,pp(0,0,-1,1)))
     &                     +sub14_2v*m14_2vx(pp(0,0,-1,1)))

          msq(22,j,k)=msq(22,j,k)+(xn-one/xn)*sub25_1(qg)*aveqg/aveqq*(
     &     (m25_1x(0,pp(1,1,1,1))+m25_1x(1,pp(1,1,1,1))+m25_1x(2,pp(1,1,1,1)))
     &    +2._dp*(m25_1x(0,pp(1,3,3,1))+m25_1x(1,pp(1,3,3,1))+m25_1x(2,pp(1,3,3,1)))
     &    +2._dp*(m25_1x(0,pp(1,2,2,1))+m25_1x(1,pp(1,2,2,1))+m25_1x(2,pp(1,2,2,1))))
 
          msq(23,j,k)=msq(23,j,k)+(xn-one/xn)*sub26_1(qg)*aveqg/aveqq
     &   *(0.5_dp
     &   *(m26_1x(0,pp(1,-1,1,-1))+m26_1x(1,pp(1,-1,1,-1))+m26_1x(2,pp(1,-1,1,-1)))
     &   +2._dp
     &   *(m26_1x(0,pp(1,-1,2,-2))+m26_1x(1,pp(1,-1,2,-2))+m26_1x(2,pp(1,-1,2,-2)))
     &   +2._dp
     &   *(m26_1x(0,pp(1,-1,3,-3))+m26_1x(1,pp(1,-1,3,-3))+m26_1x(2,pp(1,-1,3,-3))))


          msq(24,j,k)=msq(24,j,k)+(xn-one/xn)*sub24_1(qg)*aveqg/aveqq
     &   *(0.5_dp
     &   *(m24_1x(0,pp(1,-1,-1,1))+m24_1x(1,pp(1,-1,-1,1))+m24_1x(2,pp(1,-1,-1,1)))
     &   +2._dp
     &   *(m24_1x(0,pp(1,-2,-2,1))+m24_1x(1,pp(1,-2,-2,1))+m24_1x(2,pp(1,-2,-2,1)))
     &   +2._dp
     &   *(m24_1x(0,pp(1,-3,-3,1))+m24_1x(1,pp(1,-3,-3,1))+m24_1x(2,pp(1,-3,-3,1))))

           endif
          if (jj(j) == 2) then
          msq(20,j,k)=msq(20,j,k)+(aveqg/avegg)*(
     &    +3._dp*(sub16_2(gq)
     &    *(m16_2x(0,pp(0,0,1,-1))+m16_2x(1,pp(0,0,1,-1))+m16_2x(2,pp(0,0,1,-1)))
     &    +sub16_2v*m16_2vx(pp(0,0,1,-1)))
     &    +1.5_dp*(sub16_2(gq)
     &    *(m16_2x(0,pp(0,0,2,-2))+m16_2x(1,pp(0,0,2,-2))+m16_2x(2,pp(0,0,2,-2)))
     &    +sub16_2v*m16_2vx(pp(0,0,2,-2))))

          msq(21,j,k)=msq(21,j,k)+0.5_dp*(aveqg/avegg)*(sub14_2(gq)
     &    *(m14_2x(0,pp(0,0,-2,2))+m14_2x(1,pp(0,0,-2,2))+m14_2x(2,pp(0,0,-2,2)))
     &                     +sub14_2v*m14_2vx(pp(0,0,-2,2)))

          msq(22,j,k)=msq(22,j,k)+(xn-one/xn)*sub25_1(qg)*aveqg/aveqq*(
     &     3._dp*(m25_1x(0,pp(2,1,1,2))+m25_1x(1,pp(2,1,1,2))+m25_1x(2,pp(2,1,1,2)))
     &    +(m25_1x(0,pp(2,4,4,2))+m25_1x(1,pp(2,4,4,2))+m25_1x(2,pp(2,4,4,2)))
     &    +(m25_1x(0,pp(2,2,2,2))+m25_1x(1,pp(2,2,2,2))+m25_1x(2,pp(2,2,2,2))))
 
          msq(23,j,k)=msq(23,j,k)+(xn-one/xn)*sub26_1(qg)*aveqg/aveqq
     &   *(0.5_dp
     &   *(m26_1x(0,pp(2,-2,2,-2))+m26_1x(1,pp(2,-2,2,-2))+m26_1x(2,pp(2,-2,2,-2)))
     &   +3._dp
     &   *(m26_1x(0,pp(2,-2,1,-1))+m26_1x(1,pp(2,-2,1,-1))+m26_1x(2,pp(2,-2,1,-1)))
     &   +1._dp
     &   *(m26_1x(0,pp(2,-2,4,-4))+m26_1x(1,pp(2,-2,4,-4))+m26_1x(2,pp(2,-2,4,-4))))

          msq(24,j,k)=msq(24,j,k)+(xn-one/xn)*sub24_1(qg)*aveqg/aveqq
     &   *(3._dp
     &   *(m24_1x(0,pp(2,-1,-1,2))+m24_1x(1,pp(2,-1,-1,2))+m24_1x(2,pp(2,-1,-1,2)))
     &   +0.5_dp
     &   *(m24_1x(0,pp(2,-2,-2,2))+m24_1x(1,pp(2,-2,-2,2))+m24_1x(2,pp(2,-2,-2,2)))
     &   +1._dp
     &   *(m24_1x(0,pp(2,-4,-4,2))+m24_1x(1,pp(2,-4,-4,2))+m24_1x(2,pp(2,-4,-4,2))))

          endif
        elseif (j < 0) then
          msq( 7,j,k)=msq( 7,j,k)+(real(nf-1,dp)+0.5_dp)*half
     &    *(sub54_1(gq)
     &    *(m14_5x(0,pp(jj(j),0,0,jj(j)))+m14_5x(1,pp(jj(j),0,0,jj(j))))
     &    -sub54_1v
     &    *(m54_1vg(0,jj(j),0)+m54_1vg(1,jj(j),0)))
          msq(18,j,k)=msq(18,j,k)+(real(nf-1,dp)+0.5_dp)*half
     &    *(sub54_2(gq)
     &    *(m24_5x(1,pp(jj(j),0,0,jj(j)))+m24_5x(2,pp(jj(j),0,0,jj(j))))
     &    -sub54_2v
     &    *(m54_2vg(1,jj(j),0)+m54_2vg(2,jj(j),0)))
          msq( 5,j,k)=msq( 5,j,k)+(real(nf-1,dp)+0.5_dp)*half
     &    *(sub54_6(gq)
     &    *(m54_6x(0,pp(jj(j),0,0,jj(j)))+m54_6x(2,pp(jj(j),0,0,jj(j))))
     &    -sub54_6v
     &    *(m54_6vg(0,jj(j),0)+m54_6vg(2,jj(j),0)))

          msq( 9,j,k)=msq( 9,j,k)+0.5_dp*half*sub46_1(gq)
     &    *(m16_4x(0,pp(jj(j),0,0,jj(j)))+m16_4x(1,pp(jj(j),0,0,jj(j))))
     &                           -0.5_dp*half*sub46_1v
     &    *(m46_1vg(0,jj(j),0)+m46_1vg(1,jj(j),0))
          msq(17,j,k)=msq(17,j,k)+0.5_dp*half*sub46_2(gq)
     &    *(m26_4x(1,pp(jj(j),0,0,jj(j)))+m26_4x(2,pp(jj(j),0,0,jj(j))))
     &                           -0.5_dp*half*sub46_2v
     &    *(m46_2vg(1,jj(j),0)+m46_2vg(2,jj(j),0))
          msq( 2,j,k)=msq( 2,j,k)+0.5_dp*half*sub46_5(gq)
     &    *(m46_5x(0,pp(jj(j),0,0,jj(j)))+m46_5x(2,pp(jj(j),0,0,jj(j))))
     &                           -0.5_dp*half*sub46_5v
     &    *(m46_5vg(0,jj(j),0)+m46_5vg(2,jj(j),0))

          if (jj(j) == -1) then
          msq(19,j,k)=msq(19,j,k)+0.5_dp*(aveqg/avegg)*(sub15_2(gq)
     &    *(m15_2x(0,pp(0,0,1,-1))+m15_2x(1,pp(0,0,1,-1))+m15_2x(2,pp(0,0,1,-1)))
     &                     +sub15_2v*m15_2vx(pp(0,0,1,-1)))

          msq(20,j,k)=msq(20,j,k)+(aveqg/avegg)*(
     &    +2.5_dp*(sub16_2(gq)
     &    *(m16_2x(0,pp(0,0,1,-1))+m16_2x(1,pp(0,0,1,-1))+m16_2x(2,pp(0,0,1,-1)))
     &    +sub16_2v*m16_2vx(pp(0,0,1,-1)))
     &    +2._dp*(sub16_2(gq)
     &    *(m16_2x(0,pp(0,0,2,-2))+m16_2x(1,pp(0,0,2,-2))+m16_2x(2,pp(0,0,2,-2)))
     &    +sub16_2v*m16_2vx(pp(0,0,2,-2))))

          msq(22,j,k)=msq(22,j,k)+(xn-one/xn)*sub25_1(qg)*aveqg/aveqq
     &   *(0.5_dp
     &   *(m25_1x(0,pp(-1,1,1,-1))+m25_1x(1,pp(-1,1,1,-1))+m25_1x(2,pp(-1,1,1,-1)))
     &   +2._dp
     &   *(m25_1x(0,pp(-1,2,2,-1))+m25_1x(1,pp(-1,2,2,-1))+m25_1x(2,pp(-1,2,2,-1)))
     &   +2._dp
     &   *(m25_1x(0,pp(-1,3,3,-1))+m25_1x(1,pp(-1,3,3,-1))+m25_1x(2,pp(-1,3,3,-1))))

          msq(23,j,k)=msq(23,j,k)+(xn-one/xn)*sub26_1(qg)*aveqg/aveqq
     &   *(0.5_dp
     &   *(m26_1x(0,pp(-1,1,1,-1))+m26_1x(1,pp(-1,1,1,-1))+m26_1x(2,pp(-1,1,1,-1)))
     &   +2._dp
     &   *(m26_1x(0,pp(-1,1,2,-2))+m26_1x(1,pp(-1,1,2,-2))+m26_1x(2,pp(-1,1,2,-2)))
     &   +2._dp
     &   *(m26_1x(0,pp(-1,1,3,-3))+m26_1x(1,pp(-1,1,3,-3))+m26_1x(2,pp(-1,1,3,-3))))

          msq(24,j,k)=msq(24,j,k)+(xn-one/xn)*sub24_1(qg)*aveqg/aveqq*(
     &     (m24_1x(0,pp(-1,-1,-1,-1))+m24_1x(1,pp(-1,-1,-1,-1))
     &     +m24_1x(2,pp(-1,-1,-1,-1)))
     &    +2._dp*(m24_1x(0,pp(-1,-3,-3,-1))+m24_1x(1,pp(-1,-3,-3,-1))
     &     +m24_1x(2,pp(-1,-3,-3,-1)))
     &    +2._dp*(m24_1x(0,pp(-1,-2,-2,-1))+m24_1x(1,pp(-1,-2,-2,-1))
     &     +m24_1x(2,pp(-1,-2,-2,-1))))

           endif
          if (jj(j) == -2) then
          msq(19,j,k)=msq(19,j,k)+0.5_dp*(aveqg/avegg)*(sub15_2(gq)
     &    *(m15_2x(0,pp(0,0,2,-2))+m15_2x(1,pp(0,0,2,-2))+m15_2x(2,pp(0,0,2,-2)))
     &                     +sub15_2v*m15_2vx(pp(0,0,2,-2)))

          msq(20,j,k)=msq(20,j,k)+(aveqg/avegg)*(
     &    +3._dp*(sub16_2(gq)
     &    *(m16_2x(0,pp(0,0,1,-1))+m16_2x(1,pp(0,0,1,-1))+m16_2x(2,pp(0,0,1,-1)))
     &    +sub16_2v*m16_2vx(pp(0,0,1,-1)))
     &    +1.5_dp*(sub16_2(gq)
     &    *(m16_2x(0,pp(0,0,2,-2))+m16_2x(1,pp(0,0,2,-2))+m16_2x(2,pp(0,0,2,-2)))
     &    +sub16_2v*m16_2vx(pp(0,0,2,-2))))

          msq(22,j,k)=msq(22,j,k)+(xn-one/xn)*sub25_1(qg)*aveqg/aveqq
     &   *(3._dp
     &   *(m25_1x(0,pp(-2,1,1,-2))+m25_1x(1,pp(-2,1,1,-2))+m25_1x(2,pp(-2,1,1,-2)))
     &   +0.5_dp
     &    *(m25_1x(0,pp(-2,2,2,-2))+m25_1x(1,pp(-2,2,2,-2))+m25_1x(2,pp(-2,2,2,-2)))
     &   +1._dp
     &   *(m25_1x(0,pp(-2,4,4,-2))+m25_1x(1,pp(-2,4,4,-2))+m25_1x(2,pp(-2,4,4,-2))))

          msq(23,j,k)=msq(23,j,k)+(xn-one/xn)*sub26_1(qg)*aveqg/aveqq
     &   *(0.5_dp
     &   *(m26_1x(0,pp(-2,2,2,-2))+m26_1x(1,pp(-2,2,2,-2))+m26_1x(2,pp(-2,2,2,-2)))
     &   +3._dp
     &   *(m26_1x(0,pp(-2,2,1,-1))+m26_1x(1,pp(-2,2,1,-1))+m26_1x(2,pp(-2,2,1,-1)))
     &   +1._dp
     &   *(m26_1x(0,pp(-2,2,4,-4))+m26_1x(1,pp(-2,2,4,-4))+m26_1x(2,pp(-2,2,4,-4))))

          msq(24,j,k)=msq(24,j,k)+(xn-one/xn)*sub24_1(qg)*aveqg/aveqq*(
     &     3._dp*(m24_1x(0,pp(-2,-1,-1,-2))+m24_1x(1,pp(-2,-1,-1,-2))
     &     +m24_1x(2,pp(-2,-1,-1,-2)))
     &    +(m24_1x(0,pp(-2,-4,-4,-2))+m24_1x(1,pp(-2,-4,-4,-2))
     &     +m24_1x(2,pp(-2,-4,-4,-2)))
     &    +(m24_1x(0,pp(-2,-2,-2,-2))+m24_1x(1,pp(-2,-2,-2,-2))
     &     +m24_1x(2,pp(-2,-2,-2,-2))))

          endif
        endif
      endif


      enddo
      enddo
      
      endif 

      return
      end

