!---  CW June 16
!-----New version of tree-level process for gamma + 2j, written in the style of
!---- Z2jet
!---- this version returns the various color basis MEs
      subroutine qqb_gam2jx_new(p,msq,mqq,ppmsqx,msqx_cs)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'zcouple.f'
      include 'ewcharge.f'
      include 'nflav.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      include 'flags.f'
      include 'lc.f'
      include 'kpart.f'
      include 'pp.f'
      real(dp) :: p(mxpart,4),msq(-nf:nf,-nf:nf)
      real(dp) :: qqb_gagg,qqb_gagg_cs(0:2)
      real(dp) :: qbq_gagg,qbq_gagg_cs(0:2)
      real(dp) :: qg_gaqg,qg_gaqg_cs(0:2),qg_gagq,qg_gagq_cs(0:2)
      real(dp) :: qbg_gaqbg,qbg_gaqbg_cs(0:2),qbg_gagqb,qbg_gagqb_cs(0:2)
      real(dp) :: gq_gaqg,gq_gaqg_cs(0:2),gq_gagq,gq_gagq_cs(0:2)
      real(dp) :: gqb_gaqbg,gqb_gaqbg_cs(0:2),gqb_gagqb,gqb_gagqb_cs(0:2)
      real(dp) :: gg_gaqqb,gg_gaqqb_cs(0:2),gg_gaqbq,gg_gaqbq_cs(0:2)
      real(dp):: tup,tdo,fac,faclo
    
      complex(dp):: qRb_a(2,2,2),qRb_b(2,2,2)
      complex(dp):: qqb_a(2,2,2),qqb_b(2,2,2)

      complex(dp):: qbq_a(2,2,2),qbq_b(2,2,2)
      complex(dp):: qbR_a(2,2,2),qbR_b(2,2,2)

      complex(dp):: qq_a(2,2,2),qq_b(2,2,2)
      complex(dp):: qR_a(2,2,2),qR_b(2,2,2)

      complex(dp):: qbRb_a(2,2,2),qbRb_b(2,2,2)
      complex(dp):: qbqb_a(2,2,2),qbqb_b(2,2,2)
      
      complex(dp):: qRb_ax(2,2,2),qRb_bx(2,2,2),qqb_ax(2,2,2),qqb_bx(2,2,2)
      complex(dp):: qbR_ax(2,2,2),qbR_bx(2,2,2),qbq_ax(2,2,2),qbq_bx(2,2,2)
      
      real(dp):: mqq(0:2,fn:nf,fn:nf),msq0,msq1,msq2
      real(dp):: ppmsqx(0:2,ppmax),tmpmsqx(0:2,0:ppmax)
      real(dp):: msqx_cs(0:2,-nf:nf,-nf:nf)

      integer:: rcolourchoice
      logical:: rGflag
      integer,parameter::swap(2)=(/2,1/),swap1(0:2)=(/0,2,1/)
     
      integer:: i,j,k,f,j1,j2,j3,nquark,nup,ndo,icol
      real(dp):: Qsum_all

c--- if we're calculating the REAL or VIRT matrix elements, we
c--- need all the colour structures, but want to preserve
c--- the actual value of colourchoice
      if ((kpart==kreal) .or. (kpart==kvirt)) then
        rcolourchoice=colourchoice
        colourchoice=0
      endif     
c--- if we're calculating the REAL matrix elements with Qflag=TRUE,
c    the subtraction terms involve the (Gflag=TRUE) matrix elements
      if ((kpart==kreal) .and. (Qflag .eqv. .true.)) then
        rGflag=Gflag
        Gflag=.true.
      endif

      msq(:,:)=zip

      Qsum_all=zip
      do j=1,nf
         Qsum_all=Qsum_all+Q(j)**2
      enddo
!-----DEBUG FOR KC
!      include 'kinpoint_qqbggga_kc.f'
!      include 'kinpoint_qqbQQga_kc.f'
!      call spinorz(5,p,za,zb)
!      call ga_qqbQQb(1,2,3,4,za,zb,qqb_a,qqb_b)


!     call qqbgg_ga(1,2,3,4,5,za,zb,qqb_gagg)
!     ---- END DEBUG
      
      call spinoru(5,p,za,zb)

!=====CW note I had to swap 4<->5 to get pole cancellation,
!===== I havent updated the names yet so this is a bit misleading
      if(Gflag) then
!-------calculate 2-quark, 2-gluon amplitudes
      call qqbgg_ga(2,1,4,5,3,za,zb,qqb_gagg)
      call storecs_ga(qqb_gagg_cs)
      call qqbgg_ga(4,1,2,5,3,za,zb,qg_gaqg)
      call storecs_ga(qg_gaqg_cs)
      call qqbgg_ga(4,2,1,5,3,za,zb,gq_gaqg)
      call storecs_ga(gq_gaqg_cs)

      call qqbgg_ga(5,1,2,4,3,za,zb,qg_gagq)
      call storecs_ga(qg_gagq_cs)
      call qqbgg_ga(5,2,1,4,3,za,zb,gq_gagq)
      call storecs_ga(gq_gagq_cs)

      call qqbgg_ga(5,4,1,2,3,za,zb,gg_gaqbq)
      call storecs_ga(gg_gaqbq_cs)

      do i=0,2
       qbq_gagg_cs(i)=qqb_gagg_cs(swap1(i))
       qbg_gaqbg_cs(i)=qg_gaqg_cs(swap1(i))
       gqb_gaqbg_cs(i)=gq_gaqg_cs(swap1(i))
       qbg_gagqb_cs(i)=qg_gagq_cs(swap1(i))
       gqb_gagqb_cs(i)=gq_gagq_cs(swap1(i))
       gg_gaqqb_cs(i)=gg_gaqbq_cs(swap1(i))
      enddo
      
      fac=2._dp*esq*gsq**2*V*xn

      do i=0,2
         qqb_gagg_cs(i)=half*aveqq*fac*qqb_gagg_cs(i)
         qbq_gagg_cs(i)=half*aveqq*fac*qbq_gagg_cs(i)
         gq_gaqg_cs(i)=aveqg*fac*gq_gaqg_cs(i)
         qg_gaqg_cs(i)=aveqg*fac*qg_gaqg_cs(i)
         gqb_gaqbg_cs(i)=aveqg*fac*gqb_gaqbg_cs(i)
         qbg_gaqbg_cs(i)=aveqg*fac*qbg_gaqbg_cs(i)
         gq_gagq_cs(i)=aveqg*fac*gq_gagq_cs(i)
         qg_gagq_cs(i)=aveqg*fac*qg_gagq_cs(i)
         gqb_gagqb_cs(i)=aveqg*fac*gqb_gagqb_cs(i)
         qbg_gagqb_cs(i)=aveqg*fac*qbg_gagqb_cs(i)
         gg_gaqqb_cs(i)=avegg*fac*gg_gaqqb_cs(i)
         gg_gaqbq_cs(i)=avegg*fac*gg_gaqbq_cs(i)
      enddo

      
      qqb_gagg=qqb_gagg_cs(0)+qqb_gagg_cs(1)+qqb_gagg_cs(2)
      qbq_gagg=qqb_gagg_cs(0)+qbq_gagg_cs(1)+qbq_gagg_cs(2)
      gq_gaqg=gq_gaqg_cs(0)+gq_gaqg_cs(1)+gq_gaqg_cs(2)
      qg_gaqg=qg_gaqg_cs(0)+qg_gaqg_cs(1)+qg_gaqg_cs(2)
      gqb_gaqbg=gqb_gaqbg_cs(0)+gqb_gaqbg_cs(1)+gqb_gaqbg_cs(2)
      qbg_gaqbg=qbg_gaqbg_cs(0)+qbg_gaqbg_cs(1)+qbg_gaqbg_cs(2)
      gg_gaqqb=gg_gaqqb_cs(0)+gg_gaqqb_cs(1)+gg_gaqqb_cs(2)
      gg_gaqbq=gg_gaqbq_cs(0)+gg_gaqbq_cs(1)+gg_gaqbq_cs(2)


      endif
!---- i1-> i1
!-----i2-> i4
!-----i4--> i2
!---- i5--> i5

!-----four-quark amplitudes
      if (Qflag) then 

!----- ga_qqbQQB(u(1),ub(2),d(4),db(5),ga(3))

!      q Qb -> q Qb
        call ga_qqbQQb(1,4,5,2,za,zb,qRb_a,qRb_b)
!      q Q -> q Q
! instead of calling call ga_qqbQQb(1,4,2,5,za,zb,qR_a,qR_b)
        do j2=1,2
        do j3=1,2
        qR_a(:,j2,j3)=+qRb_a(:,swap(j2),j3)
        qR_b(:,j2,j3)=-qRb_b(:,swap(j2),j3)
        enddo
        enddo

!------ qb Q -> qb Q
        call ga_qqbQQb(5,1,2,4,za,zb,qbR_a,qbR_b)
!------ qb Qb -> qb Qb
! instead of calling call ga_qqbQQb(4,1,5,2,za,zb,qbRb_a,qbRb_b)
        do j1=1,2
        qbRb_a(j1,:,:)=-qRb_a(swap(j1),:,:)
        qbRb_b(j1,:,:)=+qRb_b(swap(j1),:,:)
        enddo

!-----q qb -> QQb
        call ga_qqbQQb(1,2,5,4,za,zb,qqb_a,qqb_b)
!-----qb q -> QQb
! instead of calling call ga_qqbQQb(2,1,5,4,za,zb,qbq_a,qbq_b)
        do j1=1,2
        qbq_a(j1,:,:)=-qqb_a(swap(j1),:,:)
        qbq_b(j1,:,:)=+qqb_b(swap(j1),:,:)
        enddo

!-----qq->qq
        call ga_qqbQQb(1,5,2,4,za,zb,qq_a,qq_b)
!-----qbqb->qbqb
! instead of calling ga_qqbQQb(5,1,4,2,za,zb,qbqb_a,qbqb_b)
        do j1=1,2
        do j2=1,2
        qbqb_a(j1,j2,:)=-qq_a(swap(j1),swap(j2),:)
        qbqb_b(j1,j2,:)=-qq_b(swap(j1),swap(j2),:)
        enddo
        enddo

c--- four extra amplitudes added for 4-label matrix elements

c--- qRb->Rbq
c instead of calling ga_qqbQQb(1,5,4,2,za,zb,qRb_ax,qRb_bx)
      do j1=1,2
        do j2=1,2
        do j3=1,2
        qRb_ax(j1,j2,j3)=-qbqb_a(swap(j1),j2,j3)
        qRb_bx(j1,j2,j3)=+qbqb_b(swap(j1),j2,j3)
        enddo
        enddo
        enddo
        
c--- qqb->qbq
c instead of calling ga_qqbQQb(1,2,4,5,za,zb,qqb_ax,qqb_bx)
        do j1=1,2
        do j2=1,2
        do j3=1,2
        qqb_ax(j1,j2,j3)=+qqb_a(j1,swap(j2),j3)
        qqb_bx(j1,j2,j3)=-qqb_b(j1,swap(j2),j3)
        enddo
        enddo
        enddo

c--- qbR->qbR    
c instead of calling ga_qqbQQb(4,1,2,5,za,zb,qbR_ax,qbR_bx)
        do j1=1,2
        do j2=1,2
        do j3=1,2
        qbR_ax(j1,j2,j3)=-qRb_a(swap(j1),swap(j2),j3)
        qbR_bx(j1,j2,j3)=-qRb_b(swap(j1),swap(j2),j3)
        enddo
        enddo
        enddo

c--- qbq->qbq    
c instead of calling ga_qqbQQb(2,1,4,5,za,zb,qbq_ax,qbq_bx)
        do j1=1,2
        do j2=1,2
        do j3=1,2
        qbq_ax(j1,j2,j3)=+qbq_a(j1,swap(j2),j3)
        qbq_bx(j1,j2,j3)=-qbq_b(j1,swap(j2),j3)
        enddo
        enddo
        enddo

        faclo=2._dp*V*gsq**2*esq*aveqq

      endif

      
c--- initialize temporary error to zero; note that second dimension
c--- runs from 0 to ppmax to catch irrelevant entries
      tmpmsqx(:,:)=zip

      if (Gflag) then
      do j=-nf,nf
      do k=-nf,nf
      do icol=0,2
        msqx_cs(icol,j,k)=0._dp
      enddo
      
      if( j .ne. 0 .and. k .ne. 0 .and. j .ne. -k) goto 19

      if     ((j == 0) .and. (k == 0)) then

         do icol=0,2
            do nquark=1,2
               tmpmsqx(icol,pp(j,k,-nquark,nquark))=
     &              +abs(Q(nquark))**2*gg_gaqbq_cs(icol)
               tmpmsqx(icol,pp(j,k,nquark,-nquark))=
     &              +abs(Q(nquark))**2*gg_gaqqb_cs(icol)
            enddo
           msqx_cs(icol,j,k)=+3._dp*tmpmsqx(icol,pp(j,k,-1,1))
     &                       +2._dp*tmpmsqx(icol,pp(j,k,-2,2))
! The following commented-out code is required for pole cancellation
! in the virtual contribution with Gflag=.true.
!           msqx_cs(icol,j,k)=+3._dp*tmpmsqx(icol,pp(j,k,1,-1))
!     &                       +2._dp*tmpmsqx(icol,pp(j,k,2,-2))
        enddo
      elseif ((j > 0) .and. (k < 0)) then
          do icol=0,2
             msqx_cs(icol,j,k)=
     &       +abs(Q(j))**2*qqb_gagg_cs(icol)
          enddo

c---Statistical factor already included above
      elseif ((j < 0) .and. (k > 0)) then
          do icol=0,2
             msqx_cs(icol,j,k)=
     &       +abs(Q(k))**2*qbq_gagg_cs(icol)
          enddo
      elseif ((j > 0) .and. (k == 0)) then
          do icol=0,2
c--- normal case
             msqx_cs(icol,j,k)=
     &            +abs(Q(j))**2*qg_gaqg_cs(icol)
             if (j <= 4) then
             tmpmsqx(icol,pp(j,k,j,k))=msqx_cs(icol,j,k)
             tmpmsqx(icol,pp(j,k,k,j))=
     &            +abs(Q(j))**2*qg_gagq_cs(icol)
          endif
          enddo
      elseif ((j < 0) .and. (k == 0)) then
          do icol=0,2
             msqx_cs(icol,j,k)=
     &       +abs(Q(-j))**2*qbg_gaqbg_cs(icol)
             if (j >= -4) then
             tmpmsqx(icol,pp(j,k,j,k))=msqx_cs(icol,j,k)
             tmpmsqx(icol,pp(j,k,k,j))=
     &       +abs(Q(-j))**2*qbg_gagqb_cs(icol)
             endif
          enddo
      elseif ((j == 0) .and. (k > 0)) then
          do icol=0,2
             msqx_cs(icol,j,k)=
     &       +abs(Q(k))**2*gq_gaqg_cs(icol)
             if (k <= 4) then
             tmpmsqx(icol,pp(j,k,k,j))=msqx_cs(icol,j,k)
             tmpmsqx(icol,pp(j,k,j,k))=
     &       +abs(Q(k))**2*gq_gagq_cs(icol)
             endif
          enddo
      elseif ((j == 0) .and. (k < 0)) then
          do icol=0,2
             msqx_cs(icol,j,k)=
     &       +abs(Q(-k))**2*gqb_gaqbg_cs(icol)
             if (k >= -4) then
             tmpmsqx(icol,pp(j,k,k,j))=msqx_cs(icol,j,k)
             tmpmsqx(icol,pp(j,k,j,k))=
     &       +abs(Q(-k))**2*gqb_gagqb_cs(icol)
             endif
          enddo
      endif
      msq(j,k)=msqx_cs(0,j,k)+msqx_cs(1,j,k)+msqx_cs(2,j,k)
   19 continue
      enddo
      enddo
      endif

!---- assemble qq pieces

      mqq(:,:,:)=zip
      
      if (Qflag) then
      
        do j=-nf,nf
        do k=-nf,nf
      
        mqq(:,j,k)=zip
         
        if ((j > 0) .and. (k > 0)) then
c----QQ case
          call msq_gam2jqq(j,k,qR_a,qR_b,qq_a,qq_b,msq0,msq1,msq2)
                      
          mqq(0,j,k)=faclo*msq0
          mqq(1,j,k)=faclo*msq1
          mqq(2,j,k)=faclo*msq2
          
          if ((j <= 4) .and. (k <= 4)) then
            tmpmsqx(:,pp(j,k,j,k))=mqq(:,j,k)
            if (j .ne. k) then
              call msq_gam2jqq(j,k,qq_a,qq_b,qR_a,qR_b,msq0,msq1,msq2)
              tmpmsqx(0,pp(j,k,k,j))=faclo*msq0
              tmpmsqx(1,pp(j,k,k,j))=faclo*msq1
              tmpmsqx(2,pp(j,k,k,j))=faclo*msq2
            endif
          endif

          elseif ((j < 0) .and. (k < 0)) then
c----QbQb case
          call msq_gam2jqq(-j,-k,qbRb_a,qbRb_b,qbqb_a,qbqb_b,
     &                msq0,msq1,msq2)
                      
          mqq(0,j,k)=faclo*msq0
          mqq(1,j,k)=faclo*msq1
          mqq(2,j,k)=faclo*msq2
            
          if ((j >= -4) .and. (k >= -4)) then
            tmpmsqx(:,pp(j,k,j,k))=mqq(:,j,k)
            if (j .ne. k) then
              call msq_gam2jqq(-j,-k,qbqb_a,qbqb_b,qbRb_a,qbRb_b,
     &                    msq0,msq1,msq2)
              tmpmsqx(0,pp(j,k,k,j))=faclo*msq0
              tmpmsqx(1,pp(j,k,k,j))=faclo*msq1
              tmpmsqx(2,pp(j,k,k,j))=faclo*msq2
            endif
          endif
            
        elseif ((j > 0) .and. (k < 0)) then
C---Q-Qb case
          call msq_gam2jqqb(j,k,qRb_a,qRb_b,qqb_a,qqb_b,
     &                msq0,msq1,msq2,tup,tdo)
                      
          mqq(0,j,k)=faclo*msq0
          mqq(1,j,k)=faclo*msq1
          mqq(2,j,k)=faclo*msq2
          
          if ((j <= 4) .and. (k >= -4)) then
            tmpmsqx(:,pp(j,k,j,k))=mqq(:,j,k)
            if (j == -k) then
              do f=1,3,2
                if (f .ne. j) then
                  tmpmsqx(0,pp(j,k,f,-f))=zip
                  tmpmsqx(1,pp(j,k,f,-f))=faclo*tdo
                  tmpmsqx(2,pp(j,k,f,-f))=zip 
                endif
              enddo
              do f=2,4,2
                if (f .ne. j) then
                  tmpmsqx(0,pp(j,k,f,-f))=zip
                  tmpmsqx(1,pp(j,k,f,-f))=faclo*tup
                  tmpmsqx(2,pp(j,k,f,-f))=zip
                endif
              enddo
            endif
          endif

          if (j == -k) then  
            if ((j==1).or.(j==3).or.(j==5)) then
              nup=2
              ndo=nf-3
            else
              nup=1
              ndo=nf-2
            endif
            mqq(1,j,k)=mqq(1,j,k)+faclo
     &      *(real(nup,dp)*tup+real(ndo,dp)*tdo)
          endif
          
          if ((j <= 4) .and. (k >= -4)) then

            call msq_gam2jqqb(j,k,qRb_ax,qRb_bx,qqb_ax,qqb_bx,
     &                msq0,msq1,msq2,tup,tdo)
                    
            tmpmsqx(0,pp(j,k,k,j))=faclo*msq0
            tmpmsqx(1,pp(j,k,k,j))=faclo*msq1
            tmpmsqx(2,pp(j,k,k,j))=faclo*msq2
          
            if (j == -k) then
              do f=1,3,2
                if (f .ne. j) then
                  tmpmsqx(0,pp(j,k,-f,f))=zip
                  tmpmsqx(1,pp(j,k,-f,f))=faclo*tdo
                  tmpmsqx(2,pp(j,k,-f,f))=zip
                endif
              enddo
              do f=2,4,2
                if (f .ne. j) then
                  tmpmsqx(0,pp(j,k,-f,f))=zip
                  tmpmsqx(1,pp(j,k,-f,f))=faclo*tup
                  tmpmsqx(2,pp(j,k,-f,f))=zip
                endif
              enddo
            endif
          endif
 
        elseif ((j < 0) .and. (k > 0)) then
C---Qb-Q case
          call msq_gam2jqqb(-j,-k,qbR_a,qbR_b,qbq_a,qbq_b,
     &                msq0,msq1,msq2,tup,tdo)
                      
          mqq(0,j,k)=faclo*msq0
          mqq(1,j,k)=faclo*msq1
          mqq(2,j,k)=faclo*msq2

          if ((j >= -4) .and. (k <= 4)) then
            tmpmsqx(:,pp(j,k,k,j))=mqq(:,j,k)
            if (j == -k) then
              do f=1,3,2
                if (f .ne. k) then
                  tmpmsqx(0,pp(j,k,f,-f))=zip
                  tmpmsqx(1,pp(j,k,f,-f))=faclo*tdo
                  tmpmsqx(2,pp(j,k,f,-f))=zip 
                endif
              enddo
              do f=2,4,2
                if (f .ne. k) then
                  tmpmsqx(0,pp(j,k,f,-f))=zip
                  tmpmsqx(1,pp(j,k,f,-f))=faclo*tup
                  tmpmsqx(2,pp(j,k,f,-f))=zip
                endif
              enddo
            endif
          endif

          if (j == -k) then 
            if ((k==1).or.(k==3).or.(k==5)) then
              nup=2
              ndo=nf-3
            else
              nup=1
              ndo=nf-2
            endif
            mqq(1,j,k)=mqq(1,j,k)
     &      +faclo*(real(nup,dp)*tup+real(ndo,dp)*tdo)
          endif
          
          if ((j >= -4) .and. (k <= 4)) then

           call msq_gam2jqqb(-j,-k,qbR_ax,qbR_bx,qbq_ax,qbq_bx,
     &               msq0,msq1,msq2,tup,tdo)
                     
           tmpmsqx(0,pp(j,k,j,k))=faclo*msq0
           tmpmsqx(1,pp(j,k,j,k))=faclo*msq1
           tmpmsqx(2,pp(j,k,j,k))=faclo*msq2
           
           
            if (j == -k) then
              do f=1,3,2
                if (f .ne. k) then
                  tmpmsqx(0,pp(j,k,-f,f))=zip
                  tmpmsqx(1,pp(j,k,-f,f))=faclo*tdo
                  tmpmsqx(2,pp(j,k,-f,f))=zip
                endif
              enddo
              do f=2,4,2
                if (f .ne. k) then
                  tmpmsqx(0,pp(j,k,-f,f))=zip
                  tmpmsqx(1,pp(j,k,-f,f))=faclo*tup
                  tmpmsqx(2,pp(j,k,-f,f))=zip
                endif
              enddo
            endif
          endif
      endif

      msq(j,k)=msq(j,k)+mqq(0,j,k)+mqq(1,j,k)+mqq(2,j,k)

      enddo
      enddo

      endif
      
 999  continue
c--- restore proper colourchoice if necessary
      if ((kpart==kreal) .or. (kpart==kvirt)) then
        colourchoice=rcolourchoice
      endif
c--- restore proper parton sub-process selection, if necessary
      if ((kpart==kreal) .and. (Qflag .eqv. .true.)) then
        Gflag=rGflag
      endif

c--- fill ppmsqx array to be returned
      ppmsqx(:,1:ppmax)=tmpmsqx(:,1:ppmax)

      return
      end
