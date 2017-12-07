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
      real(dp) :: qg_gaqg,qg_gaqg_cs(0:2)
      real(dp) :: qbg_gaqbg,qbg_gaqbg_cs(0:2)
      real(dp) :: gq_gaqg,gq_gaqg_cs(0:2)
      real(dp) :: gqb_gaqbg,gqb_gaqbg_cs(0:2)
      real(dp) :: gg_gaqqb, gg_gaqqb_cs(0:2)
      real(dp):: tup,tdo,fac,faclo
    
      complex(dp):: a111,a112,a121,a211,a122,a212,a221,a222
      complex(dp):: b111,b112,b121,b211,b122,b212,b221,b222

      complex(dp):: qRb_a(2,2,2),qRb_b(2,2,2)
      complex(dp):: qqb_a(2,2,2),qqb_b(2,2,2),prop

      complex(dp):: qbq_a(2,2,2),qbq_b(2,2,2)
      complex(dp):: qbR_a(2,2,2),qbR_b(2,2,2)

      complex(dp):: qq_a(2,2,2),qq_b(2,2,2)
      complex(dp):: qR_a(2,2,2),qR_b(2,2,2)

      complex(dp):: qbRb_a(2,2,2),qbRb_b(2,2,2)
      complex(dp):: qbqb_a(2,2,2),qbqb_b(2,2,2)
      
      real(dp):: coupqe(nf,2,2)
      real(dp):: mqq(0:2,fn:nf,fn:nf),msq0,msq1,msq2
      real(dp):: ppmsqx(0:2,ppmax),tmpmsqx(0:2,0:ppmax)
      real(dp):: msqx_cs(0:2,-nf:nf,-nf:nf)

      integer:: rcolourchoice
      logical:: rGflag
      integer,parameter::swap(2)=(/2,1/),swap1(0:2)=(/0,2,1/)
     
      integer:: i,j,k,pq,pl,nquark,nup,ndo,j1,j2,j3,icol
      real(dp):: Qsum_all,Qsum_noj
      real(dp):: test,Bigagam


      
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
!     _---- END DEBUG
      
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

      do i=0,2
       qbq_gagg_cs(i)=qqb_gagg_cs(swap1(i))
       qbg_gaqbg_cs(i)=qg_gaqg_cs(swap1(i))
       gqb_gaqbg_cs(i)=gq_gaqg_cs(swap1(i))
      enddo
      
      call qqbgg_ga(5,4,1,2,3,za,zb,gg_gaqqb)
      call storecs_ga(gg_gaqqb_cs)

      fac=2._dp*esq*gsq**2*V*xn

      do i=0,2
         qqb_gagg_cs(i)=half*aveqq*fac*qqb_gagg_cs(i)
         qbq_gagg_cs(i)=half*aveqq*fac*qbq_gagg_cs(i)
         gq_gaqg_cs(i)=aveqg*fac*gq_gaqg_cs(i)
         qg_gaqg_cs(i)=aveqg*fac*qg_gaqg_cs(i)
         gqb_gaqbg_cs(i)=aveqg*fac*gqb_gaqbg_cs(i)
         qbg_gaqbg_cs(i)=aveqg*fac*qbg_gaqbg_cs(i)
         gg_gaqqb_cs(i)=avegg*fac*gg_gaqqb_cs(i)
      enddo

      
      qqb_gagg=qqb_gagg_cs(0)+qqb_gagg_cs(1)+qqb_gagg_cs(2)
      qbq_gagg=qqb_gagg_cs(0)+qbq_gagg_cs(1)+qbq_gagg_cs(2)
      gq_gaqg=gq_gaqg_cs(0)+gq_gaqg_cs(1)+gq_gaqg_cs(2)
      qg_gaqg=qg_gaqg_cs(0)+qg_gaqg_cs(1)+qg_gaqg_cs(2)
      gqb_gaqbg=gqb_gaqbg_cs(0)+gqb_gaqbg_cs(1)+gqb_gaqbg_cs(2)
      qbg_gaqbg=qbg_gaqbg_cs(0)+qbg_gaqbg_cs(1)+qbg_gaqbg_cs(2)
      gg_gaqqb=gg_gaqqb_cs(0)+gg_gaqqb_cs(1)+gg_gaqqb_cs(2)


      endif
!---- i1-> i1
!-----i2-> i4
!-----i4--> i2
!---- i5--> i5

      if(Qflag) then 
!-----four-quark amplitudes
!----- ga_qqbQQB(u(1),ub(2),d(4),db(5),ga(3))
!      q Q -> q Q
       call ga_qqbQQb(1,4,2,5,za,zb,qR_a,qR_b)

!------ qb Q -> qb Q
       call ga_qqbQQb(5,1,2,4,za,zb,qbR_a,qbR_b)
!-----qb q -> QQb
      call ga_qqbQQb(2,1,5,4,za,zb,qbq_a,qbq_b)

!     q Qb -> q Qb
       call ga_qqbQQb(1,5,4,2,za,zb,qRb_a,qRb_b)

!-----q qb -> QQb
      call ga_qqbQQb(1,2,4,5,za,zb,qqb_a,qqb_b)


      
!------ qb Qb -> qb Qb
       call ga_qqbQQb(4,1,5,2,za,zb,qbRb_a,qbRb_b)      

!-----qq->qq
      call ga_qqbQQb(1,5,2,4,za,zb,qq_a,qq_b)
!-----qbqb->qbqb
      call ga_qqbQQb(5,1,4,2,za,zb,qbqb_a,qbqb_b)

      faclo=2._dp*V*gsq**2*esq*aveqq

      endif

      
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
     &              +abs(Q(nquark))**2*gg_gaqqb_cs(icol)
               tmpmsqx(icol,pp(j,k,nquark,-nquark))=
     &              +abs(Q(nquark))**2*gg_gaqqb_cs(icol)
            enddo
           msqx_cs(icol,j,k)=+3._dp*tmpmsqx(icol,pp(j,k,-1,1))
     &                       +2._dp*tmpmsqx(icol,pp(j,k,-2,2))
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
     &            +abs(Q(j))**2*qg_gaqg_cs(icol)
          endif
          enddo
      elseif ((j < 0) .and. (k == 0)) then
          do icol=0,2
             msqx_cs(icol,j,k)=
     &       +abs(Q(-j))**2*qbg_gaqbg_cs(icol)
             if (j >= -4) then
             tmpmsqx(icol,pp(j,k,j,k))=msqx_cs(icol,j,k)
             tmpmsqx(icol,pp(j,k,k,j))=
     &       +abs(Q(-j))**2*qbg_gaqbg_cs(icol)
             endif
          enddo
      elseif ((j == 0) .and. (k > 0)) then
          do icol=0,2
             msqx_cs(icol,j,k)=
     &       +abs(Q(k))**2*gq_gaqg_cs(icol)
             if (k <= 4) then
             tmpmsqx(icol,pp(j,k,k,j))=msqx_cs(icol,j,k)
             tmpmsqx(icol,pp(j,k,j,k))=
     &       +abs(Q(k))**2*gq_gaqg_cs(icol)
             endif
          enddo
      elseif ((j == 0) .and. (k < 0)) then
          do icol=0,2
             msqx_cs(icol,j,k)=
     &       +abs(Q(-k))**2*gqb_gaqbg_cs(icol)
             if (k >= -4) then
             tmpmsqx(icol,pp(j,k,k,j))=msqx_cs(icol,j,k)
             tmpmsqx(icol,pp(j,k,j,k))=
     &       +abs(Q(-k))**2*gqb_gaqbg_cs(icol)
             endif
          enddo
      endif
      msq(j,k)=msqx_cs(0,j,k)+msqx_cs(1,j,k)+msqx_cs(2,j,k)
   19 continue
      enddo
      enddo
      endif

!---- assemble qq pieces MAKE MORE LIKE Z2JET FOR INTF. 

      mqq(:,:,:)=zip
      if(Qflag) then
         do j=-nf,nf
            do k=-nf,nf
      
               do icol=0,2
                  mqq(icol,j,k)=zip
               enddo
      
          if ((j > 0) .and. (k > 0)) then
c----QQ case
            if (j .ne. k) then
            a111=(Q(j))*qR_a(1,1,1)
     &          +(Q(k))*qR_b(1,1,1)
            a121=(Q(j))*qR_a(1,2,1)
     &          +(Q(k))*qR_b(1,2,1)
            a112=(Q(j))*qR_a(1,1,2)
     &          +(Q(k))*qR_b(1,1,2)
            a122=(Q(j))*qR_a(1,2,2)
     &          +(Q(k))*qR_b(1,2,2)
            a211=(Q(j))*qR_a(2,1,1)
     &          +(Q(k))*qR_b(2,1,1)
            a221=(Q(j))*qR_a(2,2,1)
     &          +(Q(k))*qR_b(2,2,1)
            a212=(Q(j))*qR_a(2,1,2)
     &          +(Q(k))*qR_b(2,1,2)
            a222=(Q(j))*qR_a(2,2,2)
     &          +(Q(k))*qR_b(2,2,2)
            mqq(0,j,k)=zip
            mqq(1,j,k)=faclo*(
     &        real(a111*conjg(a111),dp)+real(a112*conjg(a112),dp)
     &       +real(a221*conjg(a221),dp)+real(a222*conjg(a222),dp)
     &       +real(a122*conjg(a122),dp)+real(a212*conjg(a212),dp)
     &       +real(a121*conjg(a121),dp)+real(a211*conjg(a211),dp))
            mqq(2,j,k)=zip
            elseif (j == k) then
            a111=(Q(j))*(qR_a(1,1,1)+qR_b(1,1,1))
            b111=(Q(j))*(qq_a(1,1,1)+qq_b(1,1,1))
            a112=(Q(j))*(qR_a(1,1,2)+qR_b(1,1,2))
            b112=(Q(j))*(qq_a(1,1,2)+qq_b(1,1,2))
            a221=(Q(j))*(qR_a(2,2,1)+qR_b(2,2,1))
            b221=(Q(j))*(qq_a(2,2,1)+qq_b(2,2,1))
            a222=(Q(j))*(qR_a(2,2,2)+qR_b(2,2,2))
            b222=(Q(j))*(qq_a(2,2,2)+qq_b(2,2,2))

            a121=(Q(j))*qR_a(1,2,1)
     &          +(Q(k))*qR_b(1,2,1)
            b121=(Q(j))*qq_a(1,2,1)
     &          +(Q(k))*qq_b(1,2,1)
            a122=(Q(j))*qR_a(1,2,2)
     &          +(Q(k))*qR_b(1,2,2)
            b122=(Q(j))*qq_a(1,2,2)
     &          +(Q(k))*qq_b(1,2,2)
            a211=(Q(j))*qR_a(2,1,1)
     &          +(Q(k))*qR_b(2,1,1)
            b211=(Q(j))*qq_a(2,1,1)
     &          +(Q(k))*qq_b(2,1,1)
            a212=(Q(j))*qR_a(2,1,2)
     &          +(Q(k))*qR_b(2,1,2)
            b212=(Q(j))*qq_a(2,1,2)
     &          +(Q(k))*qq_b(2,1,2)

            mqq(0,j,k)=half*faclo*(
     &      +real(a111*conjg(b111),dp)+real(a112*conjg(b112),dp)
     &      +real(a221*conjg(b221),dp)+real(a222*conjg(b222),dp))*two/xn
            mqq(1,j,k)=half*faclo*(
     &       real(a111*conjg(a111),dp)+real(a112*conjg(a112),dp)
     &      +real(a221*conjg(a221),dp)+real(a222*conjg(a222),dp)
     &      +real(a122*conjg(a122),dp)+real(a212*conjg(a212),dp)
     &      +real(a121*conjg(a121),dp)+real(a211*conjg(a211),dp))
            mqq(2,j,k)=half*faclo*(
     &       real(b111*conjg(b111),dp)+real(b112*conjg(b112),dp)
     &      +real(b221*conjg(b221),dp)+real(b222*conjg(b222),dp)
     &      +real(b122*conjg(b122),dp)+real(b212*conjg(b212),dp)
     &      +real(b121*conjg(b121),dp)+real(b211*conjg(b211),dp))
     
            endif
          elseif ((j < 0) .and. (k < 0)) then
c----QbQb case
            if (j .ne. k) then
            a111=(Q(-j))*qbRb_a(1,1,1)
     &          +(Q(-k))*qbRb_b(1,1,1)
            a121=(Q(-j))*qbRb_a(1,2,1)
     &          +(Q(-k))*qbRb_b(1,2,1)

            a112=(Q(-j))*qbRb_a(1,1,2)
     &          +(Q(-k))*qbRb_b(1,1,2)
            a122=(Q(-j))*qbRb_a(1,2,2)
     &          +(Q(-k))*qbRb_b(1,2,2)

            a211=(Q(-j))*qbRb_a(2,1,1)
     &          +(Q(-k))*qbRb_b(2,1,1)
            a221=(Q(-j))*qbRb_a(2,2,1)
     &          +(Q(-k))*qbRb_b(2,2,1)

            a212=(Q(-j))*qbRb_a(2,1,2)
     &          +(Q(-k))*qbRb_b(2,1,2)
            a222=(Q(-j))*qbRb_a(2,2,2)
     &          +(Q(-k))*qbRb_b(2,2,2)
            mqq(0,j,k)=zip
            mqq(1,j,k)=faclo*(
     &        real(a111*conjg(a111),dp)+real(a112*conjg(a112),dp)
     &       +real(a221*conjg(a221),dp)+real(a222*conjg(a222),dp)
     &       +real(a122*conjg(a122),dp)+real(a212*conjg(a212),dp)
     &       +real(a121*conjg(a121),dp)+real(a211*conjg(a211),dp))
            mqq(2,j,k)=zip

            elseif (j == k) then

            a111=(Q(-j))*(qbRb_a(1,1,1)+qbRb_b(1,1,1))
            b111=(Q(-j))*(qbqb_a(1,1,1)+qbqb_b(1,1,1))
            a112=(Q(-j))*(qbRb_a(1,1,2)+qbRb_b(1,1,2))
            b112=(Q(-j))*(qbqb_a(1,1,2)+qbqb_b(1,1,2))
            a221=(Q(-j))*(qbRb_a(2,2,1)+qbRb_b(2,2,1))
            b221=(Q(-j))*(qbqb_a(2,2,1)+qbqb_b(2,2,1))
            a222=(Q(-j))*(qbRb_a(2,2,2)+qbRb_b(2,2,2))
            b222=(Q(-j))*(qbqb_a(2,2,2)+qbqb_b(2,2,2))


            a121=(Q(-j))*qbRb_a(1,2,1)
     &          +(Q(-k))*qbRb_b(1,2,1)
            a122=(Q(-j))*qbRb_a(1,2,2)
     &          +(Q(-k))*qbRb_b(1,2,2)
            a211=(Q(-j))*qbRb_a(2,1,1)
     &          +(Q(-k))*qbRb_b(2,1,1)
            a212=(Q(-j))*qbRb_a(2,1,2)
     &          +(Q(-k))*qbRb_b(2,1,2)

            b121=(Q(-j))*qbqb_a(1,2,1)
     &          +(Q(-k))*qbqb_b(1,2,1)
            b122=(Q(-j))*qbqb_a(1,2,2)
     &          +(Q(-k))*qbqb_b(1,2,2)
            b211=(Q(-j))*qbqb_a(2,1,1)
     &          +(Q(-k))*qbqb_b(2,1,1)
            b212=(Q(-j))*qbqb_a(2,1,2)
     &          +(Q(-k))*qbqb_b(2,1,2)


            mqq(0,j,k)=half*faclo*(
     &      +real(a111*conjg(b111),dp)+real(a112*conjg(b112),dp)
     &      +real(a221*conjg(b221),dp)+real(a222*conjg(b222),dp))*two/xn
            mqq(1,j,k)=half*faclo*(
     &       real(a111*conjg(a111),dp)+real(a112*conjg(a112),dp)
     &      +real(a221*conjg(a221),dp)+real(a222*conjg(a222),dp)
     &      +real(a122*conjg(a122),dp)+real(a212*conjg(a212),dp)
     &      +real(a121*conjg(a121),dp)+real(a211*conjg(a211),dp))
            mqq(2,j,k)=half*faclo*(
     &       real(b111*conjg(b111),dp)+real(b112*conjg(b112),dp)
     &      +real(b221*conjg(b221),dp)+real(b222*conjg(b222),dp)
     &      +real(b122*conjg(b122),dp)+real(b212*conjg(b212),dp)
     &      +real(b121*conjg(b121),dp)+real(b211*conjg(b211),dp))
     
            endif
C---q-qb case
         elseif ((j > 0) .and. (k < 0)) then
             if (j .ne. -k) then 
            a111=(Q(+j))*qRb_a(1,1,1)
     &          +(Q(-k))*qRb_b(1,1,1)
            a112=(Q(+j))*qRb_a(1,1,2)
     &          +(Q(-k))*qRb_b(1,1,2)
            a221=(Q(+j))*qRb_a(2,2,1)
     &          +(Q(-k))*qRb_b(2,2,1)
            a222=(Q(+j))*qRb_a(2,2,2)
     &          +(Q(-k))*qRb_b(2,2,2)

            a121=(Q(+j))*qRb_a(1,2,1)
     &          +(Q(-k))*qRb_b(1,2,1)
            a122=(Q(+j))*qRb_a(1,2,2)
     &          +(Q(-k))*qRb_b(1,2,2)
            a211=(Q(+j))*qRb_a(2,1,1)
     &          +(Q(-k))*qRb_b(2,1,1)
            a212=(Q(+j))*qRb_a(2,1,2)
     &          +(Q(-k))*qRb_b(2,1,2)
            mqq(0,j,k)=zip
            mqq(1,j,k)=zip
            mqq(2,j,k)=faclo*(
     &        real(a111*conjg(a111),dp)+real(a112*conjg(a112),dp)
     &       +real(a221*conjg(a221),dp)+real(a222*conjg(a222),dp)
     &       +real(a122*conjg(a122),dp)+real(a212*conjg(a212),dp)
     &       +real(a121*conjg(a121),dp)+real(a211*conjg(a211),dp))

            elseif (j == -k) then
c--case where final state from annihilation diagrams is the same quark
            a111=(Q(j))*(qRb_a(1,1,1)+qRb_b(1,1,1))
            b111=(Q(j))*(qqb_a(1,1,1)+qqb_b(1,1,1))

            a112=(Q(j))*(qRb_a(1,1,2)+qRb_b(1,1,2))
            b112=(Q(j))*(qqb_a(1,1,2)+qqb_b(1,1,2))

            a221=(Q(j))*(qRb_a(2,2,1)+qRb_b(2,2,1))
            b221=(Q(j))*(qqb_a(2,2,1)+qqb_b(2,2,1))

            a222=(Q(j))*(qRb_a(2,2,2)+qRb_b(2,2,2))
            b222=(Q(j))*(qqb_a(2,2,2)+qqb_b(2,2,2))

            a121=(Q(+j))*qRb_a(1,2,1)
     &          +(Q(-k))*qRb_b(1,2,1)
            a122=(Q(+j))*qRb_a(1,2,2)
     &          +(Q(-k))*qRb_b(1,2,2)
            a211=(Q(+j))*qRb_a(2,1,1)
     &          +(Q(-k))*qRb_b(2,1,1)
            a212=(Q(+j))*qRb_a(2,1,2)
     &          +(Q(-k))*qRb_b(2,1,2)

            b121=(Q(+j))*qqb_a(1,2,1)
     &          +(Q(-k))*qqb_b(1,2,1)
            b122=(Q(+j))*qqb_a(1,2,2)
     &          +(Q(-k))*qqb_b(1,2,2)
            b211=(Q(+j))*qqb_a(2,1,1)
     &          +(Q(-k))*qqb_b(2,1,1)
            b212=(Q(+j))*qqb_a(2,1,2)
     &           +(Q(-k))*qqb_b(2,1,2)
!            goto 200
            mqq(0,j,k)=faclo*(
     &      +real(a111*conjg(b111),dp)+real(a112*conjg(b112),dp)
     &      +real(a221*conjg(b221),dp)+real(a222*conjg(b222),dp))*two/xn
            mqq(1,j,k)=faclo*(
     &       real(a111*conjg(a111),dp)+real(a112*conjg(a112),dp)
     &      +real(a221*conjg(a221),dp)+real(a222*conjg(a222),dp)
     &      +real(a122*conjg(a122),dp)+real(a212*conjg(a212),dp)
     &      +real(a121*conjg(a121),dp)+real(a211*conjg(a211),dp))
            mqq(2,j,k)=faclo*(
     &       real(b111*conjg(b111),dp)+real(b112*conjg(b112),dp)
     &      +real(b221*conjg(b221),dp)+real(b222*conjg(b222),dp)
     &      +real(b122*conjg(b122),dp)+real(b212*conjg(b212),dp)
     &      +real(b121*conjg(b121),dp)+real(b211*conjg(b211),dp))
 200        continue
       if ((j==1).or.(j==3).or.(j==5)) then
           nup=2
           ndo=nf-3
       else
           nup=1
           ndo=nf-2
       endif
       if (nflav <= 4) ndo=ndo-1
       if (nflav <= 3) nup=nup-1
            b111=(Q(+j))*qqb_a(1,1,1)
     &          +(Q(+1))*qqb_b(1,1,1)
            b112=(Q(+j))*qqb_a(1,1,2)
     &          +(Q(+1))*qqb_b(1,1,2)
            b221=(Q(+j))*qqb_a(2,2,1)
     &          +(Q(+1))*qqb_b(2,2,1)
            b222=(Q(+j))*qqb_a(2,2,2)
     &          +(Q(+1))*qqb_b(2,2,2)
            b121=(Q(+j))*qqb_a(1,2,1)
     &          +(Q(+1))*qqb_b(1,2,1)
            b122=(Q(+j))*qqb_a(1,2,2)
     &          +(Q(+1))*qqb_b(1,2,2)
            b211=(Q(+j))*qqb_a(2,1,1)
     &          +(Q(+1))*qqb_b(2,1,1)
            b212=(Q(+j))*qqb_a(2,1,2)
     &          +(Q(+1))*qqb_b(2,1,2)
            
      tdo=faclo*(
     &     real(b111*conjg(b111),dp)+real(b112*conjg(b112),dp)
     &    +real(b221*conjg(b221),dp)+real(b222*conjg(b222),dp)
     &    +real(b122*conjg(b122),dp)+real(b212*conjg(b212),dp)
     &    +real(b121*conjg(b121),dp)+real(b211*conjg(b211),dp))

            b111=(Q(+j))*qqb_a(1,1,1)
     &          +(Q(+2))*qqb_b(1,1,1)
            b112=(Q(+j))*qqb_a(1,1,2)
     &          +(Q(+2))*qqb_b(1,1,2)
            b221=(Q(+j))*qqb_a(2,2,1)
     &          +(Q(+2))*qqb_b(2,2,1)
            b222=(Q(+j))*qqb_a(2,2,2)
     &          +(Q(+2))*qqb_b(2,2,2)
            b121=(Q(+j))*qqb_a(1,2,1)
     &          +(Q(+2))*qqb_b(1,2,1)
            b122=(Q(+j))*qqb_a(1,2,2)
     &          +(Q(+2))*qqb_b(1,2,2)
            b211=(Q(+j))*qqb_a(2,1,1)
     &          +(Q(+2))*qqb_b(2,1,1)
            b212=(Q(+j))*qqb_a(2,1,2)
     &          +(Q(+2))*qqb_b(2,1,2)
            
      tup=faclo*(
     &     real(b111*conjg(b111),dp)+real(b112*conjg(b112),dp)
     &    +real(b221*conjg(b221),dp)+real(b222*conjg(b222),dp)
     &    +real(b122*conjg(b122),dp)+real(b212*conjg(b212),dp)
     &    +real(b121*conjg(b121),dp)+real(b211*conjg(b211),dp))

      mqq(1,j,k)=mqq(1,j,k)+real(nup,dp)*tup+real(ndo,dp)*tdo
      endif
      elseif ((j < 0) .and. (k > 0)) then
C---Qb-q case
            if (j .ne. -k) then
            a111=(Q(-j))*qbR_a(1,1,1)
     &          +(Q(+k))*qbR_b(1,1,1)
            a121=(Q(-j))*qbR_a(1,2,1)
     &          +(Q(+k))*qbR_b(1,2,1)
            a112=(Q(-j))*qbR_a(1,1,2)
     &          +(Q(+k))*qbR_b(1,1,2)
            a122=(Q(-j))*qbR_a(1,2,2)
     &          +(Q(+k))*qbR_b(1,2,2)
            a211=(Q(-j))*qbR_a(2,1,1)
     &          +(Q(+k))*qbR_b(2,1,1)
            a221=(Q(-j))*qbR_a(2,2,1)
     &          +(Q(+k))*qbR_b(2,2,1)
            a212=(Q(-j))*qbR_a(2,1,2)
     &          +(Q(+k))*qbR_b(2,1,2)
            a222=(Q(-j))*qbR_a(2,2,2)
     &          +(Q(+k))*qbR_b(2,2,2)

            mqq(0,j,k)=zip
            mqq(1,j,k)=zip
            mqq(2,j,k)=faclo*(
     &        real(a111*conjg(a111),dp)+real(a112*conjg(a112),dp)
     &       +real(a221*conjg(a221),dp)+real(a222*conjg(a222),dp)
     &       +real(a122*conjg(a122),dp)+real(a212*conjg(a212),dp)
     &       +real(a121*conjg(a121),dp)+real(a211*conjg(a211),dp))
            elseif (j == -k) then

            a111=(Q(-j))*(qbR_a(1,1,1)+qbR_b(1,1,1))
            b111=(Q(-j))*(qbq_a(1,1,1)+qbq_b(1,1,1))
            a112=(Q(-j))*(qbR_a(1,1,2)+qbR_b(1,1,2))
            b112=(Q(-j))*(qbq_a(1,1,2)+qbq_b(1,1,2))
            a221=(Q(-j))*(qbR_a(2,2,1)+qbR_b(2,2,1))
            b221=(Q(-j))*(qbq_a(2,2,1)+qbq_b(2,2,1))
            a222=(Q(-j))*(qbR_a(2,2,2)+qbR_b(2,2,2))
            b222=(Q(-j))*(qbq_a(2,2,2)+qbq_b(2,2,2))

            a121=(Q(-j))*qbR_a(1,2,1)
     &          +(Q(+k))*qbR_b(1,2,1)
            a122=(Q(-j))*qbR_a(1,2,2)
     &          +(Q(+k))*qbR_b(1,2,2)
            a211=(Q(-j))*qbR_a(2,1,1)
     &          +(Q(+k))*qbR_b(2,1,1)
            a212=(Q(-j))*qbR_a(2,1,2)
     &          +(Q(+k))*qbR_b(2,1,2)

            b121=(Q(-j))*qbq_a(1,2,1)
     &          +(Q(+k))*qbq_b(1,2,1)
            b122=(Q(-j))*qbq_a(1,2,2)
     &          +(Q(+k))*qbq_b(1,2,2)
            b211=(Q(-j))*qbq_a(2,1,1)
     &          +(Q(+k))*qbq_b(2,1,1)
            b212=(Q(-j))*qbq_a(2,1,2)
     &          +(Q(+k))*qbq_b(2,1,2)
            mqq(0,j,k)=faclo*(
     &      +real(a111*conjg(b111),dp)+real(a112*conjg(b112),dp)
     &      +real(a221*conjg(b221),dp)+real(a222*conjg(b222),dp))*two/xn
            mqq(2,j,k)=faclo*(
     &       real(a111*conjg(a111),dp)+real(a112*conjg(a112),dp)
     &      +real(a221*conjg(a221),dp)+real(a222*conjg(a222),dp)
     &      +real(a122*conjg(a122),dp)+real(a212*conjg(a212),dp)
     &      +real(a121*conjg(a121),dp)+real(a211*conjg(a211),dp))
            mqq(1,j,k)=faclo*(
     &       real(b111*conjg(b111),dp)+real(b112*conjg(b112),dp)
     &      +real(b221*conjg(b221),dp)+real(b222*conjg(b222),dp)
     &      +real(b122*conjg(b122),dp)+real(b212*conjg(b212),dp)
     &      +real(b121*conjg(b121),dp)+real(b211*conjg(b211),dp))
c--Here we must also add the contribution of other final state quarks
c  unequal to initial annihilating quarks
       if ((k==1).or.(k==3).or.(k==5)) then
           nup=2
           ndo=nf-3
       else
           nup=1
           ndo=nf-2
       endif
       if (nflav <= 4) ndo=ndo-1
       if (nflav <= 3) nup=nup-1
            b111=(Q(-j))*qbq_a(1,1,1)
     &          +(Q(+3))*qbq_b(1,1,1)
            b112=(Q(-j))*qbq_a(1,1,2)
     &          +(Q(+3))*qbq_b(1,1,2)
            b221=(Q(-j))*qbq_a(2,2,1)
     &          +(Q(+3))*qbq_b(2,2,1)
            b222=(Q(-j))*qbq_a(2,2,2)
     &          +(Q(+3))*qbq_b(2,2,2)
            b121=(Q(-j))*qbq_a(1,2,1)
     &          +(Q(+3))*qbq_b(1,2,1)
            b122=(Q(-j))*qbq_a(1,2,2)
     &          +(Q(+3))*qbq_b(1,2,2)
            b211=(Q(-j))*qbq_a(2,1,1)
     &          +(Q(+3))*qbq_b(2,1,1)
            b212=(Q(-j))*qbq_a(2,1,2)
     &          +(Q(+3))*qbq_b(2,1,2)
      tdo=faclo*(
     &     real(b111*conjg(b111),dp)+real(b112*conjg(b112),dp)
     &    +real(b221*conjg(b221),dp)+real(b222*conjg(b222),dp)
     &    +real(b122*conjg(b122),dp)+real(b212*conjg(b212),dp)
     &    +real(b121*conjg(b121),dp)+real(b211*conjg(b211),dp))

            b111=(Q(-j))*qbq_a(1,1,1)
     &          +(Q(+2))*qbq_b(1,1,1)
            b112=(Q(-j))*qbq_a(1,1,2)
     &          +(Q(+2))*qbq_b(1,1,2)
            b221=(Q(-j))*qbq_a(2,2,1)
     &          +(Q(+2))*qbq_b(2,2,1)
            b222=(Q(-j))*qbq_a(2,2,2)
     &          +(Q(+2))*qbq_b(2,2,2)
            b121=(Q(-j))*qbq_a(1,2,1)
     &          +(Q(+2))*qbq_b(1,2,1)
            b122=(Q(-j))*qbq_a(1,2,2)
     &          +(Q(+2))*qbq_b(1,2,2)
            b211=(Q(-j))*qbq_a(2,1,1)
     &          +(Q(+2))*qbq_b(2,1,1)
            b212=(Q(-j))*qbq_a(2,1,2)
     &          +(Q(+2))*qbq_b(2,1,2)
      tup=faclo*(
     &     real(b111*conjg(b111),dp)+real(b112*conjg(b112),dp)
     &    +real(b221*conjg(b221),dp)+real(b222*conjg(b222),dp)
     &    +real(b122*conjg(b122),dp)+real(b212*conjg(b212),dp)
     &    +real(b121*conjg(b121),dp)+real(b211*conjg(b211),dp))

      mqq(1,j,k)=mqq(1,j,k)+real(nup,dp)*tup+real(ndo,dp)*tdo

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
