!---  CW June 16
!-----New version of tree-level process for gamma + 2j, written in the style of
!---- Z2jet
      subroutine qqb_gam2j(p,msq)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'constants.f'
      include 'nf.f'
      include 'qcdcouple.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      include 'nflav.f'
      include 'ewcharge.f'
      include 'ewcouple.f'
      include 'msq_cs.f'
      include 'flags.f'
      
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
      complex(dp):: qqb_a(2,2,2),qqb_b(2,2,2)

      complex(dp):: qbq_a(2,2,2),qbq_b(2,2,2)
      complex(dp):: qbR_a(2,2,2),qbR_b(2,2,2)

      complex(dp):: qq_a(2,2,2),qq_b(2,2,2)
      complex(dp):: qR_a(2,2,2),qR_b(2,2,2)

      complex(dp):: qbRb_a(2,2,2),qbRb_b(2,2,2)
      complex(dp):: qbqb_a(2,2,2),qbqb_b(2,2,2)
      
      real(dp):: mqq(0:2,fn:nf,fn:nf)
      common/mqq/mqq
!$omp threadprivate(/mqq/)
      integer,parameter::swap(2)=(/2,1/),swap1(0:2)=(/0,2,1/)
      integer:: i,j,k,nup,ndo,icol,j1,j2,j3
      real(dp):: Qsum_all
      real(dp):: test,Bigagam
      
      msq(:,:)=zip
      msq_cs(:,:,:)=zip

!---  debug
!      Q(-4)=zip
!      Q(-5)=zip
!      Q(4)=zip
!      Q(5)=zip
!      write(6,*) Q
!      Test=Bigagam(4,5,1,2,3,-5,-4)
 
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
      
!-------calculate 2-quark, 2-gluon amplitudes
      if (Gflag) then

        call qqbgg_ga(2,1,4,5,3,za,zb,qqb_gagg)
        call storecs_ga(qqb_gagg_cs)
        call qqbgg_ga(4,1,2,5,3,za,zb,qg_gaqg)
        call storecs_ga(qg_gaqg_cs)
        call qqbgg_ga(4,2,1,5,3,za,zb,gq_gaqg)
        call storecs_ga(gq_gaqg_cs)

        do icol=0,2
        qbq_gagg_cs(icol)=qqb_gagg_cs(swap1(icol))
        qbg_gaqbg_cs(icol)=qg_gaqg_cs(swap1(icol))
        gqb_gaqbg_cs(icol)=gq_gaqg_cs(swap1(icol))
        enddo
       
        call qqbgg_ga(5,4,1,2,3,za,zb,gg_gaqqb)
        call storecs_ga(gg_gaqqb_cs)

        fac=2._dp*esq*gsq**2*V*xn

        qqb_gagg_cs(:)=half*aveqq*fac*qqb_gagg_cs(:)
        qbq_gagg_cs(:)=half*aveqq*fac*qbq_gagg_cs(:)
        gq_gaqg_cs(:)=aveqg*fac*gq_gaqg_cs(:)
        qg_gaqg_cs(:)=aveqg*fac*qg_gaqg_cs(:)
        gqb_gaqbg_cs(:)=aveqg*fac*gqb_gaqbg_cs(:)
        qbg_gaqbg_cs(:)=aveqg*fac*qbg_gaqbg_cs(:)
        gg_gaqqb_cs(:)=avegg*fac*gg_gaqqb_cs(:)
        
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

        faclo=2._dp*V*gsq**2*esq*aveqq
        
      endif

!---- assemble gg pieces
      if (Gflag) then
      
        do j=-nf,nf
        do k=-nf,nf

          if( j .ne. 0 .and. k .ne. 0 .and. j .ne. -k) cycle
          if     ((j == 0) .and. (k == 0)) then
             msq(j,k)=Qsum_all*gg_gaqqb
             msq_cs(:,j,k)=Qsum_all*gg_gaqqb_cs(:)
          elseif ((j > 0) .and. (k < 0)) then
             msq(j,k)=qqb_gagg*Q(j)**2
             msq_cs(:,j,k)=qqb_gagg_cs(:)*Q(j)**2
          elseif ((j < 0) .and. (k > 0)) then
             msq(j,k)=qbq_gagg*Q(-j)**2
             msq_cs(:,j,k)=qbq_gagg_cs(:)*Q(-j)**2
          elseif ((j > 0) .and. (k == 0)) then
             msq(j,k)=qg_gaqg*Q(j)**2
             msq_cs(:,j,k)=qg_gaqg_cs(:)*Q(j)**2
          elseif ((j < 0) .and. (k == 0)) then
             msq(j,k)=qbg_gaqbg*Q(-j)**2
              msq_cs(:,j,k)=qbg_gaqbg_cs(:)*Q(-j)**2
          elseif ((j == 0) .and. (k > 0)) then
             msq(j,k)=gq_gaqg*Q(k)**2
             msq_cs(:,j,k)=gq_gaqg_cs(:)*Q(k)**2
          elseif ((j == 0) .and. (k < 0)) then
             msq(j,k)=gqb_gaqbg*Q(-k)**2
             msq_cs(:,j,k)=gqb_gaqbg_cs(:)*Q(-k)**2
          endif
            
        enddo
        enddo
      endif

      if (Qflag) then

        mqq(:,:,:)=zip

        do j=-nf,nf
        do k=-nf,nf
            
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
     &          +(Q(-k))*qqb_b(2,1,2)

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

!      write(6,*) test,msq(-5,-4)/faclo,msq(-5,-4)/test

      return
      end
