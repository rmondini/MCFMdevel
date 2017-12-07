!---  CW Virtual routine for Gamma + 2 j
!---- June 16,

      subroutine qqb_gam2j_v(p,msqv)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'ewcharge.f'
      include 'zcouple.f'
      include 'masses.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'scheme.f'
      include 'scale.f'
      include 'flags.f'
      include 'msq_cs.f'
      include 'lc.f'
      include 'first.f'
      include 'ppmax.f'
      include 'mpicommon.f'
      real(dp):: p(mxpart,4)
      real(dp):: msq0(-nf:nf,-nf:nf),msqv(-nf:nf,-nf:nf)
      real(dp):: mqq(0:2,-nf:nf,-nf:nf),msqx_cs(0:2,-nf:nf,-nf:nf)
      real(dp):: subuv(0:2)
      real(dp):: ppmsqx(0:2,ppmax)
      integer i
      common/mqq/mqq
      integer, parameter :: jj(-nf:nf)=(/-1,-2,-1,-2,-1,0,1,2,1,2,1/)
      integer:: nup,ndo,rvcolourchoice
      common/rvcolourchoice/rvcolourchoice
!$omp threadprivate(/mqq/,/rvcolourchoice/)
      parameter (nup=2,ndo=nf-nup)
      real(dp):: Qsum_all
      real(dp):: fac,facnlo
      real(dp):: ga4q_noid,ga4q_id
      real(dp) :: qqb_gagg,qqb_gagg_cs(0:2)
      real(dp) :: qbq_gagg,qbq_gagg_cs(0:2)
      real(dp) :: qg_gaqg,qg_gaqg_cs(0:2)
      real(dp) :: qbg_gaqbg,qbg_gaqbg_cs(0:2)
      real(dp) :: gq_gaqg,gq_gaqg_cs(0:2)
      real(dp) :: gqb_gaqbg,gqb_gaqbg_cs(0:2)
      real(dp) :: gg_gaqqb, gg_gaqqb_cs(0:2)
      integer j,k,icol,cs

      
      complex(dp):: qR_v1_a(2,2,2),qR_v2_a(2,2,2)
     &     ,qR_v1_b(2,2,2),qR_v2_b(2,2,2)
      complex(dp):: qR_a(2,2,2),qR_b(2,2,2)
      complex(dp):: qq_v1_a(2,2,2),qq_v2_a(2,2,2)
     &     ,qq_v1_b(2,2,2),qq_v2_b(2,2,2)
      complex(dp):: qq_a(2,2,2),qq_b(2,2,2)
      
      complex(dp):: qbR_v1_a(2,2,2),qbR_v2_a(2,2,2)
     &     ,qbR_v1_b(2,2,2),qbR_v2_b(2,2,2)
      complex(dp):: qbR_a(2,2,2),qbR_b(2,2,2)
      complex(dp):: qbq_v1_a(2,2,2),qbq_v2_a(2,2,2)
     &     ,qbq_v1_b(2,2,2),qbq_v2_b(2,2,2)
      complex(dp):: qbq_a(2,2,2),qbq_b(2,2,2)

      complex(dp):: qRb_v1_a(2,2,2),qRb_v2_a(2,2,2)
     &     ,qRb_v1_b(2,2,2),qRb_v2_b(2,2,2)
      complex(dp):: qRb_a(2,2,2),qRb_b(2,2,2)
      complex(dp):: qqb_v1_a(2,2,2),qqb_v2_a(2,2,2)
     &     ,qqb_v1_b(2,2,2),qqb_v2_b(2,2,2)
      complex(dp):: qqb_a(2,2,2),qqb_b(2,2,2)

      complex(dp):: qbRb_v1_a(2,2,2),qbRb_v2_a(2,2,2)
     &     ,qbRb_v1_b(2,2,2),qbRb_v2_b(2,2,2)
      complex(dp):: qbRb_a(2,2,2),qbRb_b(2,2,2)
      complex(dp):: qbqb_v1_a(2,2,2),qbqb_v2_a(2,2,2)
     &     ,qbqb_v1_b(2,2,2),qbqb_v2_b(2,2,2)
      complex(dp):: qbqb_a(2,2,2),qbqb_b(2,2,2)

      scheme = 'dred'

      if (Qflag .and. Gflag) then
        write(6,*) 'Both Qflag and Gflag cannot be true'
        write(6,*) 'They are set in file options.DAT'
        write(6,*) 'Failed in qqb_z2jet_v.f'
        stop
      endif
      if (first) then
        first=.false.
        if (rank == 0) then
        if ((Gflag) .or. (QandGflag)) then
          write(*,*) 'Using QQGG (VIRTUAL) matrix elements'
!          write(*,*) '[LC is     N   ]'
!          write(*,*) '[SLC is   1/N  ]'
!          write(*,*) '[SSLC is 1/N**3]'
        endif
        if ((Qflag) .or. (QandGflag)) then
          write(*,*) 'Using QQBQQB (VIRTUAL) matrix elements'
!          write(*,*) '[LC is   1 ]'
!          write(*,*) '[SLC is 1/N]'
        endif
        if     (rvcolourchoice == 1) then
          write(*,*) 'Leading colour only in VIRTUAL'
        elseif (rvcolourchoice == 2) then
          write(*,*) 'Sub-leading colour only in VIRTUAL'
        elseif (rvcolourchoice == 3) then
          write(*,*) 'Sub-sub-leading colour only in VIRTUAL'
        elseif (rvcolourchoice == 0) then
          write(*,*) 'Total of all colour structures in VIRTUAL'
        else
          write(*,*) 'Bad colourchoice'
          stop
        endif
        endif
      endif
      
      Qsum_all=zip
      do j=1,nf
         Qsum_all=Qsum_all+Q(j)**2
      enddo


 

      do j=-nf,nf
         do k=-nf,nf
            msqv(j,k)=0._dp
         enddo
      enddo
!      musq=1._dp
!      epinv=1._dp
!      include 'kinpoint_qqbggga_kc.f'
!      call spinorz(5,p,za,zb)

      call spinoru(5,p,za,zb)      
      call qqb_gam2jx_new(p,msq0,mqq,ppmsqx,msqx_cs)
 

      if (Gflag) then
c----UV counterterm contains the finite renormalization to arrive
c----at MS bar scheme. 
      if     (colourchoice == 1) then
        subuv(1)=2._dp*xn*(epinv*(11._dp-2._dp*real(nf,dp)/xn)-1._dp)/6._dp
        subuv(2)=subuv(1)
      elseif (colourchoice == 2) then
        subuv(0)=2._dp*xn*(epinv*(11._dp-2._dp*real(nf,dp)/xn)-1._dp)/6._dp
      elseif (colourchoice == 3) then
c--- all zero already
      elseif (colourchoice == 0) then
        subuv(1)=2._dp*xn*(epinv*(11._dp-2._dp*real(nf,dp)/xn)-1._dp)/6._dp
        subuv(2)=subuv(1)
        subuv(0)=subuv(1)
      endif


c--- transfer lowest order matrix elements
c--- NB: this breaks the routine if Qflag = Gflag = .true.

      do cs=0,2
        do j=-nf,nf
        do k=-nf,nf
        msq_cs(cs,j,k)=msqx_cs(cs,j,k)
        enddo
        enddo
      enddo

!----- CW like LO, should adjust these names
      call qqbgg_ga_v(2,1,4,5,3,za,zb,qqb_gagg)
      call qqbgg_ga_v(4,1,2,5,3,za,zb,qg_gaqg)
      call qqbgg_ga_v(4,2,1,5,3,za,zb,gq_gaqg)
     
      call qqbgg_ga_v(5,4,1,2,3,za,zb,gg_gaqqb)

      
       qbq_gagg=qqb_gagg
       qbg_gaqbg=qg_gaqg
       gqb_gaqbg=gq_gaqg
       
      fac=2._dp*esq*gsq**2*V*ason2pi


      qqb_gagg=fac*half*aveqq*qqb_gagg
      qbq_gagg=fac*half*aveqq*qbq_gagg
      gq_gaqg=aveqg*fac*gq_gaqg
      qg_gaqg=aveqg*fac*qg_gaqg
      
      gqb_gaqbg=aveqg*fac*gqb_gaqbg
      qbg_gaqbg=aveqg*fac*qbg_gaqbg
      gg_gaqqb=avegg*fac*gg_gaqqb

      
      endif
      
************************************************************************
*     Endpoint contributions from QQQQ matrix elements                 *
************************************************************************            
      if (Qflag) then
        subuv(1)=2._dp*xn*(epinv*(11._dp-2._dp*real(nf,dp)/xn)-1._dp)/6._dp
        subuv(2)=subuv(1)
        subuv(0)=subuv(1)
      
c--- transfer lowest order matrix elements
c--- NB: this breaks the routine if Qflag = Gflag = .true.

      do cs=0,2
        do j=-nf,nf
        do k=-nf,nf
        msq_cs(cs,j,k)=mqq(cs,j,k)
        enddo
        enddo
      enddo
   
      endif

      if(Qflag)then
!      q Qb -> q Qb
       call ga_qqbQQb(1,5,4,2,za,zb,qRb_a,qRb_b)
       call ga_qqbQQb_v(1,5,4,2,3,za,zb,qRb_v1_a,qRb_v2_a
     &      ,qRb_v1_b,qRb_v2_b)

!-----q qb -> QQb
      call ga_qqbQQb(1,2,4,5,za,zb,qqb_a,qqb_b)
      call ga_qqbQQb_v(1,2,4,5,3,za,zb,qqb_v1_a,qqb_v2_a
     &     ,qqb_v1_b,qqb_v2_b)


!------ qb Q -> qb Q
       call ga_qqbQQb(5,1,2,4,za,zb,qbR_a,qbR_b)
       call ga_qqbQQb_v(5,1,2,4,3,za,zb,qbR_v1_a,qbR_v2_a
     &      ,qbR_v1_b,qbR_v2_b)
  
!-----qb q -> QQb
      call ga_qqbQQb(2,1,5,4,za,zb,qbq_a,qbq_b)
      call ga_qqbQQb_v(2,1,5,4,3,za,zb,qbq_v1_a,qbq_v2_a
     &     ,qbq_v1_b,qbq_v2_b)

            
!      q Q -> q Q 
      call ga_qqbQQb(1,4,2,5,za,zb,qR_a,qR_b)
      call ga_qqbQQb_v(1,4,2,5,3,za,zb,qR_v1_a,qR_v2_a,qR_v1_b,qR_v2_b)
      
!------ qb Qb -> qb Qb
       call ga_qqbQQb(4,1,5,2,za,zb,qbRb_a,qbRb_b)      
       call ga_qqbQQb_v(4,1,5,2,3,za,zb,qbRb_v1_a,qbRb_v2_a
     &      ,qbRb_v1_b,qbRb_v2_b)



!-----qq->qq
      call ga_qqbQQb(1,5,2,4,za,zb,qq_a,qq_b)      
      call ga_qqbQQb_v(1,5,2,4,3,za,zb,qq_v1_a,qq_v2_a,qq_v1_b,qq_v2_b)

!-----qbqb->qbqb
      call ga_qqbQQb(5,1,4,2,za,zb,qbqb_a,qbqb_b)
      call ga_qqbQQb_v(5,1,4,2,3,za,zb,qbqb_v1_a,qbqb_v2_a,qbqb_v1_b
     &     ,qbqb_v2_b)



      facnlo=2._dp*gsq**2*esq*aveqq*ason4pi
      endif
      
      if(Gflag) then 
       
!---- assemble gg pieces
      do j=-nf,nf
      do k=-nf,nf

         if( j .ne. 0 .and. k .ne. 0 .and. j .ne. -k) goto 19
         if     ((j == 0) .and. (k == 0)) then
            msqv(j,k)=Qsum_all*gg_gaqqb
         elseif ((j > 0) .and. (k < 0)) then
            msqv(j,k)=qqb_gagg*Q(j)**2
         elseif ((j < 0) .and. (k > 0)) then
            msqv(j,k)=qbq_gagg*Q(-j)**2
         elseif ((j > 0) .and. (k == 0)) then
            msqv(j,k)=qg_gaqg*Q(j)**2
         elseif ((j < 0) .and. (k == 0)) then
            msqv(j,k)=qbg_gaqbg*Q(-j)**2
         elseif ((j == 0) .and. (k > 0)) then
            msqv(j,k)=gq_gaqg*Q(k)**2
         elseif ((j == 0) .and. (k < 0)) then
            msqv(j,k)=gqb_gaqbg*Q(k)**2
         endif
        
         
 19      continue
      enddo
      enddo

      endif


      if(Qflag) then 
!----
     
      do j=-nf,nf
      do k=-nf,nf
      
      
          if ((j > 0) .and. (k > 0)) then
c---- QQ case
            if (j .ne. k) then
               msqv(j,k) =msqv(j,k) + ga4q_noid(j,k,qR_a,qR_b,qR_v1_a,qR_v1_b
     &              ,qR_v2_a,qR_v2_b)
            elseif(j==k) then 
!-------identical quarks
               msqv(j,k) =msqv(j,k) + ga4q_id(j,k,qR_a,qR_b,qq_a,qq_b,
     &  qR_v1_a,qR_v1_b,qR_v2_a,qR_v2_b,qq_v1_a,qq_v1_b,qq_v2_a,qq_v2_b)         
            endif

         elseif((j < 0) .and. (k<0)) then
            if(j.ne.k) then 
            msqv(j,k) =msqv(j,k) + ga4q_noid(j,k,qbRb_a,qbRb_b,qbRb_v1_a
     &             ,qbRb_v1_b ,qbRb_v2_a,qbRb_v2_b)
         elseif(j==k) then 
!-------identical quarks
               msqv(j,k) =msqv(j,k) + ga4q_id(j,k,qbRb_a,qbRb_b,qbqb_a,qbqb_b,
     &  qbRb_v1_a,qbRb_v1_b,qbRb_v2_a,qbRb_v2_b,qbqb_v1_a,qbqb_v1_b,qbqb_v2_a,qbqb_v2_b)            
            endif
         elseif((j>0).and.(k<0)) then
            if(j.ne.-k) then 
             msqv(j,k) =msqv(j,k) + ga4q_noid(j,k,qRb_a,qRb_b,qRb_v1_a
     &             ,qRb_v1_b ,qRb_v2_a,qRb_v2_b)
          elseif(j==-k) then
             do i = 1,nf
                if(i.ne.j) then 
             msqv(j,k) =msqv(j,k) + ga4q_noid(j,i,qqb_a,qqb_b,qqb_v1_a
     &                  ,qqb_v1_b ,qqb_v2_a,qqb_v2_b)
                else
              msqv(j,k)= msqv(j,k) +2._dp*ga4q_id(i,i,qRb_a,qRb_b,qqb_a,qqb_b,
     &  qRb_v1_a,qRb_v1_b,qRb_v2_a,qRb_v2_b,qqb_v1_a,qqb_v1_b,qqb_v2_a,qqb_v2_b) 
           endif
           enddo
             
          endif
       elseif((j<0).and.(k>0)) then
          if(j.ne.-k) then 
             msqv(j,k) =msqv(j,k) + ga4q_noid(j,k,qbR_a,qbR_b,qbR_v1_a
     &            ,qbR_v1_b ,qbR_v2_a,qbR_v2_b)
          elseif(j==-k) then
             
          do i=1,nf
          if(i.ne.k) then 
             msqv(j,k) =msqv(j,k) + ga4q_noid(k,i,qbq_a,qbq_b,qbq_v1_a
     &                  ,qbq_v1_b ,qbq_v2_a,qbq_v2_b)
                else
                  
              msqv(j,k)= msqv(j,k) + 2._dp*ga4q_id(j,j,qbR_a,qbR_b,qbq_a,qbq_b,
     &  qbR_v1_a,qbR_v1_b,qbR_v2_a,qbR_v2_b,qbq_v1_a,qbq_v1_b,qbq_v2_a,qbq_v2_b)  
           endif
           
           enddo
        endif
          
         endif
         
      enddo
      enddo

      endif

************************************************************************
*     UV contributions are included here                               *
C     This is the correction to put the answer in UV renormalized      *
C     dred scheme with msbar coupling                                  *
************************************************************************
      do j=-nf,nf
      do k=-nf,nf

         
      do cs=0,2
         msqv(j,k)=msqv(j,k)-ason2pi*subuv(cs)*msq_cs(cs,j,k)
      enddo

      
      enddo
      enddo
      return
      end

      function ga4q_noid(j,k,lo_a,lo_b,v_del1_a,v_del1_b
     &     ,v_del2_a,v_del2_b)
      implicit none
      include 'types.f'
      real(dp) ::  ga4q_noid
      include 'constants.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'nf.f'
      include 'ewcharge.f'
      include 'epinv2.f'
      complex(dp) :: lo_a(2,2,2),lo_b(2,2,2)
      complex(dp) :: v_del1_a(2,2,2),v_del2_a(2,2,2)
      complex(dp) :: v_del1_b(2,2,2),v_del2_b(2,2,2)
      complex(dp) :: lo_tot(2,2,2),v_del1_tot(2,2,2),v_del2_tot(2,2,2)
      integer :: j,k,a,b
      integer :: h1,h2,h3
      real(dp) :: fac
      a=abs(j)
      b=abs(k)

      fac=V*2._dp*gsq**2*esq*aveqq*ason4pi
      
!-----LO amplitude
      lo_tot(:,:,:) = lo_a(:,:,:)*Q(a)+lo_b(:,:,:)*Q(b)

!-----total del_1 piece
      v_del1_tot(:,:,:) = v_del1_a(:,:,:)*Q(a)+v_del1_b(:,:,:)*Q(b)
!----- and del_2 piece 
      v_del2_tot(:,:,:) = v_del2_a(:,:,:)*Q(a)+v_del2_b(:,:,:)*Q(b)

!---- note that v_del2 doesnt actually interfer with the LO stucture so we just get the following      
      ga4q_noid=zip

      do h1=1,2
         do h2=1,2
            do h3=1,2
               ga4q_noid = ga4q_noid 
     &            +real(conjg(lo_tot(h1,h2,h3))*v_del1_tot(h1,h2,h3),dp)
     &            +real(conjg(v_del1_tot(h1,h2,h3))*lo_tot(h1,h2,h3),dp)
               
            enddo
         enddo
      enddo
      ga4q_noid=ga4q_noid*fac
      return
      end


      
      function ga4q_id(j,k,lo_a,lo_b,lo2_a,lo2_b,v_del1_a,v_del1_b
     &     ,v_del2_a,v_del2_b,v2_del1_a,v2_del1_b
     &     ,v2_del2_a,v2_del2_b)
      implicit none
      include 'types.f'
      real(dp) ::  ga4q_id
      include 'constants.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'nf.f'
      include 'ewcharge.f'
      include 'epinv2.f'
      complex(dp) :: lo_a(2,2,2),lo_b(2,2,2)
      complex(dp) :: lo2_a(2,2,2),lo2_b(2,2,2)
      complex(dp) :: v_del1_a(2,2,2),v_del2_a(2,2,2)
      complex(dp) :: v_del1_b(2,2,2),v_del2_b(2,2,2)
      complex(dp) :: v2_del1_a(2,2,2),v2_del2_a(2,2,2)
      complex(dp) :: v2_del1_b(2,2,2),v2_del2_b(2,2,2)
      complex(dp) :: lo_tot(2,2,2),v1_del1_tot(2,2,2),v1_del2_tot(2,2,2)
      complex(dp) :: lo2_tot(2,2,2),v2_del1_tot(2,2,2),v2_del2_tot(2,2,2)
      complex(dp) :: vp1(2,2,2),vp2(2,2,2)
      integer :: j,k,a,b
      integer :: h1,h2,h3
      real(dp) :: fac
      a=abs(j)
      b=abs(k)

!----- overall color factor is still V, but with half for identical quarks
      fac=V*2._dp*gsq**2*esq*aveqq*ason4pi*half 
      
!-----LO amplitude
      lo_tot(:,:,:) = lo_a(:,:,:)*Q(a)+lo_b(:,:,:)*Q(b)
      lo2_tot(:,:,:) = (lo2_a(:,:,:)*Q(a)+lo2_b(:,:,:)*Q(b))

!-----v1, total del_1 piece
      v1_del1_tot(:,:,:) = v_del1_a(:,:,:)*Q(a)+v_del1_b(:,:,:)*Q(b)
!----- and del_2 piece 
      v1_del2_tot(:,:,:) = v_del2_a(:,:,:)*Q(a)+v_del2_b(:,:,:)*Q(b)


!-----v2, total del_1 piece
      v2_del1_tot(:,:,:) = v2_del1_a(:,:,:)*Q(a)+v2_del1_b(:,:,:)*Q(b)
!----- v2, and del_2 piece 
      v2_del2_tot(:,:,:) = v2_del2_a(:,:,:)*Q(a)+v2_del2_b(:,:,:)*Q(b)

      vp1(:,:,:)=v1_del1_tot(:,:,:)+v2_del1_tot(:,:,:)
      vp2(:,:,:)=+v1_del2_tot(:,:,:)+v2_del2_tot(:,:,:)
      
!---  answer is
!(-1 + Nc^2)  (ALOC2 (v1del2 + v2de1l) + ALOC (v1de1l + v2del2))

      ga4q_id=zip

!==== two copies of non-identical bit
      do h1=1,2
         do h2=1,2
            do h3=1,2
               ga4q_id = ga4q_id 
     &           +real(conjg(lo_tot(h1,h2,h3))*v1_del1_tot(h1,h2,h3),dp)
     &           +real(conjg(v1_del1_tot(h1,h2,h3))*lo_tot(h1,h2,h3),dp)
     &          +real(conjg(lo2_tot(h1,h2,h3))*v2_del1_tot(h1,h2,h3),dp)
     &          +real(conjg(v2_del1_tot(h1,h2,h3))*lo2_tot(h1,h2,h3),dp)
               
            enddo
         enddo
      enddo
!===== 1/xn intf
      do h1=1,2
         h2=h1
!     do h2=1,2
            do h3=1,2
               ga4q_id = ga4q_id 
     &    +(real(conjg(lo_tot(h1,h2,h3))*v2_del2_tot(h1,h2,h3),dp)
     &    +real(conjg(v2_del2_tot(h1,h2,h3))*lo_tot(h1,h2,h3),dp))*(-one) 
     &    + (real(conjg(lo2_tot(h1,h2,h3))*v1_del2_tot(h1,h2,h3),dp)
     &     +real(conjg(v1_del2_tot(h1,h2,h3))*lo2_tot(h1,h2,h3),dp))*(one)        
!            enddo
         enddo
      enddo
      
      ga4q_id=ga4q_id*fac
 !     write(6,*) epinv2,ga4q_noid
      return
      end
