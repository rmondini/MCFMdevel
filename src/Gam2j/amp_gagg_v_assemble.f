!--- CW May 16
!----- this routine fills the various virtual amplitudes for the following process

!-----q(i1)+qb(i2)+g(i3)+g(i4)+gamma(i5)

!-----the array is 4 dim and corresponds to the helicities of q(i1), g(i3), g(i4) and gamma(i5)

!---- THE PHOTON COUPLES TO THE  q(i1),q(ib) line, no charges are added.

!---  routine returns 3 arrays:  amp_34 prop to (T^a3 T^a4)_12
!---                          :  amp_43 prop to (T^a4 T^a3)_12      
!---                          :  amp_del prop to del_34 del_12 



      subroutine amp_qqbggga_v_34_assemble(i1,i2,i3,i4,i5,za,zb
     &     ,amp_34,amp_nf)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'nf.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'ewcharge.f'
      include 'epinv2.f'
      include 'epinv.f'
      integer i1,i2,i3,i4,i5
      complex(dp) :: amp_34_lc(2,2,2,2),amp_34_slc(2,2,2,2)
     &,amp_34_nf(2,2,2,2)
      complex(dp) :: amp_43_lc(2,2,2,2),amp_43_slc(2,2,2,2)
     &,amp_nf(2,2,2,2)
      complex(dp) :: amp_del_lc(2,2,2,2),amp_del_slc(2,2,2,2)
     &,amp_del_nf(2,2,2,2)
      complex(dp) :: amp_34(2,2,2,2)
      complex(dp) :: amp_34_a51_nf(2,2,2,2)

      
      complex(dp):: amp_qqbggga_lc_gaMHV, amp_qqbggga_slc_gaMHV
      complex(dp):: amp_qqbggga_lc_gMHV, amp_qqbggga_slc_gMHV
      complex(dp):: amp_qqbggga_lc_gMHVadj, amp_qqbggga_slc_gMHVadj
      complex(dp):: amp_qqbggga_nf_gaMHV
      complex(dp):: amp_qqbggga_nf_gMHV
      complex(dp):: amp_qqbggga_nf_gMHVadj
      complex(dp):: amp_qqbggga_nf_A51_gaMHV
      integer h1,h2,h3,h4,i
      complex(dp):: phase
      complex(dp) :: amp_tree(2,2,2,2)
      real(dp) :: Qsum 
!      phase=-cone*half

      amp_34_lc(:,:,:,:)=czip
      amp_34_slc(:,:,:,:)=czip
      amp_34_nf(:,:,:,:)=czip
      Qsum=zip
      do i=1,nf
         Qsum=Qsum+Q(i)
      enddo
         
!------ T^34 pieces

!----- leading color  (prop to xn)
!====== rational amplitudes which interfere with 0 LO so dont bother to compute
      amp_34_lc(1,2,2,2)=czip
      amp_34_lc(1,1,1,1)=czip
      amp_34_lc(2,1,1,1)=czip
      amp_34_lc(2,2,2,2)=czip

  
!======= gamma MHV amplitudes 
      amp_34_lc(1,2,2,1)=amp_qqbggga_lc_gaMHV(i1,i2,i3,i4,i5,za,zb)
      amp_34_lc(2,1,1,2)=-amp_qqbggga_lc_gaMHV(i1,i2,i3,i4,i5,zb,za)
      amp_34_lc(2,2,2,1)=amp_qqbggga_lc_gaMHV(i2,i1,i4,i3,i5,za,zb)
      amp_34_lc(1,1,1,2)=-amp_qqbggga_lc_gaMHV(i2,i1,i4,i3,i5,zb,za)

      
!======== adjacent MHV
      amp_34_lc(1,2,1,2)=amp_qqbggga_lc_gMHVadj(i1,i2,i3,i4,i5,za,zb)
      amp_34_lc(2,1,2,1)=-amp_qqbggga_lc_gMHVadj(i1,i2,i3,i4,i5,zb,za)
      amp_34_lc(2,1,2,2)=amp_qqbggga_lc_gMHVadj(i2,i1,i4,i3,i5,za,zb)
      amp_34_lc(1,2,1,1)=-amp_qqbggga_lc_gMHVadj(i2,i1,i4,i3,i5,zb,za)

!====== non adjacent MHV
      
      amp_34_lc(1,1,2,2)=amp_qqbggga_lc_gMHV(i1,i2,i3,i4,i5,za,zb)
      amp_34_lc(2,2,1,1)=-amp_qqbggga_lc_gMHV(i1,i2,i3,i4,i5,zb,za)
      amp_34_lc(2,2,1,2)=amp_qqbggga_lc_gMHV(i2,i1,i4,i3,i5,za,zb)
      amp_34_lc(1,1,2,1)=-amp_qqbggga_lc_gMHV(i2,i1,i4,i3,i5,zb,za)


!====== sub leading color (prop to 1/xn)
      amp_34_slc(1,2,2,2)=czip
      amp_34_slc(1,1,1,1)=czip
      amp_34_slc(2,1,1,1)=czip
      amp_34_slc(2,2,2,2)=czip

  
!======= gamma MHV amplitudes 
      amp_34_slc(1,2,2,1)=amp_qqbggga_slc_gaMHV(i1,i2,i3,i4,i5,za,zb)
      amp_34_slc(2,1,1,2)=-amp_qqbggga_slc_gaMHV(i1,i2,i3,i4,i5,zb,za)
      amp_34_slc(2,2,2,1)=amp_qqbggga_slc_gaMHV(i2,i1,i4,i3,i5,za,zb)
      amp_34_slc(1,1,1,2)=-amp_qqbggga_slc_gaMHV(i2,i1,i4,i3,i5,zb,za)

!======== adjacent MHV
      amp_34_slc(1,2,1,2)=amp_qqbggga_slc_gMHVadj(i1,i2,i3,i4,i5,za,zb)
      amp_34_slc(2,1,2,1)=-amp_qqbggga_slc_gMHVadj(i1,i2,i3,i4,i5,zb,za)
      amp_34_slc(2,1,2,2)=amp_qqbggga_slc_gMHVadj(i2,i1,i4,i3,i5,za,zb)
      amp_34_slc(1,2,1,1)=-amp_qqbggga_slc_gMHVadj(i2,i1,i4,i3,i5,zb,za)

!====== non adjacent MHV
      amp_34_slc(1,1,2,2)=amp_qqbggga_slc_gMHV(i1,i2,i3,i4,i5,za,zb)
      amp_34_slc(2,2,1,1)=-amp_qqbggga_slc_gMHV(i1,i2,i3,i4,i5,zb,za)
      amp_34_slc(2,2,1,2)=amp_qqbggga_slc_gMHV(i2,i1,i4,i3,i5,za,zb)
      amp_34_slc(1,1,2,1)=-amp_qqbggga_slc_gMHV(i2,i1,i4,i3,i5,zb,za)


!====== nf loops
      amp_34_nf(1,2,2,2)=czip
      amp_34_nf(1,1,1,1)=czip
      amp_34_nf(2,1,1,1)=czip
      amp_34_nf(2,2,2,2)=czip

     
!======= gamma MHV amplitudes 
      amp_34_nf(1,2,2,1)=amp_qqbggga_nf_gaMHV(i1,i2,i3,i4,i5,za,zb)
      amp_34_nf(2,1,1,2)=-amp_qqbggga_nf_gaMHV(i1,i2,i3,i4,i5,zb,za)
      amp_34_nf(2,2,2,1)=amp_qqbggga_nf_gaMHV(i2,i1,i4,i3,i5,za,zb)
      amp_34_nf(1,1,1,2)=-amp_qqbggga_nf_gaMHV(i2,i1,i4,i3,i5,zb,za)

!======== adjacent MHV
      amp_34_nf(1,2,1,2)=amp_qqbggga_nf_gMHVadj(i1,i2,i3,i4,i5,za,zb)
      amp_34_nf(2,1,2,1)=-amp_qqbggga_nf_gMHVadj(i1,i2,i3,i4,i5,zb,za)
      amp_34_nf(2,1,2,2)=amp_qqbggga_nf_gMHVadj(i2,i1,i4,i3,i5,za,zb)
      amp_34_nf(1,2,1,1)=-amp_qqbggga_nf_gMHVadj(i2,i1,i4,i3,i5,zb,za)

!====== non adjacent MHV
      amp_34_nf(1,1,2,2)=amp_qqbggga_nf_gMHV(i1,i2,i3,i4,i5,za,zb)
      amp_34_nf(2,2,1,1)=-amp_qqbggga_nf_gMHV(i1,i2,i3,i4,i5,zb,za)
      amp_34_nf(2,2,1,2)=amp_qqbggga_nf_gMHV(i2,i1,i4,i3,i5,za,zb)
      amp_34_nf(1,1,2,1)=-amp_qqbggga_nf_gMHV(i2,i1,i4,i3,i5,zb,za)

!-----N_f pieces which go like the tree-coupling and get added in with
!----- this bit 
      amp_34_a51_nf(:,:,:,:)=czip
!-----photon mhv
      amp_34_a51_nf(1,2,2,1)=
     &     amp_qqbggga_nf_A51_gaMHV(i1,i2,i3,i4,i5,za,zb)
      amp_34_a51_nf(2,1,1,2)=
     &     -amp_qqbggga_nf_A51_gaMHV(i1,i2,i3,i4,i5,zb,za)
      amp_34_a51_nf(2,2,2,1)=
     &     amp_qqbggga_nf_A51_gaMHV(i2,i1,i4,i3,i5,za,zb)
      amp_34_a51_nf(1,1,1,2)=
     &     -amp_qqbggga_nf_A51_gaMHV(i2,i1,i4,i3,i5,zb,za)
      
!-----BDK style notation with xn factored out of A 61
      call qqbgg_ga_amp(i1,i2,i3,i4,i5,za,zb,amp_tree)
      
!==== note that the "nf" contribution has two pieces
!==== one which goes like the sum over all charges, (box diagrams)
!---- and one (purely rational and only for gaMHV) which goes like
!---- the LO bit. (cf 8.2 in BDK four partons)
      do h1=1,2
         do h2=1,2
            do h3=1,2
               do h4=1,2
                  amp_34(h1,h2,h3,h4)=amp_34_lc(h1,h2,h3,h4)+
     &                 one/xn**2*amp_34_slc(h1,h2,h3,h4)
     &                 +nf/xn*amp_34_a51_nf(h1,h2,h3,h4)
               enddo
            enddo
         enddo
      enddo

      amp_nf(:,:,:,:)=Qsum*amp_34_nf(:,:,:,:)

      
!       amp_34(:,:,:,:)=phase*amp_34(:,:,:,:)
      
      
   
      if(1.eq.2) then 
      write(6,*) '************ T 34 LC pieces ***********'
      do h1=1,2
         do h2=1,2
            do h3=1,2
               do h4=1,2
                  write(6,*) h1,h2,h3,h4,amp_34_lc(h1,h2,h3,h4)*im
               enddo
            enddo
         enddo
      enddo
      write(6,*) '****************************************'
      write(6,*)
      write(6,*) '************ T 34 SLC pieces ***********'
      do h1=1,2
         do h2=1,2
            do h3=1,2
               do h4=1,2
                  write(6,*) h1,h2,h3,h4,amp_34_slc(h1,h2,h3,h4)*im
               enddo
            enddo
         enddo
      enddo
      write(6,*) '****************************************'
      write(6,*)
      write(6,*) '************ T 34 nf pieces ***********'
      do h1=1,2
         do h2=1,2
            do h3=1,2
               do h4=1,2
                  write(6,*) h1,h2,h3,h4,amp_34_nf(h1,h2,h3,h4)*im
               enddo
            enddo
         enddo
      enddo
      write(6,*) '****************************************'
      write(6,*)
      write(6,*) '************ T 34 nf: A_s pieces ***********'
      do h1=1,2
         do h2=1,2
            do h3=1,2
               do h4=1,2
                  write(6,*) h1,h2,h3,h4,amp_34_A51_nf(h1,h2,h3,h4)*im
               enddo
            enddo
         enddo
      enddo
      write(6,*) '****************************************'
      write(6,*)
      stop
   
      endif
      
     


      return
      end
      subroutine amp_qqbggga_v_43_assemble(i1,i2,i3,i4,i5,za,zb
     &     ,amp_43)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'nf.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer i1,i2,i3,i4,i5
      complex(dp) :: amp_43(2,2,2,2)
      complex(dp) :: amp_34_lc(2,2,2,2),amp_34_slc(2,2,2,2)
     &,amp_34_nf(2,2,2,2)
      complex(dp) :: amp_43_lc(2,2,2,2),amp_43_slc(2,2,2,2)
     &,amp_43_nf(2,2,2,2)
      complex(dp) :: amp_del_lc(2,2,2,2),amp_del_slc(2,2,2,2)
     &,amp_del_nf(2,2,2,2)
      complex temp(2,2,2,2),temp_slc(2,2,2,2),temp_nf(2,2,2,2)
      integer h1,h2,h3,h4


      call amp_qqbggga_v_34_assemble(i1,i2,i4,i3,i5,za,zb
     &     ,amp_43)

!---- helicities of swap 
!      temp(:,:,:,:)= amp_43(:,:,:,:)
    
      
!      do h1=1,2
!         do h2=1,2
!            do h3=1,2
!               do h4=1,2
!                  amp_43(h1,h2,h3,h4)=temp(h1,h3,h2,h4)
!                enddo
!            enddo
!         enddo
!      enddo
                  
                  
      if(1.eq.2) then 
      write(6,*) '************ T 43 LC pieces ***********'
      do h1=1,2
         do h2=1,2
            do h3=1,2
               do h4=1,2
                  write(6,*) h1,h2,h3,h4,amp_43_lc(h1,h2,h3,h4)*im
               enddo
            enddo
         enddo
      enddo
      write(6,*) '****************************************'
      write(6,*)
      write(6,*) '************ T 43 SLC pieces ***********'
      do h1=1,2
         do h2=1,2
            do h3=1,2
               do h4=1,2
                  write(6,*) h1,h2,h3,h4,amp_43_slc(h1,h2,h3,h4)*im
               enddo
            enddo
         enddo
      enddo
      write(6,*) '****************************************'
      write(6,*)
      write(6,*) '************ T 43 nf pieces ***********'
      do h1=1,2
         do h2=1,2
            do h3=1,2
               do h4=1,2
                  write(6,*) h1,h2,h3,h4,amp_43_nf(h1,h2,h3,h4)*im
               enddo
            enddo
         enddo
      enddo
      write(6,*) '****************************************'
      write(6,*)
      endif
     


      return
      end



      subroutine amp_qqbggga_v_del34_assemble(i1,i2,i3,i4,i5,za,zb
     &     ,amp_del34)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'nf.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer i1,i2,i3,i4,i5
      complex(dp) :: amp_del34(2,2,2,2)
    
      complex(dp):: amp_qqbggga_del34_gaMHV
      complex(dp):: amp_qqbggga_del34_gMHVadj
      complex(dp):: amp_qqbggga_del34_gMHV
      integer h1,h2,h3,h4


      amp_del34(:,:,:,:)=czip
     
      
!------ del34 pieces

!----- leading color  (prop to xn)
!====== rational amplitudes which interfere with 0 LO so dont bother to compute
      amp_del34(1,2,2,2)=czip
      amp_del34(1,1,1,1)=czip
      amp_del34(2,1,1,1)=czip
      amp_del34(2,2,2,2)=czip

  
!======= gamma MHV amplitudes 
      amp_del34(1,2,2,1)=amp_qqbggga_del34_gaMHV(i1,i2,i3,i4,i5,za,zb)
      amp_del34(2,1,1,2)=-amp_qqbggga_del34_gaMHV(i1,i2,i3,i4,i5,zb,za)
      amp_del34(2,2,2,1)=amp_qqbggga_del34_gaMHV(i2,i1,i4,i3,i5,za,zb)
      amp_del34(1,1,1,2)=-amp_qqbggga_del34_gaMHV(i2,i1,i4,i3,i5,zb,za)

!======== adjacent MHV
      amp_del34(1,1,2,2)=amp_qqbggga_del34_gMHV(i1,i2,i4,i3,i5,za,zb)
      amp_del34(2,2,1,1)=-amp_qqbggga_del34_gMHV(i1,i2,i4,i3,i5,zb,za)
      amp_del34(2,2,1,2)=amp_qqbggga_del34_gMHV(i2,i1,i3,i4,i5,za,zb)
      amp_del34(1,1,2,1)=-amp_qqbggga_del34_gMHV(i2,i1,i3,i4,i5,zb,za)

!====== non adjacent MHV
      amp_del34(1,2,1,2)=amp_qqbggga_del34_gMHV(i1,i2,i3,i4,i5,za,zb)
      amp_del34(2,1,2,1)=-amp_qqbggga_del34_gMHV(i1,i2,i3,i4,i5,zb,za)
      amp_del34(2,1,2,2)=amp_qqbggga_del34_gMHV(i2,i1,i4,i3,i5,za,zb)
      amp_del34(1,2,1,1)=-amp_qqbggga_del34_gMHV(i2,i1,i4,i3,i5,zb,za)

      if(1.eq.2) then 
      write(6,*) '************ del 34 pieces ***********'
      do h1=1,2
         do h2=1,2
            do h3=1,2
               do h4=1,2
                  write(6,*) h1,h2,h3,h4,amp_del34(h1,h2,h3,h4)*im
               enddo
            enddo
         enddo
      enddo
      write(6,*) '****************************************'
      endif


      return
      end
