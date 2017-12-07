!-----routine to fill the (LO) q(i1)+qb(i2)=>g(i3)+g(i4)+gamma(i5) msq
!-----using amplitudes (to setup similar virtual routine)


      subroutine qqbgg_ga_v(i1,i2,i3,i4,i5,za,zb,qqb_gagg,qqb_gagg_nf)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'mmsq_cs_ga.f'
      include 'scale.f'
      include 'epinv.f'
      integer i1,i2,i3,i4,i5
      real(dp):: qqb_gagg,qqb_gagg_nf
      complex(dp) :: amp_34_v(2,2,2,2),amp_43_v(2,2,2,2),
     &     amp_del34(2,2,2,2), amp_del43(2,2,2,2)
      complex(dp) :: amp34(2,2,2,2),amp43(2,2,2,2),
     &     ampQED(2,2,2,2),
     &     qqbgg_gaMHV,qqbgg_gMHVadj,qqbgg_gMHV
      complex(dp) :: amp_comb_1(2,2,2,2),amp_comb_2(2,2,2,2)
      complex(dp) :: amp34_nf(2,2,2,2),amp43_nf(2,2,2,2)
      integer h1,h2,h3,h4
      real(dp) test,BigCgam
      complex(dp) phase

      phase=-cone*half
      
      qqb_gagg=zip
!---- LO amplitudes
      call qqbgg_ga_amp(i1,i2,i3,i4,i5,za,zb,amp34)
      call qqbgg_ga_amp(i1,i2,i4,i3,i5,za,zb,amp43)
!----- virtual amplitudes 
      call amp_qqbggga_v_34_assemble(i1,i2,i3,i4,i5,za,zb,amp_34_v,amp34_nf)
      call amp_qqbggga_v_34_assemble(i1,i2,i4,i3,i5,za,zb,amp_43_v,amp43_nf)
!      call amp_qqbggga_v_43_assemble(i1,i2,i3,i4,i5,za,zb,amp_43_v)
      call amp_qqbggga_v_del34_assemble(i1,i2,i3,i4,i5,za,zb,amp_del34)
      call amp_qqbggga_v_del34_assemble(i1,i2,i4,i3,i5,za,zb,amp_del43)
!      call qqbgg_ga_amp(i1,i2,i3,i4,i5,za,zb,amp_34_v)
!      call qqbgg_ga_amp(i1,i2,i4,i3,i5,za,zb,amp_43_v)

       amp_34_v(:,:,:,:)=amp_34_v(:,:,:,:)*phase
       amp_43_v(:,:,:,:)=amp_43_v(:,:,:,:)*phase
       amp_del34(:,:,:,:)=amp_del34(:,:,:,:)*phase
       amp34_nf(:,:,:,:)=amp34_nf(:,:,:,:)*phase*(xn-four/xn)
       amp43_nf(:,:,:,:)=amp43_nf(:,:,:,:)*phase*(xn-four/xn)
       
 !----- sum over polarizations
       qqb_gagg=zip
       qqb_gagg_nf=zip
       do h1=1,2
          do h2=1,2
             do h3=1,2
                do h4=1,2
                   amp_comb_1(h1,h2,h3,h4)=V*amp_34_v(h1,h2,h3,h4)
     &                 -amp_43_v(h1,h3,h2,h4)+amp_del34(h1,h2,h3,h4)
                  amp_comb_2(h1,h3,h2,h4)=V*amp_43_v(h1,h3,h2,h4)
     &                 -amp_34_v(h1,h2,h3,h4)+amp_del34(h1,h2,h3,h4)
               enddo
            enddo
         enddo
      enddo
      do h1=1,2
         do h2=1,2
            do h3=1,2
               do h4=1,2
      qqb_gagg=qqb_gagg    
     & +2._dp*real(conjg(amp34(h1,h2,h3,h4))*amp_comb_1(h1,h2,h3,h4),dp)
     & +2._dp*real(conjg(amp43(h1,h2,h3,h4))*amp_comb_2(h1,h2,h3,h4),dp)
      qqb_gagg_nf=qqb_gagg_nf    
     & +2._dp*real(conjg(amp34(h1,h2,h3,h4))*amp34_nf(h1,h2,h3,h4),dp)
     & +2._dp*real(conjg(amp43(h1,h2,h3,h4))*amp43_nf(h1,h2,h3,h4),dp)
                 enddo
            enddo
         enddo
      enddo
      return
      end
