!--- CW May 16
!----- this routine fills the various virtual amplitudes for the following process

!-----q(i1)+qb(i2)+Q(i3)+QB(i4)+gamma(i5)

!-----the array is 3 dim and corresponds to the helicities of q(i1), Q(i3) and gamma(i5)

!---- THE PHOTON COUPLES TO THE  q(i1),q(ib) line, no charges are added.

!---  routine returns 2 arrays:  lc prop to del_1  = delta(i1,i3)*delta(i2,i4)
!---                          :  slc prop to del_2 = delta(i1,i2)*delta(i3,i4)


      subroutine amp_qqbQQbga_v_assemble(i1,i2,i3,i4,i5,za,zb
     &     ,amp_del1,amp_del2)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'nf.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer i1,i2,i3,i4,i5
      complex(dp) :: amp_del1(2,2,2),amp_del2(2,2,2)
      complex(dp) :: amp_del1_lc(2,2,2),amp_del1_slc(2,2,2)
     &     ,amp_del1_nf(2,2,2)
      complex(dp) :: amp_del2_lc(2,2,2),amp_del2_slc(2,2,2)
     &     ,amp_del2_nf(2,2,2)

      complex(dp):: amp_qqbQQbga_del1_lc,amp_qqbQQbga_del1_slc
      complex(dp):: amp_qqbQQbga_del2_lc,amp_qqbQQbga_del2_slc
      complex(dp):: amp_qqbQQbga_mhvalt_del1_lc
     &     ,amp_qqbQQbga_mhvalt_del1_slc
      complex(dp):: amp_qqbQQbga_mhvalt_del2_lc
     &     ,amp_qqbQQbga_mhvalt_del2_slc
      complex(dp):: amp_qqbQQbga_mhvalt_del2_nf
     &     ,amp_qqbQQbga_del2_nf
      complex(dp):: amp_qqbQQbga_mhvalt_del1_nf
     &     ,amp_qqbQQbga_del1_nf
      integer h1,h2,h3

      
!------lc del1 pieces

      amp_del1_lc(1,1,2) = amp_qqbQQbga_del1_lc(i1,i2,i3,i4,i5,za,zb)
      amp_del1_lc(1,2,2) =
     &     amp_qqbQQbga_mhvalt_del1_lc(i1,i2,i3,i4,i5,za,zb)

      amp_del1_lc(2,2,1) = -amp_qqbQQbga_del1_lc(i1,i2,i3,i4,i5,zb,za)
      amp_del1_lc(2,1,1) =
     &     -amp_qqbQQbga_mhvalt_del1_lc(i1,i2,i3,i4,i5,zb,za)


      amp_del1_lc(1,1,1) = amp_qqbQQbga_del1_lc(i2,i1,i4,i3,i5,zb,za)
      amp_del1_lc(2,2,2) = -amp_qqbQQbga_del1_lc(i2,i1,i4,i3,i5,za,zb)

      amp_del1_lc(1,2,1) =
     &     amp_qqbQQbga_mhvalt_del1_lc(i2,i1,i4,i3,i5,zb,za)
      amp_del1_lc(2,1,2) =
     &     -amp_qqbQQbga_mhvalt_del1_lc(i2,i1,i4,i3,i5,za,zb)
    

!------slc del1 pieces

      amp_del1_slc(1,1,2) = amp_qqbQQbga_del1_slc(i1,i2,i3,i4,i5,za,zb)
      amp_del1_slc(1,2,2) =
     &     amp_qqbQQbga_mhvalt_del1_slc(i1,i2,i3,i4,i5,za,zb)

      amp_del1_slc(2,2,1) = -amp_qqbQQbga_del1_slc(i1,i2,i3,i4,i5,zb,za)
      amp_del1_slc(2,1,1) =
     &     -amp_qqbQQbga_mhvalt_del1_slc(i1,i2,i3,i4,i5,zb,za)


      amp_del1_slc(1,1,1) = amp_qqbQQbga_del1_slc(i2,i1,i4,i3,i5,zb,za)
      amp_del1_slc(2,2,2) = -amp_qqbQQbga_del1_slc(i2,i1,i4,i3,i5,za,zb)

      amp_del1_slc(1,2,1) =
     &     amp_qqbQQbga_mhvalt_del1_slc(i2,i1,i4,i3,i5,zb,za)
      amp_del1_slc(2,1,2) =
     &     -amp_qqbQQbga_mhvalt_del1_slc(i2,i1,i4,i3,i5,za,zb)
    

!------nf del1 pieces

      amp_del1_nf(1,1,2) = amp_qqbQQbga_del1_nf(i1,i2,i3,i4,i5,za,zb)
      amp_del1_nf(1,2,2) =
     &     amp_qqbQQbga_mhvalt_del1_nf(i1,i2,i3,i4,i5,za,zb)

      amp_del1_nf(2,2,1) = -amp_qqbQQbga_del1_nf(i1,i2,i3,i4,i5,zb,za)
      amp_del1_nf(2,1,1) =
     &     -amp_qqbQQbga_mhvalt_del1_nf(i1,i2,i3,i4,i5,zb,za)


      amp_del1_nf(1,1,1) = amp_qqbQQbga_del1_nf(i2,i1,i4,i3,i5,zb,za)
      amp_del1_nf(2,2,2) = -amp_qqbQQbga_del1_nf(i2,i1,i4,i3,i5,za,zb)

      amp_del1_nf(1,2,1) =
     &     amp_qqbQQbga_mhvalt_del1_nf(i2,i1,i4,i3,i5,zb,za)
      amp_del1_nf(2,1,2) =
     &     -amp_qqbQQbga_mhvalt_del1_nf(i2,i1,i4,i3,i5,za,zb)
    

            
!------lc del2 pieces

      amp_del2_lc(1,1,2) = amp_qqbQQbga_del2_lc(i1,i2,i3,i4,i5,za,zb)
      amp_del2_lc(1,2,2) =
     &     amp_qqbQQbga_mhvalt_del2_lc(i1,i2,i3,i4,i5,za,zb)

      amp_del2_lc(2,2,1) = -amp_qqbQQbga_del2_lc(i1,i2,i3,i4,i5,zb,za)
      amp_del2_lc(2,1,1) =
     &     -amp_qqbQQbga_mhvalt_del2_lc(i1,i2,i3,i4,i5,zb,za)


      amp_del2_lc(1,1,1) = amp_qqbQQbga_del2_lc(i2,i1,i4,i3,i5,zb,za)
      amp_del2_lc(2,2,2) = -amp_qqbQQbga_del2_lc(i2,i1,i4,i3,i5,za,zb)

      amp_del2_lc(1,2,1) =
     &     amp_qqbQQbga_mhvalt_del2_lc(i2,i1,i4,i3,i5,zb,za)
      amp_del2_lc(2,1,2) =
     &     -amp_qqbQQbga_mhvalt_del2_lc(i2,i1,i4,i3,i5,za,zb)
    

!------slc del2 pieces

      amp_del2_slc(1,1,2) = amp_qqbQQbga_del2_slc(i1,i2,i3,i4,i5,za,zb)
      amp_del2_slc(1,2,2) =
     &     amp_qqbQQbga_mhvalt_del2_slc(i1,i2,i3,i4,i5,za,zb)

      amp_del2_slc(2,2,1) = -amp_qqbQQbga_del2_slc(i1,i2,i3,i4,i5,zb,za)
      amp_del2_slc(2,1,1) =
     &     -amp_qqbQQbga_mhvalt_del2_slc(i1,i2,i3,i4,i5,zb,za)


      amp_del2_slc(1,1,1) = amp_qqbQQbga_del2_slc(i2,i1,i4,i3,i5,zb,za)
      amp_del2_slc(2,2,2) = -amp_qqbQQbga_del2_slc(i2,i1,i4,i3,i5,za,zb)

      amp_del2_slc(1,2,1) =
     &     amp_qqbQQbga_mhvalt_del2_slc(i2,i1,i4,i3,i5,zb,za)
      amp_del2_slc(2,1,2) =
     &     -amp_qqbQQbga_mhvalt_del2_slc(i2,i1,i4,i3,i5,za,zb)
    

!------nf del2 pieces

      amp_del2_nf(1,1,2) = amp_qqbQQbga_del2_nf(i1,i2,i3,i4,i5,za,zb)
      amp_del2_nf(1,2,2) =
     &     amp_qqbQQbga_mhvalt_del2_nf(i1,i2,i3,i4,i5,za,zb)

      amp_del2_nf(2,2,1) = -amp_qqbQQbga_del2_nf(i1,i2,i3,i4,i5,zb,za)
      amp_del2_nf(2,1,1) =
     &     -amp_qqbQQbga_mhvalt_del2_nf(i1,i2,i3,i4,i5,zb,za)


      amp_del2_nf(1,1,1) = amp_qqbQQbga_del2_nf(i2,i1,i4,i3,i5,zb,za)
      amp_del2_nf(2,2,2) = -amp_qqbQQbga_del2_nf(i2,i1,i4,i3,i5,za,zb)

      amp_del2_nf(1,2,1) =
     &     amp_qqbQQbga_mhvalt_del2_nf(i2,i1,i4,i3,i5,zb,za)
      amp_del2_nf(2,1,2) =
     &     -amp_qqbQQbga_mhvalt_del2_nf(i2,i1,i4,i3,i5,za,zb)


!=======final assemble
      amp_del1(:,:,:)=czip
      amp_del2(:,:,:)=czip

      do h1=1,2
         do h2=1,2
            do h3=1,2
               amp_del1(h1,h2,h3) = xn*amp_del1_lc(h1,h2,h3)
     &      +1._dp/xn*amp_del1_slc(h1,h2,h3)+nf/xn*amp_del1_nf(h1,h2,h3)
               amp_del2(h1,h2,h3) = +amp_del2_lc(h1,h2,h3)
     &              +1._dp/xn**2*amp_del2_slc(h1,h2,h3)
     &              +nf/xn**2*amp_del2_nf(h1,h2,h3)
            enddo
         enddo
      enddo

        if(1.eq.1) then
         write(6,*) '*********** del 1 pieces ****************'
         write(6,*) '* all pieces                           *'
         do h1=1,2
            do h2=1,2
               do h3=1,2
                  write(6,*) h1,h2,h3,amp_del1(h1,h2,h3)*im
               enddo
            enddo
         enddo
 !        write(6,*) '*********** del 2 pieces ****************'
 !        write(6,*) '* all pieces                           *'
 !        do h1=1,2
 !           do h2=1,2
 !              do h3=1,2
 !                 write(6,*) h1,h2,h3,amp_del2(h1,h2,h3)*im
 !              enddo
 !           enddo
 !        enddo
      endif

!------test with the Kirill Check
      if(1.eq.2) then
         write(6,*) '*********** del 1 pieces ****************'
         write(6,*) '* xn lc pieces                          *'
         do h1=1,2
            do h2=1,2
               do h3=1,2
                  write(6,*) h1,h2,h3,amp_del1_lc(h1,h2,h3)*im
               enddo
            enddo
         enddo
         write(6,*) '*****************************************'
         write(6,*)
         write(6,*) '*********** del 1 pieces ****************'
         write(6,*) '* 1/xn slc pieces                          *'
         do h1=1,2
            do h2=1,2
               do h3=1,2
                  write(6,*) h1,h2,h3,amp_del1_slc(h1,h2,h3)*im
               enddo
            enddo
         enddo
         write(6,*) '*****************************************'
         write(6,*)
         write(6,*) '*********** del 1 pieces ****************'
         write(6,*) '* nf  pieces                          *'
         do h1=1,2
            do h2=1,2
               do h3=1,2
                  write(6,*) h1,h2,h3,amp_del1_nf(h1,h2,h3)*im
               enddo
            enddo
         enddo
         write(6,*) '*****************************************'
         write(6,*)
         

                  write(6,*) '*********** del 2 pieces ****************'
         write(6,*) '* 1 lc pieces                          *'
         do h1=1,2
            do h2=1,2
               do h3=1,2
                  write(6,*) h1,h2,h3,amp_del2_lc(h1,h2,h3)*im
               enddo
            enddo
         enddo
         write(6,*) '*****************************************'
         write(6,*)
         write(6,*) '*********** del 2 pieces ****************'
         write(6,*) '* 1/xn**2 slc pieces                          *'
         do h1=1,2
            do h2=1,2
               do h3=1,2
                  write(6,*) h1,h2,h3,amp_del2_slc(h1,h2,h3)*im
               enddo
            enddo
         enddo
         write(6,*) '*****************************************'
         write(6,*)
         write(6,*) '*********** del 2 pieces ****************'
         write(6,*) '* nf/xn  pieces                          *'
         do h1=1,2
            do h2=1,2
               do h3=1,2
                  write(6,*) h1,h2,h3,amp_del2_nf(h1,h2,h3)*im
               enddo
            enddo
         enddo
         write(6,*) '*****************************************'
         write(6,*)
         stop
      endif

      return
      end
