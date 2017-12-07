!-----routine to fill the (LO) q(i1)+qb(i2)=>g(i3)+g(i4)+gamma(i5) msq
!-----using amplitudes (to setup similar virtual routine)


      subroutine qqbgg_ga(i1,i2,i3,i4,i5,za,zb,qqb_gagg)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'mmsq_cs_ga.f'

      integer i1,i2,i3,i4,i5
      real(dp):: qqb_gagg
      complex(dp) :: amp34(2,2,2,2),amp43(2,2,2,2),
     &     ampQED(2,2,2,2),
     &     qqbgg_gaMHV,qqbgg_gMHVadj,qqbgg_gMHV
      integer h1,h2,h3,h4
      real(dp) test,BigCgam


      
      qqb_gagg=zip
      call qqbgg_ga_amp(i1,i2,i3,i4,i5,za,zb,amp34)
      call qqbgg_ga_amp(i1,i2,i4,i3,i5,za,zb,amp43)

      do h1=1,2
         do h2=1,2
            do h3=1,2
               do h4=1,2
                  ampQED(h1,h2,h3,h4)=amp34(h1,h2,h3,h4)+
     &                 amp43(h1,h3,h2,h4)
               enddo
            enddo
         enddo
      enddo

!----- sum over polarizations
      mmsq_cs_ga(:)=zip
      do h1=1,2
         do h2=1,2
            do h3=1,2
               do h4=1,2
       mmsq_cs_ga(1)= mmsq_cs_ga(1) +real(amp34(h1,h2,h3,h4)
     &                  *conjg(amp34(h1,h2,h3,h4)),dp)
       mmsq_cs_ga(2)= mmsq_cs_ga(2) +real(amp43(h1,h2,h3,h4)
     &                  *conjg(amp43(h1,h2,h3,h4)),dp)
       mmsq_cs_ga(0)= mmsq_cs_ga(0) +real(ampQED(h1,h2,h3,h4)
     &                  *conjg(ampQED(h1,h2,h3,h4)),dp)
      enddo
      enddo
      enddo
      enddo

    
      mmsq_cs_ga(0)=-1._dp/xn**2*mmsq_cs_ga(0)
      qqb_gagg=mmsq_cs_ga(1)+mmsq_cs_ga(2)+mmsq_cs_ga(0)
!      write(6,*) 'lord: amp 34 bit',mmsq_cs_ga(1)
!      write(6,*) 'lord: amp 43 bit',mmsq_cs_ga(2)
!      write(6,*) 'lord: QED bit',mmsq_cs_ga(0)
!      write(6,*) 'total ',qqb_gagg
      return
      end
      
      subroutine qqbgg_ga_amp(i1,i2,i3,i4,i5,za,zb,amp)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'

      integer i1,i2,i3,i4,i5
      real(dp):: qqb_gagg
      complex(dp) :: amp(2,2,2,2),qqbgg_gaMHV,qqbgg_gMHVadj,qqbgg_gMHV
      complex(dp) :: fullamp(2,2,2,2),sign
      integer h1,h2,h3,h4

 !     qqb_gagg=zip

      amp(2,2,2,2)=czip
      amp(1,1,1,1)=czip
      amp(1,2,2,2)=czip
      amp(2,1,1,1)=czip

      amp(1,2,2,1)=qqbgg_gaMHV(i1,i2,i3,i4,i5,za,zb)
      amp(1,2,1,2)=qqbgg_gMHVadj(i1,i2,i3,i4,i5,za,zb)
      amp(1,1,2,2)=qqbgg_gMHV(i1,i2,i3,i4,i5,za,zb)

      amp(1,1,1,2)=-qqbgg_gaMHV(i2,i1,i4,i3,i5,zb,za)
      amp(1,2,1,1)=-qqbgg_gMHVadj(i2,i1,i4,i3,i5,zb,za)
      amp(1,1,2,1)=-qqbgg_gMHV(i2,i1,i4,i3,i5,zb,za)
      
!      amp(2,1,1,2)=-qqbgg_gaMHV(i1,i2,i3,i4,i5,zb,za)
!      amp(2,1,2,1)=-qqbgg_gMHVadj(i1,i2,i3,i4,i5,zb,za)
!      amp(2,2,1,1)=-qqbgg_gMHV(i1,i2,i3,i4,i5,zb,za)

!      amp(2,2,2,1)=qqbgg_gaMHV(i2,i1,i4,i3,i5,za,zb)
!      amp(2,1,2,2)=qqbgg_gMHVadj(i2,i1,i4,i3,i5,za,zb)
!      amp(2,2,1,2)=qqbgg_gMHV(i2,i1,i4,i3,i5,za,zb)

!      fullamp(:,:,:,:)=amp(:,:,:,:)

! for speed, use trick to get remaining amplitudes; fix sign by explicit check of crossing
      if (((i1<3).and.(i2<3)) .or. ((i1>2).and.(i2>2))) then
        sign=+cone
      else
        sign=-cone
      endif
      
      do h2=1,2
      do h3=1,2
      do h4=1,2
      amp(2,h2,h3,h4)=sign*conjg(amp(1,3-h2,3-h3,3-h4))
      enddo
      enddo
      enddo
      
!      qqb_gagg=zip
!      do h1=1,2
!         do h2=1,2
!            do h3=1,2
!               do h4=1,2
!                  write(6,*) h1,h2,h3,h4,amp(h1,h2,h3,h4),fullamp(h1,h2,h3,h4),amp(h1,h2,h3,h4)/fullamp(h1,h2,h3,h4)
!               enddo
!            enddo
!         enddo
!      enddo
!      pause

!      qqb_gagg=zip
!      do h1=1,2
!         do h2=1,2
!            do h3=1,2
!               do h4=1,2
!                  write(6,*) h1,h2,h3,h4,amp(h1,h2,h3,h4)
!                  qqb_gagg=qqb_gagg
!     &                +real(amp(h1,h2,h3,h4)*conjg(amp(h1,h2,h3,h4)),dp)
!               enddo
!            enddo
!         enddo
!      enddo
!      pause
      return
      end

      function qqbgg_gaMHV(i1,i2,i3,i4,i5,za,zb)
      implicit none
      include 'types.f'
      complex(dp) qqbgg_gaMHV
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      integer i1,i2,i3,i4,i5

      qqbgg_gaMHV=za(i1,i5)**2/(za(i2,i3)*za(i3,i4)*za(i4,i1))
      return
      end


      function qqbgg_gMHV(i1,i2,i3,i4,i5,za,zb)
      implicit none
      include 'types.f'
      complex(dp) qqbgg_gMHV
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      integer i1,i2,i3,i4,i5

      qqbgg_gMHV=-(za(i1,i3)**3
     &     /(za(i1,i4)*za(i1,i5)*za(i2,i5)*za(i3,i4)))
      return
      end

      function qqbgg_gMHVadj(i1,i2,i3,i4,i5,za,zb)
      implicit none
      include 'types.f'
      complex(dp)  qqbgg_gMHVadj
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      integer i1,i2,i3,i4,i5

      qqbgg_gMHVadj=(za(i4,i1)**2*za(i4,i2))
     & /(za(i1,i5)*za(i2,i3)*za(i2,i5)*za(i3,i4))
      return
      end
