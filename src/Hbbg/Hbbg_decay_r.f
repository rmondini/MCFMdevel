      subroutine Hbbg_decay_r(p,i1,i2,i3,i4,msqr)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'ewcouple.f'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'hbbparams.f'
      include 'nf.f'
      real(dp) :: p(mxpart,4),msqr
      complex(dp) :: ampsgg(0:2,2,2,2),ampsQQ(2,2,2)
      integer i1,i2,i3,i4
      real(dp) :: facQQ,facgg
      real(dp) :: msqQQ,msqgg,msqQQint
      integer imax,h1,h2,h3
      real(dp) babisHbbQQ_msq,babisHbbgg_msq,babis,babisgg
      logical nfonly
      common/nfonly/nfonly
!$omp threadprivate(/nfonly/)

!-----debug test against KC
!       include 'phasespace_KC_r.f'
!        call spinorz(4,p,za,zb)      
!-----end debug

      imax=max(i1,i2,i3,i4)
      call spinoru(imax,p,za,zb)
      call Hbbg_realamps(i1,i2,i3,i4,za,zb,ampsgg,ampsQQ)
      facQQ=V*mb_eff**2*gsq**2*gwsq/wmass**2/four
      facgg=V*xn*mb_eff**2*gsq**2*gwsq/wmass**2/eight

!----- H-> bbgg contribution 
      msqgg=zip
      do h1=1,2
         do h2=1,2
            do h3=1,2
               msqgg=msqgg
     &            +real(ampsgg(1,h1,h2,h3)*conjg(ampsgg(1,h1,h2,h3)),dp)
     &            +real(ampsgg(2,h1,h2,h3)*conjg(ampsgg(2,h1,h2,h3)),dp)
     &-1._dp/xn**2*real(ampsgg(0,h1,h2,h3)*conjg(ampsgg(0,h1,h2,h3)),dp)
            enddo
         enddo
      enddo

      msqgg=msqgg*facgg
!      if(nfonly) then
!         msqgg=zip
!      endif
      
!----- H->bb QQ contribution (x nf flavors)
      msqQQ=zip
      do h1=1,2
         do h2=1,2
            msqQQ=msqQQ
     &            +real(ampsQQ(1,h1,h2)*conjg(ampsQQ(1,h1,h2)),dp)
            enddo
         enddo

         msqQQint=zip
         do h1=1,2
            do h2=1,2
               msqQQint=msqQQint
     &            +(ampsQQ(2,h1,h2)*conjg(ampsQQ(1,h1,h2))
     &            +ampsQQ(1,h1,h2)*conjg(ampsQQ(2,h1,h2)))
            enddo
         enddo

!----- code below is for checking with literature MEs         
!      babis=babisHbbQQ_msq(i1,i2,i3,i4)
!      babisgg=babisHbbgg_msq(i1,i2,i3,i4)
!      write(6,*) 'b check 4Q'
!      write(6,*) babis,msqQQ*facQQ,babis/(msqQQ*facQQ)
!      write(6,*) '****'
!      write(6,*) 'b check gg'
!      write(6,*) babisgg,msqgg,babisgg/(msqgg)
!      write(6,*) '****'
!      pause


      msqQQ=msqQQ*facQQ*nf+msqQQint*facQQ/xn
!---- no nf option
!      msqQQ=zip
  
      
      msqr=msqQQ+msqgg
      return
      end
