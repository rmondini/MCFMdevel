
      function n21(p,psoft)
!===== p is array of parton momenta (so before soft drop)
!===== psoft is array of momenta after soft drop
!===== firstjet is the position of the first jet momentum in the array psoft
!===== nj represents the number of jets in psoft
!===== so jet momenta are located in psoft at positions firstjet,...,firstjet+nj-1
      implicit none
      include 'types.f'
      real(dp):: n21
c----calculate variable N_2^1 for boosted H(bb) CMS analysis (arxiv 1709.05543)
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'softdrop_var.f' 
      real(dp):: p(mxpart,4),psoft(mxpart,4)
      real(dp) :: e23corr,e2corr
     
      if(nplj.lt.3) then
         n21=0._dp
      else
         n21=e23corr(p)/(e2corr(p)**2)
      endif
      
      return
      end




      function zpt(i,p)
      implicit none
      include 'types.f'
      real(dp):: zpt
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'softdrop_var.f'
      real(dp):: p(mxpart,4)
      integer :: i,j
      real(dp) pt,ptsum
    
      ptsum=0._dp

      do j=1,nplj
         ptsum=ptsum+pt(pinSDj(idlptjet,j),p)
      enddo

      zpt=pt(pinSDj(idlptjet,i),p)/ptsum

      return
      end




      function e2corr(p)
      implicit none
      include 'types.f'
      real(dp):: e2corr
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      real(dp):: p(mxpart,4)
      integer :: i,j
      real(dp):: zpt,r
      include 'softdrop_var.f'

      e2corr=0._dp

      do i=1,nplj-1
         do j=i+1,nplj
            e2corr=e2corr+zpt(i,p)*zpt(j,p)*r(p,pinSDj(idlptjet,i),pinSDj(idlptjet,j))
         enddo
      enddo
      
      return
      end




      function e23corr(p)
      implicit none
      include 'types.f'
      real(dp):: e23corr
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'softdrop_var.f'
      real(dp):: p(mxpart,4)
      integer :: i,j,k
      real(dp):: zpt,e3rmin12

      e23corr=0._dp

      do i=1,nplj-2
         do j=i+1,nplj-1
            do k=j+1,nplj
               e23corr=e23corr+zpt(i,p)*zpt(j,p)*zpt(k,p)*e3rmin12(p,i,j,k)
            enddo
         enddo
      enddo

!     write(6,*) e23corr
!   write(6,*) e3rmin12(p,1,2,3)
!        pause
      
      return
      end



      
      function e3rmin12(p,ii,jj,kk)
      implicit none
      include 'types.f'
      real(dp):: e3rmin12
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'softdrop_var.f'
      real(dp):: p(mxpart,4)
      integer :: i,ii,jj,kk
      real(dp):: r
      real(dp) e3rmin1val,e3rmin2val,e3rmin3val,rvec(3)
   
      rvec(1)=r(p,pinSDj(idlptjet,ii),pinSDj(idlptjet,jj))
      rvec(2)=r(p,pinSDj(idlptjet,ii),pinSDj(idlptjet,kk))
      rvec(3)=r(p,pinSDj(idlptjet,jj),pinSDj(idlptjet,kk))

      e3rmin1val=minval(rvec)
      e3rmin3val=maxval(rvec)

      do i=1,3
         if((rvec(i).gt.e3rmin1val) .and. (rvec(i).lt.e3rmin3val)) e3rmin2val=rvec(i)
      enddo

      e3rmin12=e3rmin1val*e3rmin2val
!----- uncomment this to do 3e1 
!      e3rmin12=rvec(1)*rvec(2)*rvec(3)
      
      return
      end




