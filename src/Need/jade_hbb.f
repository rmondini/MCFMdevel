!---  test routine for H=>bb
!---- this routine clusters particles using the jade algorhtim
!---- comparing y_ij = (pi+pj)^2  to y_cut
!---- adding the particles four-momenta if true

      subroutine jade_hbb(p,ycut,jets,pjet,ip,Qsq)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'masses.f' 
      real(dp) :: p(mxpart,4)
      integer npart,jets,i,j,ip
      real(dp) :: ycut,yij
      real(dp) :: phad(mxpart,4),pjet(mxpart,4)
      integer maxjet,iter
      real(dp) :: temp(4),yij_det,yijmin,Qsq,yij_dur   
      integer nmin1,nmin2,isub
      logical is_hadronic
      
!     Qsq=125._dp**2
      pjet(:,:)=0._dp
      phad(:,:)=0._dp
!      maxjet=0
      jets=0
!      ycut=0.01_dp

      phad(:,:)=p(:,:)
      maxjet=ip
!---- first get hadrons 
!      write(6,*) npart+2-isub
!      do i=3,npart+2-isub
!         if(is_hadronic(i)) then
!            maxjet=maxjet+1
!            phad(maxjet,:)=p(i,:)
!         endif
!      enddo

 !     write(6,*) maxjet

!      if(maxjet >= 3) then 
!      write(6,*) maxjet
!      call writeout(phad)
!      pause
!      endif
      iter=0
!-----begin clustering
 1    iter=iter+1

      yijmin=zip
!---- find min y ij
      do i=iter,maxjet
         do j=i+1,maxjet
            yij=yij_det(phad,i,j,Qsq)          
            if((i==iter).and.(j==i+1)) then
               yijmin=yij
               nmin1=i
               nmin2=j
            elseif(yij < yijmin) then
               yijmin=yij
               nmin1=i
               nmin2=j
            endif
         enddo
      enddo
 !--- if yij < y_cut then combine particles and write in 
      if(yijmin < ycut) then
         phad(nmin1,:)=phad(nmin1,:)+phad(nmin2,:)
         temp(:)=phad(nmin2,:)
         phad(nmin2,:)=phad(maxjet,:)
         phad(maxjet,:)=temp(:)
         maxjet=maxjet-1
         iter=iter-1         
      else
!-----got a jet
         jets=jets+1
      endif

      if(iter < maxjet-1) goto 1
 2    continue
      jets=jets+1
!      write(6,*) jets
      
      do i=1,jets
         pjet(i,:)=phad(i,:)
      enddo
 
      
      return
      end

      function yij_det(p,i,j,Qsq)
      implicit none
      include 'types.f'
      real(dp):: yij_det
      include 'constants.f'
      include 'mxpart.f'
      real(dp) :: p(mxpart,4),Qsq
      integer i,j,nu

      yij_det=(p(i,4)+p(j,4))**2

      do nu=1,3
         yij_det=yij_det-((p(i,nu)+p(j,nu))**2)
      enddo

      yij_det=yij_det/Qsq

      
      return
      end
      
      function yij_dur(p,i,j,Qsq)
      implicit none
      include 'types.f'
      real(dp):: yij_dur
      include 'constants.f'
      include 'mxpart.f'
      real(dp) :: p(mxpart,4),Qsq
      integer i,j,nu
      real(dp) :: pidpj,Emin,cosij

      Emin=min(p(i,4),p(j,4))
      
!      yij_det=(p(i,4)+p(j,4))**2
      pidpj=zip
      do nu=1,3
         pidpj=pidpj+p(i,nu)*p(j,nu)
      enddo
      cosij=pidpj/abs(sqrt(p(i,1)**2+p(i,2)**2+p(i,3)**2))

      cosij=cosij/abs(sqrt(p(j,1)**2+p(j,2)**2+p(j,3)**2))

      
      yij_dur=two*Emin**2*(one-cosij)/Qsq
    
      return
      end
      
      
      subroutine jade_hbb_fullp(p,ycut,jets,pjet,isub,Qsq)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'masses.f'
      include 'npart.f'
      real(dp) :: p(mxpart,4)
      integer jets,i,j,ip
      real(dp) :: ycut,yij
      real(dp) :: phad(mxpart,4),pjet(mxpart,4)
      integer maxjet,iter
      real(dp) :: temp(4),yij_det,yijmin,Qsq,yij_dur   
      integer nmin1,nmin2,isub
      logical is_hadronic
      
      pjet(:,:)=0._dp
      phad(:,:)=0._dp
      jets=0
      maxjet=0
!---- first get hadrons 
!      write(6,*) npart+2-isub
      do i=3,npart+2-isub
         if(is_hadronic(i)) then
            maxjet=maxjet+1
            phad(maxjet,:)=p(i,:)
         endif
      enddo

!       write(6,*) maxjet

!      if(maxjet >= 3) then 
!       write(6,*) maxjet
!       call writeout(phad)
!       pause
!      endif
      iter=0
!-----begin clustering
 1    iter=iter+1

      yijmin=zip
!---- find min y ij
      do i=iter,maxjet
         do j=i+1,maxjet
            yij=yij_det(phad,i,j,Qsq)
            
            if((i==iter).and.(j==i+1)) then
               yijmin=yij
               nmin1=i
               nmin2=j
            elseif(yij < yijmin) then
               yijmin=yij
               nmin1=i
               nmin2=j
            endif
!            write(6,*) i,j,yij,yijmin,nmin1,nmin2
       
!         if(maxjet==3) then
!            write(6,*) i,j,yij
!         endif
         enddo
      enddo

!      write(6,*) 'found ymin =',yijmin
!      write(6,*) 'comp to ',ycut
!      call writeout(phad)
!      pause
 !--- if yij < y_cut then combine particles and write in 
      if(yijmin < ycut) then
 !        if(maxjet ==3) call writeout(phad)
         phad(nmin1,:)=phad(nmin1,:)+phad(nmin2,:)
         temp(:)=phad(nmin2,:)
         phad(nmin2,:)=phad(maxjet,:)
         phad(maxjet,:)=temp(:)
!         if(maxjet ==3)write(6,*) 'clustered,  ',nmin1,nmin2
!         if(maxjet ==3) call writeout(phad)
         maxjet=maxjet-1
         iter=iter-1         
      else
!-----got a jet
         jets=jets+1
      endif

      if(iter < maxjet-1) goto 1
 2    continue
      jets=jets+1
!      write(6,*) jets
      
      do i=1,jets
         pjet(i,:)=phad(i,:)
      enddo
 
      
      return
      end
