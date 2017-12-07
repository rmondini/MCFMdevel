!---- rotuine for determining an array of groomed-jets
!---- using soft-drop, as a by product the anti-kt clustered
!---- events are also returned here so the regular clustering need not occur

!      subroutine gen_softdrop(pparton,pkt,psoftd,ktjets,sdjets,isub)
      subroutine gen_softdrop(isub)
      implicit none
      include 'types.f' 
      include 'constants.f'
      include 'mxpart.f'
      include 'npart.f'
      include 'jetcuts.f'
      include 'jetlabel.f'
      include 'plabel.f'
      include 'is_functions_com.f'
      real(dp) :: pparton(mxpart,4)
      real(dp) :: pkt(mxpart,4),psoftd(mxpart,4)
      integer ktjets,sdjets,startj(mxpart)
      integer i,isub,maxjet
     
      sdjets=0
      ktjets=0
      
      startj(:)=0
      psoftd(:,:)=zip
      pkt(:,:)=zip


      maxjet=0
!---- first copy over  EW particles
      do i=3,npart+2-isub
!         psoft(i,:)=pparton(i,:)
!         pkt(i,:)=pparton(i,:)
      enddo

      write(6,*) jetindex
      write(6,*) jets
      pause
      
!-----cluster with standard jet algorithim, but keep tabs on which partons entered jet

!----- strip to subjets using CA to make clustering trees

!--- check soft drop condition 
      
      return
      end 

      
