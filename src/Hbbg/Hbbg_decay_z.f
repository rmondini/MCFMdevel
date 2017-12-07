!---- integrated dipoles for H=>b(i1)+bbar(i2)+g(i3)
!---- CW April 17
!---- note that there is no convolution for this piece, only ff dipoles
!---- so we just return a real number to re-weight the LO 
      subroutine Hbbg_decay_z(p,i1,i2,i3,msq_z)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'scale.f'
      include 'qcdcouple.f' 
      include 'nf.f'
      include 'epinv.f'
      include 'epinv2.f'
      real(dp) :: p(mxpart,4),msq_z
      integer i1,i2,i3,i
      real(dp)::ff_qq,ff_gg,ff_gq,fi_gq
      real(dp)::hbbgidips(5),hdecay
      real(dp):: xl12,xl13,xl23
      real(dp):: dipHD_qq,dipHD_qg,dipHD_gg
      real(dp):: dot
      logical nfonly
      common/nfonly/nfonly
!$omp threadprivate(/nfonly/)

      msq_z=zip
      call Hbbg_decay(p,i1,i2,i3,hdecay)
      
      hbbgidips(:)=zip

    
      
      xl12=log(musq/(+two*dot(p,i1,i2)))
      xl13=log(musq/(+two*dot(p,i1,i3)))
      xl23=log(musq/(+two*dot(p,i2,i3)))

      
      hbbgidips(1)=-(ca-2._dp*cf)*dipHD_qq(xl12)*ason2pi
      hbbgidips(2)=+(ca)*dipHD_qq(xl13)*ason2pi*half
      hbbgidips(3)=+(ca)*dipHD_qq(xl23)*ason2pi*half
      hbbgidips(4)=+(ca*dipHD_gg(xl13)*half
     &      +tr*nf*dipHD_qg(xl13))*ason4pi
      hbbgidips(5)=+(ca*dipHD_gg(xl23)*half
     &     +tr*nf*dipHD_qg(xl23))*ason4pi

!----- this bit for only nf 
!       if(nfonly) then
!          hbbgidips(:)=zip
!          hbbgidips(4)=+nf*dipHD_qg(xl13)*ason4pi*tr
!          hbbgidips(5)=+nf*dipHD_qg(xl23)*ason4pi*tr
!       endif

!---- this bit no nf
!          hbbgidips(4)= hbbgidips(4)-nf*dipHD_qg(xl13)*ason4pi*tr
!          hbbgidips(5)= hbbgidips(5)-nf*dipHD_qg(xl23)*ason4pi*tr
      
      
      do i=1,5
         msq_z=msq_z+hbbgidips(i)*hdecay
      enddo

      return
      write(6,*) 'musq = ',musq
      write(6,*) 's12  = ',+two*dot(p,i1,i2)
      write(6,*) 's13  = ',+two*dot(p,i1,i3)
      write(6,*) 's23  = ',+two*dot(p,i2,i3)
      write(6,*) 'epinv = ',epinv,epinv2
      write(6,*) 'ason2pi = ',ason2pi
      write(6,*) 'ason4pi = ',ason4pi
      write(6,*) 'hbbgdips(1) = ',hbbgidips(1)/ason4pi
      write(6,*) 'hbbgdips(2) = ',hbbgidips(2)/ason4pi
      write(6,*) 'hbbgdips(3) = ',hbbgidips(3)/ason4pi
      write(6,*) 'hbbgdips(4) = ',hbbgidips(4)/ason4pi
      write(6,*) 'hbbgdips(5) = ',hbbgidips(5)/ason4pi

      write(6,*) 'int dip = ',msq_z/hdecay/ason2pi
      pause
!      msq_z=msq_z+hdecay*epinv/3._dp*ason2pi*nf
      return 

      end

      function dipHD_gg(L)
      implicit none
      include 'types.f'
      real(dp) :: dipHD_gg
      include 'scheme.f'
      include 'epinv.f'
      include 'constants.f'
      include 'alfacut.f'
  
      real(dp) :: L 
        
      dipHD_gg=100._dp/9._dp + (11._dp*epinv)/3._dp
     &     + 2*epinv**2 + (11._dp*L)/3._dp + 2*epinv*L + L**2 - pisq
        dipHD_gg=dipHD_gg-two*log(aff)**2+two*11._dp/6._dp*(aff-1._dp-log(aff))
  
      if(scheme.eq.'dred') then
         dipHD_gg=dipHD_gg-1._dp/3._dp
      elseif(scheme.ne.'tH-V') then
         write(6,*) 'scheme not recognized in dipHD_gg',scheme
         stop
      endif
      
      return
      end

      



      function dipHD_qg(L)
      implicit none
      include 'types.f'
      real(dp) :: dipHD_qg
      include 'scheme.f'
      include 'epinv.f'
      include 'constants.f'
      include 'alfacut.f'
      real(dp) :: L 
!------ single flavor g=>qq dipole
      dipHD_qg=-(16._dp/9._dp) - (2*epinv)/3._dp - (2*L)/3._dp
!------aff corrections
       dipHD_qg=dipHD_qg-2._dp/3._dp*(aff-1._dp-log(aff))  
      return
      end
      
      function dipHD_qq(L)
      implicit none
      include 'types.f'
      real(dp) :: dipHD_qq
      include 'scheme.f'
      include 'epinv.f'
      include 'constants.f'
      include 'alfacut.f'
      real(dp) :: L 
        
      dipHD_qq=epinv**2 + epinv*(3._dp/2._dp + L) + 1._dp/2._dp*(10._dp
     &     + 3._dp*L + L**2 - pisq)
!      write(6,*) 'L',L
!      write(6,*)  dipHD_qq
!------ aff corrections
       dipHD_qq=dipHD_qq+1.5_dp*(aff-1._dp-log(aff))-log(aff)**2
         
      if(scheme.eq.'dred') then
         dipHD_qq=dipHD_qq-half
      elseif(scheme.ne.'tH-V') then
         write(6,*) 'scheme not recognized in dipHD_qq',scheme
         stop
      endif
      
      return
      end

   
