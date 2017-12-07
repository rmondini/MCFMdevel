*****************
* LHAPDF version*
*****************
      subroutine fdist_select(ih,ibeam,x,xmu,fx)
! This should be identical to the usual fdist routine, but it masks
! pdfs according to selectpdfs(ibeam), passed via common block
      implicit none
      include 'selectpdfs.f'
      double precision fx(-5:5),x,xmu,fPDF(-6:6),fphoton
      integer Iprtn,ih,ibeam
      logical has_photon
c---  ih1=+1 proton 
c---  ih1=-1 pbar

C---set to zero if x out of range
      if (x .ge. 1d0) then
          do Iprtn=-5,5
             fx(Iprtn)=0d0
          enddo
          return
      endif
!$omp critical(LHApdf) 
      if (has_photon()) then
        call evolvePDFphoton(x,xmu,fPDF,fphoton)
      else
        call evolvePDF(x,xmu,fPDF)
      endif
!$omp end critical(LHApdf) 
      if (ih.eq.1) then
        do Iprtn=-5,5
          fx(+Iprtn)=fPDF(+Iprtn)/x
        enddo
      elseif(ih.eq.-1) then
        do Iprtn=-5,5
          fx(+Iprtn)=fPDF(-Iprtn)/x
        enddo
      endif

c--- explicitly set some pdfs to zero according to selectpdfs
      if     (selectpdfs(ibeam) == 0) then
        fx(1:5)=0d0    ! remove quark pdfs
        fx(-5:-1)=0d0  ! remove antiquark pdfs
      elseif (selectpdfs(ibeam) == +1) then
        fx(-5:0)=0d0  ! remove antiquark and gluon pdfs
      elseif (selectpdfs(ibeam) == -1) then
        fx(0:5)=0d0  ! remove quark and gluon pdfs
      endif
      
      return
      end

  

