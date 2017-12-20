      integer :: firstjet,npjr(mxpart,mxpart)
      common/softdrop_variables/firstjet,npjr
      integer nj,noj,idlptjet,nplj
      common/njets_boost/nj,noj,idlptjet,nplj
!$omp threadprivate(/softdrop_variables/)
!$omp threadprivate(/njets_boost/)
