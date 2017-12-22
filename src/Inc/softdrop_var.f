      integer :: firstjet,npjr(mxpart,mxpart),pinSDj(mxpart,mxpart)
      common/softdrop_variables/firstjet,npjr,pinSDj
      integer nj,noj,idlptjet,nplj
      common/njets_boost/nj,noj,idlptjet,nplj
!$omp threadprivate(/softdrop_variables/)
!$omp threadprivate(/njets_boost/)
