      real(dp) pincj(mxpart,4)
      logical scet_inc_jets
      integer scetjets
      common/scetincj/pincj
      common/scetjets/scetjets
!$omp threadprivate(/scetincj/,/scetjets/)
      common/scetinc_log/scet_inc_jets
