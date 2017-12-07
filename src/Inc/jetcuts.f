      real(dp):: ptjetmin,etajetmin,etajetmax,ptbjetmin,etabjetmax,ptleadjet
      common/jetcuts/ptjetmin,etajetmin,etajetmax,ptbjetmin,etabjetmax,ptleadjet
!$omp threadprivate(/jetcuts/)
