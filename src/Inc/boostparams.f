
      real(dp) :: zcut,beta,R0
      logical :: doboostanal
      common/boostparams/zcut,beta,R0
      common/boostparamslog/doboostanal
!$omp threadprivate(/boostparams/) 
