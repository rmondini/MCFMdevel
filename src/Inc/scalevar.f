      logical doscalevar
      logical, save:: foundpow=.false.
      integer alphaspow,maxscalevar
      real(dp):: scalereweight(6)
      common/doscalevarcommon/doscalevar,maxscalevar
      common/scalevar/alphaspow,scalereweight
!$omp threadprivate(foundpow,/scalevar/)
