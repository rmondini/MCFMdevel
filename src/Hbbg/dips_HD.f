!---- modification of regular dipolesub.f to allow for decay routines
!---- routine is for H=>4 parton decays, could be generalized later if needed
!---  CW
!----- Many (only thus far) modification is to change msq into constants from arrays
      subroutine dips_HD(nd,p,ip,jp,kp,sub,subv,msq,msqv,
     & subr_born,subr_corr)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'qcdcouple.f'
      include 'qqgg.f'
      include 'ptilde.f'
      include 'alfacut.f'
      include 'kprocess.f'
      include 'dynamicscale.f'
      include 'initialscales.f'
      include 'dipolescale.f'
      include 'facscale.f'
      include 'incldip.f'
      real(dp):: p(mxpart,4),ptrans(mxpart,4),sub(4),subv,vecsq
      real(dp):: x,omx,z,omz,y,omy,u,omu,sij,sik,sjk,dot,vec(4)
      real(dp):: msq,msqv,vtilde
      integer:: nd,ip,jp,kp,nu,j,k,ipt,i
c--      logical:: includedipole
      external subr_born,subr_corr
      real(dp) :: ptemp(mxpart,4)
      integer isp_dipmin,isp_dipnpart
      common/isp_dip/isp_dipmin,isp_dipnpart
!$omp threadprivate(/isp_dip/)
      
C---  Initialize the dipoles to zero
      do j=1,4
      sub(j)=0._dp
      enddo
      subv=0._dp
      msq=0._dp
      msqv=0._dp
      incldip(nd)=.true.
      
      sij=two*dot(p,ip,jp)
      sik=two*dot(p,ip,kp)
      sjk=two*dot(p,jp,kp)
!-----DIPOLE IS ALWAYS FINAL-FINAL
      
      y=sij/(sij+sjk+sik)
     
C---  Modification so that only close to singular subtracted
      if (y > aff) then
         incldip(nd)=.false.
         return
      endif
        
      z=sik/(sjk+sik)
      omz=one-z
      omy=one-y
C---  calculate the ptrans-momenta 
      call transform(p,ptrans,y,ip,jp,kp)
      call storeptilde(nd,ptrans)
      ptemp(:,:)=zip
!------need to position the gluon correctly for ME
!-----isp does nothing in dipole, but is used to position gluon
      do i=1,isp_dipnpart
         ptemp(i,:)=ptrans(isp_dipmin+i-1,:)
      enddo
     
!      call writeout(p) 
!      call writeout(ptemp)
      
c--   Check to see if this dipole will be included
c     incldip(nd)=includedipole(nd,ptrans)
c     if (incldip(nd) .eqv. .false.) return
      
      do nu=1,4
         vec(nu)=z*p(ip,nu)-omz*p(jp,nu)
      enddo

!      write(6,*) vec(4)*ptemp(3,4)-vec(3)*ptemp(3,3)
!     &     -vec(2)*ptemp(3,2)-vec(1)*ptemp(3,1)
c---  if using a dynamic scale, set that scale with dipole kinematics
      if (dynamicscale) then
         call scaleset(initscale,initfacscale,ptrans)
         dipscale(nd)=facscale
      endif
      
      call subr_born(ptemp,msq)
      ipt=3
      
      call subr_corr(ptemp,vec,ipt,msqv)
      
      sub(qq)=gsq/sij*(two/(one-z*omy)-one-z)
      sub(gq)=gsq/sij
      sub(gg)=gsq/sij*(two/(one-z*omy)+two/(one-omz*omy)-four)
      subv   =+4._dp*gsq/sij/sij
      
        
      return
      end
      
