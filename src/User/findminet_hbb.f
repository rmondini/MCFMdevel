      subroutine findminet_hbb(p,pjet,pjetmin,pjetmax,dkmin,nk,ipow)
      implicit none
      include 'types.f'
c--- this finds the minimum dkmin for pjet indices pjetmin through pjetmax
c--- returns dijmin and indices of minimum in (nmin1,nmin2)
C--- calculate the beam proto-jet separation see NPB406(1993)187, Eqn. 7
C---  S.~Catani, Y.~L.~Dokshitzer, M.~H.~Seymour and B.~R.~Webber
C--- in practice this is just the minimum ptsq of protojets       
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      real(dp):: p(mxpart,4),pjet(mxpart,4),dkmin,dk,pt
      integer:: pjetmin,pjetmax,nk,i,ipow

      dkmin=pt(pjetmin,pjet)
      if (ipow .ne. 1) dkmin=dkmin**(ipow)
      nk=pjetmin

c--- if only one entry, this must be the minimum
      if (pjetmin+1 > pjetmax) return
            
      do i=pjetmin+1,pjetmax
!        dk=pt(i,pjet)
        dk=pjet(i,4)
        if (ipow .ne. 1) dk=dk**(ipow)
        if (dk < dkmin) then
          dkmin=dk
          nk=i
        endif
      enddo
            
      return
      end
      
       subroutine findmind_hbb(p,pjet,pjetmin,pjetmax,dijmin,nmin1,nmin2,
     &                    ipow)
      implicit none
      include 'types.f'
c--- this finds the minimum dij for pjet indices pjetmin through pjetmax
c--- returns dijmin and indices of minimum in (nmin1,nmin2)
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      real(dp):: p(mxpart,4),pjet(mxpart,4),dijmin,dij_hbb,d
      integer:: pjetmin,pjetmax,nmin1,nmin2,i,j,ipow

      do i=pjetmin,pjetmax
        do j=i+1,pjetmax
          d=dij_hbb(p,pjet,i,j,ipow)
          if ((i == pjetmin) .and. (j == i+1)) then
            dijmin=d
            nmin1=i
            nmin2=j
          elseif (d < dijmin) then
            dijmin=d
            nmin1=i
            nmin2=j
          endif
        enddo
      enddo
      
      return
      end


      function dij_hbb(p,pjet,i,j,ipow)
      implicit none
      include 'types.f'
      real(dp):: dij_hbb
C---calculate the proto-jet separation see NPB406(1993)187, Eqn. 7
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      integer:: i,j,ipow
      real(dp):: p(mxpart,4),pjet(mxpart,4),pti,ptj,pt,r
c      real(dp):: etarap,yi,yj,phii,phij
      
!      pti=pt(i,pjet)
!      ptj=pt(j,pjet)
      pti=pjet(i,4)
      ptj=pjet(j,4)
      
c--- old method - bad because (phii-phij) can be > pi       
c      yi=etarap(i,pjet)
c      yj=etarap(j,pjet)

c      phii=atan2(pjet(i,1),pjet(i,2))
c      phij=atan2(pjet(j,1),pjet(j,2))
      
c      dij=sqrt((yi-yj)**2+(phii-phij)**2)

c--- new method - r() calculates true value of 0 < (phi-phij) < pi
      dij_hbb=r(pjet,i,j)
              
      if (ipow .ne. 1) then
        pti=pti**(ipow)
        ptj=ptj**(ipow)
      endif
      dij_hbb=dij_hbb*min(pti,ptj)
      
      return
      end
      
     
