      subroutine ZZswapcoeff(Xa,Xb,Xrat,Xrats)
      implicit none
      include 'types.f'
      integer h3,h5
      complex(dp):: Xa(2,2),Xb(2,2),tmp(2,2),
     & Xrat(2,2,2,2),Xrats(2,2,2,2)
      
      tmp(:,:)=Xa(:,:)
      do h3=1,2
      do h5=1,2
      Xa(h3,h5)=tmp(3-h3,3-h5)
      enddo
      enddo

      tmp(:,:)=Xb(:,:)
      do h3=1,2
      do h5=1,2
      Xb(h3,h5)=tmp(3-h3,3-h5)
      enddo
      enddo
      
      do h3=1,2
      do h5=1,2
      Xrat(2,1,h3,h5)=Xrats(1,2,3-h3,3-h5)
      Xrat(1,1,h3,h5)=Xrats(2,2,3-h3,3-h5)
      enddo
      enddo
      
      return
      end
      
