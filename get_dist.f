c     
c     This is a simple program to compute the distance between two vectors      
c      
c      date begun:  February 1, 2013
c     last updated:  February 22, 2013
c
c 
      subroutine get_dist(r1, r2, tranm, tranmi, dis, dis2)
      implicit double precision (a-h,o-z)
c
c     passed variables
      double precision, intent(in), dimension(1:3) :: r1, r2
      double precision, intent(in), dimension(1:3,1:3) :: tranm, tranmi
c
c     returned variables
      double precision, intent(out), dimension(1:3) :: dis
      double precision, intent(out) :: dis2
c     
c     local variables
      double precision, dimension(1:3) :: dist
c
c     compute distances
c
      
      dis(1:3) = r1(1:3) - r2(1:3)
c
c     transform from simulation coords to unit cell coords
c
                  do k = 1, 3, 1
                     dist(k) = sum(tranm(k,1:3)*dis(1:3))
                  enddo
c            
c     apply minimum image convention
c					
                  do k = 1, 3, 1
                     dist(k) = dist(k) - nint(dist(k))       !slow version of code
c                      if (dist(k) .gt. 0.5d0) dist(k) = dist(k) - 1.d0
c                      if (dist(k) .lt. 0.5d0) dist(k) = dist(k) + 1.d0
                  enddo
c
c     transform from unit cell coords to simulation coords
c
                  do k = 1, 3, 1
                     dis(k) = sum(tranmi(k,1:3)*dist(1:3))
                  enddo
c
c     return the squared distance
c     
                  dis2 = sum(dis*dis)
                  
      return
      end subroutine
      
