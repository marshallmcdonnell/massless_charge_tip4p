c     
c     This is a simple program to compute the angle between two vectors
c      
c      date begun:  February 27, 2013
c     last updated:  February 27, 2013
c
c 
      subroutine get_angle(r1, r_vertex, r2, tranm, tranmi, costheta)
      implicit double precision (a-h,o-z)
c
c     passed variables
      double precision, intent(in), dimension(1:3) :: r1, r_vertex, r2
      double precision, intent(in), dimension(1:3,1:3) :: tranm, tranmi
c
c     returned variables
      double precision, intent(out) :: costheta
c     
c     local variables
      double precision, dimension(1:3) :: dis1, dis2, dist1, dist2
      double precision :: dotproduct, sqrdis1, sqrdis2, dis_1, dis_2, 
     & denom
c
c     compute distances
c
      
      dis1(1:3) = r_vertex(1:3) - r1(1:3)
      dis2(1:3) = r_vertex(1:3) - r2(1:3)
                  do k = 1, 3, 1
                     dist1(k) = sum(tranm(k,1:3)*dis1(1:3))
                     dist2(k) = sum(tranm(k,1:3)*dis2(1:3))
                  enddo
c            
c     apply minimum image convention
c					
                  do k = 1, 3, 1
                        dist1(k) = dist1(k) - nint(dist1(k))
                        dist2(k) = dist2(k) - nint(dist2(k))
                  enddo
c
c     transform from unit cell coords to simulation coords
c
                  do k = 1, 3, 1
                     dis1(k) = sum(tranmi(k,1:3)*dist1(1:3))
                     dis2(k) = sum(tranmi(k,1:3)*dist2(1:3))
                  enddo
c
c     compute the dot product of vector 1(r1 to r_vertex) and vector 2 (r2 to r_vertex)
c
      dotprod = sum( dis1(1:3)*dis2(1:3) )
c
c     compute the 2-norm of both vector 1 an vector 2
c
c
c     compute the squared distance
c     
      sqrdis1 = sum(dis1*dis1)
      sqrdis2 = sum(dis2*dis2)
c
c     compute the distance
c
      dis_1 = sqrt(sqrdis1)
      dis_2 = sqrt(sqrdis2)
c
c     multiply the 2-norms = denominator
c
      denom = dis_1*dis_2
c
c     compute and return cos(theta)
c
      costheta = dotprod/denom
                  
                  
      return
      end subroutine
      
