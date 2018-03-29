c
c     This is a simple program that computes the cross product.
c     The output vector (rout) is a normal vector to both input 
c     vectors (u and v).
c
      subroutine get_cross_prod(u, v, rout)
      implicit double precision (a-h,o-z)
c     Variables
      double precision, intent(in), dimension(1:3) :: u, v
      double precision, intent(out), dimension(1:3) ::rout
c-----------------------------------------------*
c     Cross Product of vectors u and v (u x v)
c-----------------------------------------------*
      rout(1) = u(2)*v(3) - u(3)*v(2)
      rout(2) = u(1)*v(3) - u(3)*v(1)
      rout(2) = -1.d0*rout(2)
      rout(3) = u(1)*v(2) - u(2)*v(1)
      return
      end subroutine
      