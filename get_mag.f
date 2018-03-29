      subroutine get_mag(v, mag_vv)
      implicit double precision (a-h, o-z)
c     Variables
      double precision :: dot_vv
      double precision, intent(in), dimension(1:3) :: v
      double precision, intent(out) :: mag_vv
c
      mag_vv = 0.d0
      call get_dot_prod(v, v, dot_vv)
      mag_vv = dsqrt(dot_vv)
c
      return
      end