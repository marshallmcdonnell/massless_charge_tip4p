      subroutine get_unit_vec(v, u)
      implicit double precision (a-h, o-z)
c     Variables
      double precision :: mag_vv
      double precision, intent(in),  dimension(1:3) :: v
      double precision, intent(out), dimension(1:3) :: u
c
      call get_mag(v, mag_vv)
      u = v/mag_vv
c
      return
      end