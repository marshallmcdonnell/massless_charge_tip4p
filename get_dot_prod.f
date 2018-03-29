      subroutine get_dot_prod(u, v, dot_uv)
      implicit double precision (a-h, o-z)
c     Variables
      double precision, intent(in), dimension(1:3) :: u, v
      double precision, intent(out) :: dot_uv
c
      dot_uv = 0.d0
      do i = 1, 3, 1
          dot_uv = dot_uv + u(i)*v(i)
      enddo
      return
      end