c
c     This program will perform a finite rotation of a vector (r)
c     around a given axis (axis) by a given angle (angle) and 
c     return the given vector (r_rotated)
c
c     This rotation is performed using the rotation formula
c     found in "Classical Mechanics" by Goldstein, 2nd Ed.
c     pg. 164-165
c
c     author  Marshall McDonnell
c     Department of Chemical Engineering
c     University of Tennessee, Knoxville
c     last updated  December 5, 2013
c
      subroutine rotate(r, angle, axis, r_rotated)
      implicit double precision (a-h, o-z)
c     Inputs/Outputs
      double precision, intent(in) :: angle   !radians
      double precision, intent(in), dimension(1:3) :: r, axis
      double precision, intent(out), dimension(1:3) :: r_rotated
c     Local Variables
      double precision :: dot_ar
      double precision, dimension(1:3) :: cross_ra
c-------------------------------------------------------------------*
c     Step 1: get dot product and cross product of axis w/ r vector
c
      call get_dot_prod(axis, r, dot_ar)      !axis . r
      call get_cross_prod(r, axis, cross_ra)  !r  x axis
c
c     Step 2: get rotated vector using Rodrigues rotation formula
      r_rotated = r*cos(angle) 
     &          + axis*(dot_ar)*(1 - cos(angle)) 
     &          + cross_ra*sin(angle)
c     
      return
      end
      
      