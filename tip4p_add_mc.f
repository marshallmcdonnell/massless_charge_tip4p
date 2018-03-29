      subroutine tip4p_add_mc (h1_pos, h2_pos, o_pos, om_pos)
      implicit double precision (a-h, o-z)
c     Input/output variables
      double precision, dimension(1:3),intent(in) :: h1_pos,h2_pos,o_pos
      double precision, dimension(1:3),intent(out) :: om_pos
c     Local Variables
      character (len=32) :: cin
      character (len=2)  :: ctemp
      double precision, parameter :: pi = 2.d0*dasin(1.d0)
      double precision, dimension(1:3), parameter :: origin = 0.d0
      double precision, dimension(1:3) :: unit_x, unit_y, unit_z, xyz
      double precision, dimension(1:3) :: rom
      double precision, dimension(1:3) :: norm_vec, unit_norm_vec
      double precision, dimension(1:3) :: h1_o_vec, h2_o_vec
      double precision, dimension(1:3,1:3) :: ucell
      double precision, dimension(1:3,1:3) :: tranm, tranmi
      double precision, dimension(1:3) :: tranmag, tranmagh

c---------------------------------*
c     Start of Code
c---------------------------------*
c     read in input files
      cin = 'in.data'
      open(unit=2,file=cin,form='formatted',status='unknown')
      read (2,*) ctemp 
      read (2,*) ctemp ! input filename
      read (2,*) ctemp 
      read (2,*) ctemp ! output filename
      read (2,*) ctemp 
      read (2,*) ctemp ! massless charge element
      read (2,*) ctemp 
      read (2,*) Vol
      read (2,*) ctemp
      read (2,*) ucell(1,1:3) 
      read (2,*) ucell(2,1:3) 
      read (2,*) ucell(3,1:3) 
      read (2,*) ctemp
      read (2,*) r
      close(unit=2, status='keep')
c     create massless charge vector (spherical coordinates)
      rom(1)  = r         !radius
      rom(2)  = pi/2.d0   !theta (polar angle)
      rom(3)  = 0.d0      !phi (azimuthal angle) ...will be changed in Step 2 (after angle between bonds determined)
c---------------------------------*
c     Get transformation matrix
c---------------------------------*
      adotbxc = ucell(1,1)*ucell(2,2)*ucell(3,3) 
     &        + ucell(1,2)*ucell(2,3)*ucell(3,1)
     &        + ucell(1,3)*ucell(2,1)*ucell(3,2) 
     &        - ucell(1,1)*ucell(2,3)*ucell(3,2)
     &        - ucell(1,2)*ucell(2,1)*ucell(3,3) 
     &        - ucell(1,3)*ucell(2,2)*ucell(3,1)
      vol_1uc = abs(adotbxc) 
      if (abs(vol_1uc) .le. 1.0d-14) then
         print *, ' Your volume of the unit cell matrix is zero.'
         stop
      endif
      call get_tranm(vol_1uc, Vol, ucell, 
     & tranm, tranmi, tranmag, tranmagh)
c---------------------------------------------------------------------------------------------*
c     Step 1: get bond vectors from oxygen to both hydrogens
c---------------------------------------------------------------------------------------------*
      call get_dist(h1_pos, o_pos, tranm, tranmi, h1_o_vec, dis2)
      call get_dist(h2_pos, o_pos, tranm, tranmi, h2_o_vec, dis2)
c---------------------------------------------------------------------------------------------*
c     Step 2: get angle between the bonds and between bond to massless charge
c---------------------------------------------------------------------------------------------*
      call get_angle(h1_o_vec, origin, h2_o_vec, tranm, tranmi,cosphi)
      phi_hoh     = acos(cosphi)      !angle between bonds    (radians)
      phi_hom     = phi_hoh/2.d0      !angle from a bond to massless charge 
      rom(3)      = phi_hom           !set phi angle in massless charge vector
c---------------------------------------------------------------------------------------------*
c     Step 3: create scaled "x-axis" of O-H1 bond w/ massless charge magnitude
c---------------------------------------------------------------------------------------------*
      call get_unit_vec(h1_o_vec, om_pos)
      om_pos = rom(1)*om_pos          !scale unit vector to radius of massless charge from oxygen
c---------------------------------------------------------------------------------------------*
c     Step 4: get normal vector from molecular plane (cross product of O-H bonds)
c---------------------------------------------------------------------------------------------*
      call get_cross_prod (h1_o_vec, h2_o_vec, norm_vec)
c---------------------------------------------------------------------------------------------*
c     Step 5: create unit vector for axis of rotation from normal vector
c---------------------------------------------------------------------------------------------*
      call get_unit_vec(norm_vec, unit_norm_vec)
c---------------------------------------------------------------------------------------------*
c     Step 6: rotate massless charge vector from H1 vector around the normal vector
c             Note 1: still relative to origin
c             Note 2: A COUNTERCLOCKWISE rotation of the coordinate system translates into
c                     a CLOCKWISE rotation of the vector (hence -1 multiplied by angle below)
c                         - See Goldstein, 'Classical Mechanics', 2nd ed., page 164, section 4-7
c---------------------------------------------------------------------------------------------*
      call rotate(om_pos, -1.d0*rom(3), unit_norm_vec, xyz)
      om_pos = xyz
c---------------------------------------------------------------------------------------------*
c     Step 7: set position of massless charge
c---------------------------------------------------------------------------------------------*
      om_pos = o_pos + om_pos
c
      return
      end
