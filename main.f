      program add_massless_charge
      implicit double precision (a-h, o-z)
c     Variables
      character (len=2)  :: mc_element
      character (len=32) :: cin, ctemp
      character (len=64) :: input_xyz, output_xyz
      character (len=2), allocatable :: cname(:)
      double precision, allocatable :: r(:,:)
      double precision, dimension(1:3) :: h1_pos, h2_pos, o_pos, om_pos

c-----------------------------*
c     open input control file
c     to get filenames
c-----------------------------*
      cin = 'in.data'
      open(unit=2,file=cin,form='formatted',status='unknown')
      read (2,*) ctemp
      read (2,*) input_xyz
      read (2,*) ctemp
      read (2,*) output_xyz
      read (2,*) ctemp
      read (2,*) mc_element
      close(unit=2, status='keep')
 
c-----------------------------*
c     open input coord. file
c     and read in positions
c-----------------------------*
      open(unit=1, file=input_xyz,form='formatted',status='unknown')
      read (1,*) Natoms
      allocate (cname(1:Natoms), r(1:Natoms,1:3) )
      read (1,*) 
      do i = 1, Natoms, 1
          read (1,*) cname(i), r(i,1:3)
      enddo
      close(unit=1,status='keep')
      Nmolec = Natoms/3
c-----------------------------*
c     loop back over and 
c     add massless charge 
c     and write to new output
c     xyz file
c-----------------------------*
      open(unit=3, file=output_xyz,form='formatted',status='unknown')
      write(3,*) Natoms + Nmolec
      write(3,*) ''
      do i = 1, Nmolec, 1
          i_count_h = 0
          do j = 1, 3, 1
              if( cname((i-1)*3 + j) .eq. 'O') then
                  o_pos(1:3) = r((i-1)*3 + j, 1:3)
              elseif( cname((i-1)*3 + j) .eq. 'H') then
                  if(i_count_h .eq. 0) then
                      h1_pos(1:3) = r((i-1)*3 + j, 1:3)
                      i_count_h = i_count_h + 1
                  elseif(i_count_h .eq. 1) then
                      h2_pos(1:3) = r((i-1)*3 + j, 1:3)
                      i_count_h = i_count_h + 1
                  endif
              else
                  print *, 'molecule', i, 'atom', j, 
     &             'w/ name', cname((i-1)*3 + j), 'not O or H'
              endif
              write(3,1007) cname((i-1)*3 + j), r((i-1)*3 + j, 1:3)
          enddo
          call tip4p_add_mc(h1_pos, h2_pos, o_pos, om_pos)
          write(3,1007) mc_element, om_pos(1:3)
      enddo
      close(unit=3,status='keep')
c      read *, paused
 1007 format(a2,1x,3(f12.8,1x))
      stop
      end
      
