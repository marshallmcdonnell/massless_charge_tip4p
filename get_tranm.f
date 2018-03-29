      subroutine get_tranm(vol_1uc, Vol, ucell, 
     & tranm, tranmi, tranmag, tranmagh)
      implicit double precision (a-h, o-z)
c     passed variables
      double precision, intent(in) :: vol_1uc, Vol
      double precision, intent(in), dimension(1:3,1:3) :: ucell
      double precision, intent(out), dimension(1:3,1:3) :: tranm
      double precision, intent(out), dimension(1:3,1:3) :: tranmi
      double precision, intent(out), dimension(1:3) :: tranmag
      double precision, intent(out), dimension(1:3) :: tranmagh
c
c     determine number of unit cells 
cp      print *, ' get_tranm: Vol ', Vol
cp      print *, ' get_tranm: vol_1uc ', vol_1uc
      xncell = Vol/vol_1uc 
c     because we are in a 3-d system
      xnfact = xncell**(1.d0/3.d0)
c     calculate the transformation matrix
      tranmi(1:3,1:3) = xnfact*ucell(1:3,1:3) 
c     calculate the inverse of the transformation matrix
cp      print *, ' get_tranm: xnfact ', xnfact 
cp      print *, ' get_tranm: ucell ', ucell 
cp      print *, ' get_tranm: tranmi ', tranmi
      call xinvmat (tranmi, tranm, deta, 3)
cp      print *, ' deta ', deta 
      if (abs(deta) .le. 1.0d-14) then 
         print *, ' Your determinant of the unit cell matrix is zero.'
         stop
      endif
c     calculate the magnitude of the vectors of the transformation matrix
      do i = 1, 3, 1
         tranmag(i) = sqrt(sum( tranmi(1:3,i)*tranmi(1:3,i) ) )
      enddo
      tranmagh = 0.5d0*tranmag
c
      return
      end