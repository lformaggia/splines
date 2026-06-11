C=NAME loas1d
C=RCSfile $Source: /u/flosys3d/software/cvsroot/splines/loas1d.f,v $
C=VERSION $Revision: 1.1.1.1 $ $Date: 1996/03/05 17:16:38 $
C=TYPE subroutine
C=PURPOSE genera;
C=KEYWORDS partametric_polynimial_curves Data_base_interface
C=BLOCK ABSTRACT
c
C It reads the geometrical information for a parametric  curve
c from data file structure and stores it into isn and xln1d
C
C
C=END ABSTRACT
C=BLOCK USAGE
C
c     call loas1d(isn,xln1d,is,ksio,isio,n,isg,isl1,isl2,
c     1                  ug1,ug2,xg,marc,msel,ito,ierr)
c
c  INPUT      DIMENSION     DESCRIPTION
c
c is                        curve local numbering
c ksio                      i/o unit number associated to the file of keys
c isio                      i/o unit number associated to the file of coefficents
c marc                      max. dimesnion of isn
c msel                      max dimension of xln1d
c ito                       i/o unit number associated to error i/o channel
c
c
c  OUTPUT     DIMENSION     DESCRIPTION
c
c n                         n. of knots (n. of arcs+1)
c isg                       global numbering to be attributed to the curve
c xln1d       3,msel        array of coefficents
c isn         marc          keys to array of  coefficents
c ug1,ug2                   global coordinate of end points of definition curve
c isl1,isl2                 local  numbering of the surfaces adjacent to the curve
c xg                        control parameter (optional)
c ierr                      error indicator (see DIAGNOSTICS)
c
c * DIAGNOSTICS*
c
c  ierr =0                  no error
c  ierr =1                  n. of arcs in curve exceeds marc
c  ierr =2                  n. of coefficents for curve exceeds msel
c  ierr =3                  ierr=1 + ierr=2
c
c  * NOTE *
c
c  if ierr/=0 a message is sent to the unit associated to ito. In standard FORTRAN
c  if you want the message on the standard output you should put ito=6
C
C=END USAGE
C
C=BLOCK SOURCE
C
      subroutine loas1d(isn,xln1d,is,ksio,isio,narc,isg,isl1,isl2,
     1                  ug1,ug2,xg,marc,msel,ito,ierr)
      implicit none
C
C IT LOADS CURVE COEFFICENTS FROM FILE
C
C
      integer*8 i,ierr,is,isg,isio,isl1,isl2,isn(marc),ito,j,key
      integer*8 ksio,marc,msel,narc,ndeg,ips
      real*8 ax,ay,az,ug1,ug2,xg,xl,xln1d(3,msel)
C
      ips = 0
      ierr =0
      read(ksio,rec=is)key,narc,isg,ug1,ug2,isl1,isl2,xg
      key = key-1
      do 10 i=1,narc-1
        ips=ips+1
        if(i.lt.marc)isn(i)=ips
        read(isio,rec=key+ips)ndeg,xl
        xln1d(1,ips)=dble(ndeg)
        xln1d(2,ips)=xl
        do 20 j=1,ndeg
          ips = ips+1
          read(isio,rec=key+ips)ax,ay,az
          if(ips.gt.msel)go to 20
          xln1d(1,ips)  =ax
          xln1d(2,ips)  =ay
          xln1d(3,ips)  =az
20      continue
10    continue
      if(narc-1.gt.marc)then
        write(ito,*)' ERROR IN LOAS1D : 001'
        write(ito,*)' TOO MANY ARCHS ON CURVE n.',isg
        write(ito,*)' (local numbering = ',is,' )'
        write(ito,*)' increase MARC to at least ',narc-1
        write(ito,*)' File retrieving aborted'
        ierr =1
      endif
      if(ips.gt.msel)then
        write(ito,*)' ERROR IN LOAS1D : 002'
        write(ito,*)' TOO MANY COEFFICENTS ON CURVE n.',isg
        write(ito,*)' (local numbering = ',is,' )'
        write(ito,*)' increase MSEL to at least ',IPS
        write(ito,*)' File retrieving aborted'
        ierr = ierr + 2
      endif
      return
      end
C
C=END SOURCE
