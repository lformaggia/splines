C=NAME gtimp
C=RCSfile $Source: /u/flosys3d/software/cvsroot/splines/gtimp.f,v $
C=VERSION $Revision: 1.1.1.1 $ $Date: 1996/03/05 17:16:38 $
C=TYPE subroutine
C=PURPOSE local_to psplin
C=KEYWORDS parametric_cubic_splines
C=BLOCK ABSTRACT
C
C It inserts a value into the imposed tangent vector timp
C
C=END ABSTRACT
C=BLOCK USAGE
C
C call gtimp(ndimn,timp,q,n,is)
C
C  INPUT       DIMENSION      DESCRIPTION
c
c ndimn                     N. of dimensions
c n                         n. of knots
c q           ndimn,n       Knots coordinates
c is                        Switch:
c                                         =1 1st end of curve  =2 2nd end of curve
c
C OUTPUT      DIMENSION     DESCRIPTION
c
c timp        ndimn,2       imposed tangent value
c
c
c This routine inserts into timp a value for the tangent at the end of the
c curve, to be later used for the cubic spline interpolation. The value inserted is
c the direction versor defined by the first (last) two knots on the curve, i.e.
c
c is =1 -> timp(*,1) = (q(*,2)-q(*,1)  )/||q(*,2)-q(*,1)  ||
c is =2 -> timp(*,2) = (q(*,n)-q(*,n-1))/||q(*,n)-q(*,n-1)||
c
c This imposed value is the one used if ispty=1 in psplin or psplin2
c
c
C=END USAGE
C
C=BLOCK SOURCE
C
      subroutine gtimp(ndimn,timp,q,n,is)
      parameter(zero=1.e-10)
c
      real timp(ndimn,2),q(ndimn,*)
c
      if(is.eq.1)i=1
      if(is.eq.2)i=n-1
      x1 = 0.
      do 10 id=1,ndimn
         timp(id,is)=q(id,i+1) - q(id,i)
         x1         =x1 + timp(id,is)**2
10    continue
      if(x1.lt.zero)then
         x1=1.
      else
         x1 = 1./sqrt(x1)
      endif
      do 20  id=1,ndimn
         timp(id,is)=timp(id,is)*x1
20    continue
      return
      end
C
C=END SOURCE
