C=NAME buicur
C=RCSfile $Source: /u/flosys3d/software/cvsroot/splines/buicur.f,v $
C=VERSION $Revision: 1.1.1.1 $ $Date: 1996/03/05 17:16:34 $
C=TYPE subroutine
C=AUTHOR L.formaggia
C=PURPOSE General
C=KEYWORDS Parametric_cubic_splines Parametric_polynomial_curves Data_structure_interface
C
C=BLOCK ABSTRACT
C
c It converts the rapresentation of a cubic spline curve into the general
c rapresentation for a polynomial curve
C
C=END ABSTRACT
C=BLOCK USAGE
C
C    call buicur(n,q,t,cs,isn,xln1d,marc,msel,ierr)
C
C  INPUT      DIMENSION     DESCRIPTION
C
C
C  n                        n. of knots
c  t          3,n           derivative vector at knots (see routine psplin)
C  cs         2,n-1         arc chord length  (see routine psplin)
C  q          3,n           knots coordinates
c  marc                     max. number of arc
c  msel                     max. dimension of array of coefficents
C
C OUTPUT      DIMENSION     DESCRIPTION
c
C  isn        marc          array of pointers to xln1d
C  xln1d      msel          array of coefficents
c  ierr                     error condition
C
C
C   * DIAGNOSTICS *
c
c ierr = 0    no error
c ierr = 1    number of arcs exceeded (marc<(n-1))
c ierr = 2     number of coefficents exceeded (msel<5*(n-1))
c ierr = 3    1+2
C
C= END USAGE
C
C=BLOCK SOURCE
C
      subroutine buicur(n,q,t,cs,isn,xln1d,marc,msel,ierr)
c
c builds internal rapresentation for a curve from the arrays q t and cs
c which describe the curve in terms of a cubic spline
c
      integer isn(*)
      real xln1d(3,*),t(3,n),cs(n),q(3,*)
c
      ierr =0
      inkseg  =1
c
      if(n.ge.marc)then
         ierr = ierr +1
      endif
      if(5*(n-1) .gt.msel)then
         ierr = ierr +2
      endif
      if(ierr.ne.0)return
c
      do 20 i=1,n-1
        isn(i)=inkseg
        xln1d(1,inkseg)=4
        xln1d(2,inkseg)=cs(i)
        inkseg=inkseg+1
        call patc1d(3,xln1d(1,inkseg),q,t,cs,i)
        inkseg = inkseg+4
20    continue
      return
      end
C
C=END SOURCE
