c=NAME slen
C=RCSfile $Source: /u/flosys3d/software/cvsroot/splines/slen.f,v $
C=VERSION $Revision: 1.1.1.1 $ $Date: 1996/03/05 17:16:40 $
C=TYPE subroutine
C=PURPOSE general
C=KEYWORDS  parametric_cubic_splines
C=BLOCK ABSTRACT
C
c It evaluate the lenght of the arcs of a cubic spline
c curve by using a Romberg type integration.
c
C=END ABSTRACT
C
C=BLOCK USAGE
c
c    call SLEN(ndimn,n,q,t,cs,len,eps)
c
c INPUT       DIMENSION     DESCRIPTION
c ndimn                     N. of dimensions 
c                           (IT MUST BE equal to 2 or 3 !!)
c n                         n. of knots
c q           ndimn,n       Knots coordinates
c cs          n             Arc chord length
C t           ndimn,n       Tangent at knots
c eps                       Fractional accuracy for numerical integration.
c
c  OUTPUT     DIMENSION     DESCRIPTION
c
c  len        n-1           arcs length
c
c
c SUBROUTINE CALLED: coeff -> compute coefficents
c                    lengt -> compute integral
c
c FORMULA USED:
c
c  len(i) = int[u(i),u(i+1)](||x,u|| du)  (where int[a,b] indicates an integral from
c                                          a to b)
c
C=END USAGE
C
C=BLOCK SOURCE
C
      subroutine slen(ndimn,n,q,t,cs,len,eps)
      real len(*),scoef(5)   ,cs(*)
c     real q(ndimn,*),t(ndimn,*),scoef(5)
      in=1
      u1=0
      u2=1
      do 10 i=1,n-1
c
c evaluate length**2 polinomial coeficents
c
         call coeff(ndimn,n,q,t,cs,i,scoef)
c
         eps1=eps*cs(i)
         a1 = scoef(1)
         a2 = scoef(2)
         a3 = scoef(3)
         a4 = scoef(4)
         a5 = scoef(5)
         call    lengt(in,a1,a2,a3,a4,a5,u1,u2,s,eps)
         len(i)=s
10    continue
      return
      end
C
C=END SOURCE
