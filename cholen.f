c=NAME cholen
C=RCSfile $Source: /u/flosys3d/software/cvsroot/splines/cholen.f,v $
C=VERSION $Revision: 1.1.1.1 $ $Date: 1996/03/05 17:16:35 $
C=TYPE subroutine
C=PURPOSE psplin
C=KEYWORDS  parametric_cubic_splines
C=BLOCK ABSTRACT
C
C get the chord lenght of each arc of a cubic spline
C
C=END ABSTRACT
C
C=BLOCK USAGE
C
C       call CHOLEN(ndimn,n,q,cs)
C
C=END USAGE
C
C=BLOCK SOURCE
C
      subroutine cholen(ndimn,n,q,cs)
c
      real      cs(*),q(ndimn,*)
c
c evaluates chord length cs(i) = |r(i+1)-r(i)|
c
      do 10 i=1,n-1
         cs(i)=0
         do 20 id=1,ndimn
c
            cs(i) = cs(i)+
     1             (q(id,i+1)-q(id,i))**2
c
20       continue
         cs(i) = sqrt(cs(i))
10    continue
      return
      end
C
C=END SOURCE
