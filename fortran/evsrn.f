C=NAME evsrn
C=RCSfile $Source: /u/flosys3d/software/cvsroot/splines/evsrn.f,v $
C=VERSION $Revision: 1.1.1.1 $ $Date: 1996/03/05 17:16:36 $
C=TYPE subroutine
C=PURPOSE general
C=KEYWORDS parametric_surface
C=BLOCK ABSTRACT
C
C It evaluates surface normal, given the derivatives
C
C
C
C=END ABSTRACT
C=BLOCK USAGE
C
C call evsrn(ru,rv,xnor,anor)
C
C
C  INPUT       DIMENSION      DESCRIPTION
c
c  ru         3             r,u
c  rv         3             r,v
c
C OUTPUT      DIMENSION     DESCRIPTION
c
C xnor        3             surface normal: r,u * r,v /|| r,u *r,v||
C anor                      modulus of derivative cross product:
c                           ||r,u * r,v||
C
C
C=END USAGE
C
C=BLOCK SOURCE
C
      subroutine evsrn(ru,rv,xnor,anor)
c
c
c evaluate normal to curve
c
      implicit none
      real*8 ano,anor,ru(3),rv(3),xnor(3)
c
      xnor(1) = ru(2)*rv(3) - rv(2)*ru(3)
      xnor(2) = ru(3)*rv(1) - rv(3)*ru(1)
      xnor(3) = ru(1)*rv(2) - rv(1)*ru(2)
      anor = sqrt(xnor(1)*xnor(1)+xnor(2)*xnor(2)+xnor(3)*xnor(3))
      if(anor.eq.0.d0)then
        xnor(1)=0.d0
        xnor(2)=0.d0
        xnor(3)=0.d0
        return
      endif
c
      ano = 1.d0/anor
      xnor(1) = xnor(1)*ano
      xnor(2) = xnor(2)*ano
      xnor(3) = xnor(3)*ano
      return
      end
C
C=END SOURCE
