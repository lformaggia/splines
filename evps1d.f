C=NAME evps1d
C=RCSfile $Source: /u/flosys3d/software/cvsroot/splines/evps1d.f,v $
C=VERSION $Revision: 1.1.1.1 $ $Date: 1996/03/05 17:16:36 $
C=TYPE subroutine
C=PURPOSE general
C=KEYWORDS parametric_Polynomial_curves
C=BLOCK ABSTRACT
C
C It computes  position vector and derivatives on a parametric curve in
c the tridimensional space
c
C=END ABSTRACT
C=BLOCK USAGE
C
C call evps1d(xl,u,xp,xd,xdd,xddd,sz,in)
C
C  INPUT       DIMENSION      DESCRIPTION
c
c xl          3,*           array containing the coefficents of the arc
c                           in which the point lies.
c u                         parametric coordinate local to the arc
c in                        switch which determines the output:
c                           0  -> position vector, 1st,2nd and 3rd derivative
c                           1  -> only position vector
c                           2  -> position vector and 1st derivative
c                           -1 -> only 1st derivative
c
C OUTPUT      DIMENSION     DESCRIPTION
c
c xp          3             position vector          r
c xd          3             first derivative         r,u
c xdd         3             2nd derivative           r,uu
c xddd        3             third derivative         r,uuu
c sz                        first derivative modulus ||r,u||
c
c NOTE : the degree of the polinomial is contained in the first location
c        of array xl, following the convention described in the User's Reference Manual.
c
C=END USAGE
C
C=BLOCK SOURCE
C
      subroutine evps1d(xl,u,xp,xd,xdd,xddd,sz,in)
c
c  in =0 all in=1 only point in=2 point+ 1st derivative
c            in=-1 only 1st derivative
c
      implicit none
      integer*8 id,in,k,ndeg
      real*8 sz,u,xddd(3),xdd(3),xd(3),xl(3,0:*),xm,xp(3)
c
      sz = 0.d0
      ndeg = int(xl(1,0))
      if(in.eq.-1)go to 31
c
      do 1 id=1,3
1     xp(id)=xl(id,ndeg)
      do 2 k=ndeg-1,1,-1
      xp(1) = xp(1)*u+xl(1,k)
      xp(2) = xp(2)*u+xl(2,k)
      xp(3) = xp(3)*u+xl(3,k)
2     continue
31    continue
c
      if(in.eq.1)return
c
      do 11 id=1,3
11    xd(id)=(ndeg-1)*xl(id,ndeg)
      do 12 k=ndeg-1,2,-1
      xm = dble(k-1)
      xd(1) = xd(1)*u+xm*xl(1,k)
      xd(2) = xd(2)*u+xm*xl(2,k)
      xd(3) = xd(3)*u+xm*xl(3,k)
12    continue
      sz = sqrt(xd(1)**2+xd(2)**2+xd(3)**2)
      if(in.eq.2.or.in.eq.-1)return
c
      do 21 id=1,3
21    xdd(id)=0
      do 22 k=ndeg,3,-1
      xm = dble((k-1)*(k-2))
      xdd(1) = xdd(1)*u+xm*xl(1,k)
      xdd(2) = xdd(2)*u+xm*xl(2,k)
      xdd(3) = xdd(3)*u+xm*xl(3,k)
22    continue
c
      do 33 id=1,3
33    xddd(id)=0.
      do 32 k=ndeg,4,-1
      xm = dble((k-1)*(k-2)*(k-3))
      xddd(1) = xddd(1)*u+xm*xl(1,k)
      xddd(2) = xddd(2)*u+xm*xl(2,k)
      xddd(3) = xddd(3)*u+xm*xl(3,k)
32    continue
c
      return
      end
C
C=END SOURCE
