C=NAME distps
C=RCSfile $Source: /u/flosys3d/software/cvsroot/splines/distps.f,v $
C=VERSION $Revision: 1.1.1.1 $ $Date: 1996/03/05 17:16:35 $
C=TYPE real function
C=PURPOSE Local_to locps4 Local_to xminli
C=KEYWORDS Parametric_polynimial_surfaces Line_search_Minimization
C=BLOCK ABSTRACT
c
c This routine evaluates the distance function defined as
c
c          DISTPS = || r - xp ||
c
c where xp is a point on the 3d space and r a point on a parametric surface.
c  r is here a function of a single parameter:
c
c     r = R(x) = r(u(x),v(x)) where
c
c     u(x) = ugs + x*xi(1)   v(x)  = vgs + x*xi(2)
c
c ugs,vgs are the coordinate of a starting point, while xi is a unitary vector
c which defines a direction. That is, the distance function is computed as
c is
c
c          ||R(x) -xp|| = || r( ugs+x*xi(1),vgs+x*xi(2) ) - xp||
c
c The routine may return also some derivative information
C
C
C=END ABSTRACT
C=BLOCK USAGE
C
c dist= distps(xp,ugs,vgs,xi,x,keypa,aps,n,m,z,r,der,g,ru,rv,ug,vg,ind)
C
C
c   INPUT     DIMENSION       DESCRIPTION
c
c n                         n. of knots on the surface along u-coordinate
c m                         n. of knots on the surface along v-coordinate
c aps         3,*           surface coefficent in standard format.
c keypa       (n-1)(m-1)    array of pointers to aps.
c xp          3             coordinates of given point.
c ugs,vgs                   starting value of the global parametric coordinates
c xi          2             direction
c x                         coordinate along  direction xi
c ind                         0 -> derivative information not required
c                           /=0 -> derivative information given
c
c
c   OUTPUT     DIMENSION       DESCRIPTION
c
C ug,vg                     parametric coordinates corresponding to x
c                           ug = ugs + x*xi(1) , vg=vgs + x*xi(2)
C r           3             r = r(ug,vg).
c distps                    distance = ||r(ug,vg) - xp||
c z           3             r(ug,vg) - xp
c ru          3             partial derivativo of position vector: r,u(ug,vg) (ind /=0)
c rv          3             partial derivative of position vector: r,v*ug,vg) (ind /=0)
c g           2             gradient of distance function ||r - xp|| at point (ug,vg)
c                           (ind /=0)
c der                       derivative of distance function along search direction:
c                           der = g*xi (ind /=0)
c
c  OTHERS
c
c zeps        parameter     tiny number (It prevents zero division
c
C
C=END USAGE
C
C=BLOCK SOURCE
C
      real*8 function distps(xp,ugs,vgs,xi,beta,keypa,aps,n,m,z,r,der,
     1                       g,ru,rv,ug,vg,ind)
      implicit none
      integer*8 i,id,ind,j,keypa(*),kpa,m,n,ndu,ndv
      real*8 aps(3,*),beta,der,g(2),r(3),ru(3),ruu(3),ruv(3),rv(3)
      real*8 rvv(3),u,ug,ugs,v,vg,vgs,xi(2),xp(3),z(3),xxx,zeps
      parameter(zeps=1.d-32)
c
      ug = ugs + beta*xi(1)
      vg = vgs + beta*xi(2)
c
      i = min(max(1,int(ug)),n-1)
      j = min(max(1,int(vg)),m-1)
      u = ug - dble(i)
      v = vg - dble(j)
c
c make sure u,v,ug,and vg are within correct bounds
c
      u = max(min(u,1.d0),0.d0)
      v = max(min(v,1.d0),0.d0)
      ug = dble(i) + u
      vg = dble(j) + v
c
      kpa = keypa((i-1)*(m-1)+j)
      ndu = int(aps(1,kpa))
      ndv = int(aps(2,kpa))
      if(ind.eq.0)then
        id=1
      else
        id=2
      endif
c
      call evsurg(aps(1,kpa),u,v,ndu,ndv,r,ru,rv,ruv,ruu,rvv,id)
c
      z(1) = r(1) -xp(1)
      z(2) = r(2) -xp(2)
      z(3) = r(3) -xp(3)
      distps=sqrt(z(1)*z(1)+z(2)*z(2)+z(3)*z(3))
c
      if(ind.eq.0)return
c
c evaluate derivative along xi
c
      g(1) = 2.d0*(z(1)*ru(1)+z(2)*ru(2)+z(3)*ru(3))
      g(2) = 2.d0*(z(1)*rv(1)+z(2)*rv(2)+z(3)*rv(3))
      xxx  = 1.d0/max(distps,zeps)
      g(1) = xxx*g(1)
      g(2) = xxx*g(2)
      der  = g(1)*xi(1) + g(2)*xi(2)
      return
      end
C
C=END SOURCE
