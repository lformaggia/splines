C=NAME evsurg
C=RCSfile $Source: /u/flosys3d/software/cvsroot/splines/evsurg.f,v $
C=VERSION $Revision: 1.1.1.1 $ $Date: 1996/03/05 17:16:36 $
C=TYPE subroutine
C=PURPOSE general
C=KEYWORDS parametric_Polynomial_surfaces
C=BLOCK ABSTRACT
C
C It computes  position vector and derivatives on a parametric surface in
c the tridimensional space, given the patch coefficents and the local parametric
c coordinates
c
C=END ABSTRACT
C=BLOCK USAGE
C
C  call evsurg(aps,u,v,ndu,ndv,r,ru,rv,ruv,ruu,rvv,ind)
C
C  INPUT       DIMENSION      DESCRIPTION
c
c aps         3,*           array of patch coefficents
c u,v                       local parametric coordinates
c ndu,ndu                   degree of polinomial +1 in the u and v direction
c                           repectively (i.e. ndu=5 means that the highest power of u
c                           in the polinomial is u**4)
c ind                       switch governing the kind of output given
c                           -1 -> only ru and rv
c                            0 -> all
c                            1 -> only position vector r
c                            2 -> only r and 1st derivs r,u and r,v
c
C OUTPUT      DIMENSION     DESCRIPTION
c
c r           3             position vector
c ru,rv       3             1st derivatives   r,u r,v
c ruu,rvv,ruv 3             2nd derivatives   r,uu r,vv r,uv
c
c NOTE:
c     Though the value of ndu and ndv should be already contained in the first two
c     locations of aps (following the convention described in the User's Reference Manual),
c     they must be given as input. This for compatibility with an old routine. It may
c     change in a future release, where ndu and ndv may be output quantities.
c
C=END USAGE
C
C=BLOCK SOURCE
C
      subroutine evsurg(aps,u,v,ndu,ndv,r,ru,rv,ruv,ruu,rvv,ind)
c
c  it evaluates point position and derivatives corresponding to
c  the parametric coordinates u,v of a surface patch whose coefficents
c  are stored in aps, the first record of apa contains ndu,ndv,xl and
c  it is skipped as we assume that ndu and ndv are already known.
c
c  ndu,ndv = degree of patch +1  on u and v directions repectively
c  in      = indicator :
c                       -1 only ru and rv
c                        0 all
c                        1 only r
c                        2 only r,ru,rv
c
      real aps(3,0:*)
      real r(*),ru(*),rv(*),ruv(*),rvv(*),ruu(*)
      double precision a1d,a2d,a3d,ud,vd,r1d,r2d,r3d
c
      apatch(id,i,j) = aps(id,(j-1)*ndu+i)
c
      ud = dble(u)
      vd = dble(v)
      if(ind.eq.-1)go to 100
      r1d=0.d00
      r2d=0.d00
      r3d=0.d00
      do 10 j=ndv,1,-1
          a1d=0.d00
          a2d=0.d00
          a3d=0.d00
          do 11 i=ndu,1,-1
             a1d  = a1d*ud + dble(apatch(1,i,j))
             a2d  = a2d*ud + dble(apatch(2,i,j))
             a3d  = a3d*ud + dble(apatch(3,i,j))
11        continue
          r1d = r1d*vd + a1d
          r2d = r2d*vd + a2d
          r3d = r3d*vd + a3d
10     continue
       r(1) = r1d
       r(2) = r2d
       r(3) = r3d
100    continue
       if(ind.eq.1)return
c
      ru(1)=0.
      ru(2)=0.
      ru(3)=0.
      rv(1)=0.
      rv(2)=0.
      rv(3)=0.
c
      do 20 j=ndv,1,-1
          a1=0.
          a2=0.
          a3=0.
          do 21 i=ndu,2,-1
             xm = (i-1)
             a1 = a1*u + apatch(1,i,j)*xm
             a2 = a2*u + apatch(2,i,j)*xm
             a3 = a3*u + apatch(3,i,j)*xm
21        continue
          ru(1) = ru(1)*v +a1
          ru(2) = ru(2)*v +a2
          ru(3) = ru(3)*v +a3
20     continue
c
      do 30 j=ndv,2,-1
          a1=0.
          a2=0.
          a3=0.
          do 31 i=ndu,1,-1
             a1 = a1*u + apatch(1,i,j)
             a2 = a2*u + apatch(2,i,j)
             a3 = a3*u + apatch(3,i,j)
31        continue
          xm = (j-1)
          rv(1) = rv(1)*v +a1*xm
          rv(2) = rv(2)*v +a2*xm
          rv(3) = rv(3)*v +a3*xm
30     continue
c
       if(ind.eq.2.or.ind.eq.-1 )return
       ruv(1)=0.
       ruv(2)=0.
       ruv(3)=0.
       ruu(1)=0.
       ruu(2)=0.
       ruu(3)=0.
       rvv(1)=0.
       rvv(2)=0.
       rvv(3)=0.
      do 40 j=ndv,2,-1
          a1=0.
          a2=0.
          a3=0.
          do 41 i=ndu,2,-1
             xm=(i-1)
             a1 = a1*u + apatch(1,i,j)*xm
             a2 = a2*u + apatch(2,i,j)*xm
             a3 = a3*u + apatch(3,i,j)*xm
41        continue
          xm = (j-1)
          ruv(1) = ruv(1)*v +a1*xm
          ruv(2) = ruv(2)*v +a2*xm
          ruv(3) = ruv(3)*v +a3*xm
40     continue
c
      do 50 j=ndv,1,-1
          a1=0.
          a2=0.
          a3=0.
          do 51 i=ndu,3,-1
             xm=(i-1)*(i-2)
             a1 = a1*u + apatch(1,i,j)*xm
             a2 = a2*u + apatch(2,i,j)*xm
             a3 = a3*u + apatch(3,i,j)*xm
51        continue
          ruu(1) = ruu(1)*v +a1
          ruu(2) = ruu(2)*v +a2
          ruu(3) = ruu(3)*v +a3
50     continue
c
      do 60 j=ndv,3,-1
          a1=0.
          a2=0.
          a3=0.
          do 61 i=ndu,1,-1
             a1 = a1*u + apatch(1,i,j)
             a2 = a2*u + apatch(2,i,j)
             a3 = a3*u + apatch(3,i,j)
61        continue
          xm = (j-1)*(j-2)
          rvv(1) = rvv(1)*v +a1*xm
          rvv(2) = rvv(2)*v +a2*xm
          rvv(3) = rvv(3)*v +a3*xm
60     continue
       return
       end
C
C=END SOURCE
