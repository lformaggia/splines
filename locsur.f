      subroutine LOCSUR(n,m,ndimn,tanret,choret,coor,x,
     1                  ug ,vg ,u10,u20,v10,v20,r,z,dist,imes)
      implicit none
      integer*8 i,i0,id,imes,ispa,j,j0,m,n,ndimn
      real*8 apatch(5,4,4),aps(1),beta,choret(2,n,m),coor(ndimn,*)
      real*8 da,dist,eps1,eps2,eps3,h(2),h11,h12,h22,r(*),ru(5)
      real*8 ruu(5),ruv(5),rv(5),rvv(5),tanret(ndimn,3,n,m)
      real*8 u,u1,u10,u2,u20,ug,ug0,v,v1,v10,v2,v20,vg,vg0,x(*),xd
      real*8 xden,xl,xpp,z(*)
      parameter(eps1=0.001d0,eps2=0.00001d0,eps3=0.0001d0)
      common/surpa/ispa,aps
c
c This routine locates the global coordinates ugs and vgs in the
c parametric plane which correspond to the point in the bicubic
c surface which minimize the Euclidian norm
c
c                    ||X - R||
c
c where X is a fixed point. It uses a Newton-Rapson type procedure
c associated with a steepest descent search. The search is limited in the
c prescibed interval [(u1,u2),(v1,v2)]
c
c DIST = distance                    R = found point
c UG,VG= global parametric coords    Z = (R - X)
c IMES = message :
c                 0 point R is on the surface (||r-x||<eps1*xl)
c                   where xl is an average patch linear dimension
c                -1 minimum found, but point not on the spline
c                -2 extrema found
c
       imes = 0
       u1   = max(1.d0,u10)
       v1   = max(1.d0,v10)
       u2   = min(dble(n),u20)
       v2   = min(dble(m),v20)
       ug  = min(u2 ,max(u1,ug ))
       vg  = min(v2 ,max(v1,vg ))
       i   = min(int(ug ),n-1)
       j   = min(int(vg ),m-1)
       u   = ug  -dble(i)
       v   = vg  -dble(j)
       ug0 = -10.d0
       vg0 = -10.d0
       i0 =-10
       j0 =-10
c
1000  continue
c
c start Newton iterations
c
      if(i0.ne.i.or.j0.ne.j)then
         call evapa2(n,m,ndimn,tanret,choret,coor,apatch,i,j)
         i0=i
         j0=j
         if(ispa.eq.0)then
         xl=max(choret(1,i,j),choret(2,i,j))
         else
         xl =0.d0
         do 143 id =1,ndimn
143        xl = xl + (apatch(id,1,1)-
     1         (apatch(id,1,1)+apatch(id,2,1)+apatch(id,3,1)+
     1          apatch(id,4,1)))**2
         xl =sqrt(xl)
      endif
      endif
c
      call gpsur(ndimn,u,v,r,apatch)
c
      dist = 0.d0
      do 10 id=1,ndimn
        z(id) = r(id) - x(id)
10      dist = dist + z(id)*z(id)
      dist = sqrt(dist)
      if(dist.lt.eps1*xl)return
c
c get search direction h = grad(||r-x||**2)
c
      call gders1(ndimn,ru,rv,ruv,ruu,rvv,apatch,u,v)
      h(1) = 0.d0
      h(2) = 0.d0
      do 20 id=1,ndimn
        h(1) = h(1) + z(id)*ru(id)
        h(2) = h(2) + z(id)*rv(id)
20    continue
c
      beta = sqrt(h(1)**2+h(2)**2)
c
c have I reached a minimum?
c
      xd = max(dist,xl)
      if(beta.lt.xd*eps1)then
         imes =-1
         return
      endif
c
c have I reached an extrema ?
c
      if(abs(ug0-ug).lt.eps2.and.abs(vg0-vg).lt.eps2)then
         imes = -2
         return
      endif
      ug0 =ug
      vg0 =vg
      h(1) = h(1)/beta
      h(2) = h(2)/beta
c
c  get the step da
c
      xden = 0.d0
      xpp  = 0.d0
      h11  = h(1)*h(1)
      h22  = h(2)*h(2)
      h12  = h(1)*h(2)
      do 50 id=1,ndimn
       xden = xden + h11*ru(id)*ru(id)+2.d0*h12*ru(id)*rv(id)+
     1               h22*rv(id)*rv(id)
       xpp  = xpp  +(h11*ruu(id)+h22*rvv(id)+2.d0*h12*ruv(id))*
     1               z(id)
50    continue
      if(xpp+xden.gt.eps3*xl*xl)xden = xden+xpp
      da = -beta/xden
c
c limit step
c
      if(abs(da).gt.1.d0)da=sign(1.d0,da)
      ug=min(u2,max(u1,ug+da*h(1)))
      vg=min(v2,max(v1,vg+da*h(2)))
      i =min(int(ug),n-1)
      j =min(int(vg),m-1)
      u =ug -dble(i)
      v =vg -dble(j)
      go to 1000
      return
      end
