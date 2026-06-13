      subroutine LOCSU2(n,m,ndimn,tanret,choret,coor,x,
     1                  ug ,vg ,u10,u20,v10,v20,r,z,dist,imes)
      implicit none
      integer*8 i,i0,id,imes,ispa,iter,j,j0,m,miter,n,ndimn
      real*8 a,apatch(5,4,4),aps(1),b,beta,c,choret(2,n,m)
      real*8 coor(ndimn,*),da,da1,det,dist,du,dv,eps1,eps2,eps3,gg(2)
      real*8 h(2),h0(2),r(*),ru(5),ruu(5),ruv(5),rv(5)
      real*8 rvv(5),tanret(ndimn,3,n,m),u,u1,u10,u2,u20,ug,ug0
      real*8 v,v1,v10,v2,v20,vg,vg0,x(*),xd,xgr1,xgr2,xl,z(*)
      parameter(eps1=0.001d0,eps2=0.00001d0,eps3=0.0001d0)
      parameter(miter=100)
c
c This routine locates the global coordinates ugs and vgs in the
c parametric plane which correspond to the point in the bicubic
c surface which minimize the Euclidian norm
c
c                    ||X - R||
c
c where X is a fixed point. It uses a Newton procedure
cThe search is limited in the
c presrcibed interval [(u1,u2),(v1,v2)]
c
c DIST = distance                    R = found point
c UG,VG= global parametric coords    Z = (R - X)
c IMES = message :
c                 0 point R is on the surface (||r-x||<eps1*xl)
c                   where xl is an average patch linear dimension
c                -1 minimum found, but point not on the spline
c                -2 extrema found
c                -3 mAX ITERATION EXCEEDED
c
      common/surpa/ispa,aps
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
        iter=0
c
1000  continue
c
c start Newton iterations
c
      if(iter.gt.miter)then
        imes=-3
        return
      endif
c
      iter=iter+1
      if(i0.ne.i.or.j0.ne.j)then
         call evapa2(n,m,ndimn,tanret,choret,coor,apatch,i,j)
         i0=i
         j0=j
         if(ispa.eq.0)then
         xl=max(choret(1,i,j),choret(2,i,j))
      else
         xl =0.d0
         do id =1,ndimn
           xl = xl + (apatch(id,1,1)-
     1         (apatch(id,1,1)+apatch(id,2,1)+apatch(id,3,1)+
     1          apatch(id,4,1)))**2
         end do
         xl =sqrt(xl)
      endif
      endif
c
      call gpsur(ndimn,u,v,r,apatch)
c
      dist = 0.d0
      do id=1,ndimn
        z(id) = r(id) - x(id)
        dist = dist + z(id)*z(id)
      end do
      dist = sqrt(dist)
      if(dist.lt.eps1*xl)return
c
      call gders1(ndimn,ru,rv,ruv,ruu,rvv,apatch,u,v)
      xgr1 = 0.d0
      xgr2 = 0.d0
      do id=1,ndimn
        xgr1 = xgr1 + z(id)*ru(id)
        xgr2 = xgr2 + z(id)*rv(id)
      end do
c
      beta = sqrt(xgr1**2+xgr2**2)
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
c                      | a   b |
c get the Hessian  H = |       |
c                      | b   c |
      a = 0.d0
      b = 0.d0
      c = 0.d0
      do id=1,ndimn
        a = a + ru(id)*ru(id)
        b = b + ru(id)*rv(id)
        c = c + rv(id)*rv(id)
      end do
      det = 1.d0/(a*c-b*b)
      du = det*(c*xgr1 - b*xgr2)
      dv = det*(a*xgr2 - b*xgr1)
c
c  get the step da and limit it (max(da) =1)
c
      da = sqrt(du*du+dv*dv)
      da1=1.d0/max(0.0001d0,da)
      du = du*da1
      dv = dv*da1
      da = min(1.d0,da)
c
c get new point
c
      ug=min(u2,max(u1,ug-da*du))
      vg=min(v2,max(v1,vg-da*dv))
      i =min(int(ug),n-1)
      j =min(int(vg),m-1)
      u =ug -dble(i)
      v =vg -dble(j)
      go to 1000
      return
      end
