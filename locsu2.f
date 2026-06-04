      subroutine LOCSU2(n,m,ndimn,tanret,choret,coor,x,
     1                  ug ,vg ,u10,u20,v10,v20,r,z,dist,imes)
      parameter(eps1=0.001,eps2=0.00001,eps3=0.0001)
      parameter(mdimn=5)
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
       real tanret(ndimn,3,n,m),choret(2,n,m),x(*),r(*),z(*)
       real coor(ndimn,*), apatch(mdimn,4,4),h0(2),gg(2)
       real h(2),ru(mdimn),rv(mdimn),ruu(mdimn),ruv(mdimn),rvv(mdimn)
c
c
       imes = 0
       u1   = max(1.,u10)
       v1   = max(1.,v10)
       u2   = min(float(n),u20)
       v2   = min(float(m),v20)
       ug  = min(u2 ,max(u1,ug ))
       vg  = min(v2 ,max(v1,vg ))
       i   = min(int(ug ),n-1)
       j   = min(int(vg ),m-1)
       u   = ug  -i
       v   = vg  -j
       ug0 = -10
       vg0 = -10
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
         xl =0.
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
      dist = 0
      do 10 id=1,ndimn
        z(id) = r(id) - x(id)
10      dist = dist + z(id)*z(id)
      dist = sqrt(dist)
      if(dist.lt.eps1*xl)return
c
      call gders1(ndimn,ru,rv,ruv,ruu,rvv,apatch,u,v)
      xgr1 = 0.
      xgr2 = 0.
      do 20 id=1,ndimn
        xgr1 = xgr1 + z(id)*ru(id)
        xgr2 = xgr2 + z(id)*rv(id)
20    continue
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
      a = 0
      b = 0
      c = 0
      do 90 id=1,ndimn
        a = a + ru(id)*ru(id)
        b = b + ru(id)*rv(id)
        c = c + rv(id)*rv(id)
90    continue
      det = 1./(a*c-b*b)
      du = det*(c*xgr1 - b*xgr2)
      dv = det*(a*xgr2 - b*xgr1)
c
c  get the step da and limit it (max(da) =1)
c
      da = sqrt(du*du+dv*dv)
      da1=1./max(0.0001,da)
      du = du*da1
      dv = dv*da1
      da = min(1.,da)
c
c get new point
c
      ug=min(u2,max(u1,ug-da*du))
      vg=min(v2,max(v1,vg-da*dv))
      i =min(int(ug),n-1)
      j =min(int(vg),m-1)
      u =ug -i
      v =vg -j
      go to 1000
      return
      end
