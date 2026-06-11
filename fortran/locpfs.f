      subroutine locpfs(ndimn,n,x,q,cs,t,ugs,u1,u2,ug,r,v,dist,imes)
      implicit none
      integer*8 id,imes,is,n,ndimn
      real*8 beta,cs(*),da,dist,eps1,eps2,eps3,h,q(ndimn,*),r(3)
      real*8 rdd(3),rddd(3),rd(3),s,t(ndimn,*),u,u1,u2
      real*8 ug,ug0,ugs,v(3),x(ndimn),xden,xd,xnb(3),xnn(3)
      real*8 xl,xpp
      integer*8 mit
      parameter(eps1=0.001d0,mit=6)
      parameter(eps2=0.001d0,eps3=0.001d0)
c
c this routine evaluates the point r  on a ferguson spline
c which is nearest to the fixed point x
c
c ugs  = global u coordinate of s     dist = distance
c        the point from which the     ug   = global u coord.
c        the search has to start             of found point
c        ugst = ist + u ,>n1,<n2      u1,u2= global u coord of
c                                            search limits
c imes= message giving the results        v = r - x
c       of the search:
c       imes = 0 point r on the spline
c       imes =-1 minimum found (point not on spline)
c       imes =-2 extrema reached
c
       imes = 0
       is   =min(int(ugs),n-1)
        u   = ugs - dble(is)
       ug  = ugs
       ug0   = -10.d0
c
c start Newton iterations
c
1000   continue
       xl = cs(is)
c
c get point and gradients
c
       call getp2(ndimn,n,q,t,cs,is,u,r,rd,rdd,rddd,s,0)
       h = 0.d0
       dist = 0.d0
       do 10 id=1,ndimn
        v(id) = (r(id) - x(id))
         h = h +   v(id)*rd(id)
         dist = dist + v(id)*v(id)
10     continue
       dist = sqrt(dist)
       if(dist .le.eps1*xl)return
c
c point found on the spline - imes =0
c
       beta = abs(h)
c
        xd=max(xl,dist)*eps1
       if(abs(h).lt.xd)then
c
c gradient null at the point, maximum reached, imes=-1
c
        imes = -1
        return
       endif
       if(abs(ug0-ug).lt.eps2)then
c
c I have reached an extrema imes=-2
c
           imes =-2
           return
       endif
       ug0= ug
       h = h/beta
       da = 0.d0
c      do 30 nit = 1,mit
       xden =0.d0
       xpp  =0.d0
       do 40 id=1,ndimn
          xden = xden + rd(id)*rd(id)
          xpp  = xpp  + v(id)*rdd(id)
40     continue
       if((xpp+xden).gt.eps3*xl*xl)xden = xpp+xden
       da =-beta/xden
c
c dont let ug step over more than 1 segment
c
       if(abs(da).gt.1.d0)da = sign(1.d0,da)
       ug = ug +da*h
       ug = min(u2      ,max(u1,ug ))
       is = min(int(ug),n-1)
       u  = ug - dble(is)
       go to 1000
       return
       end
