       subroutine locpfs(ndimn,n,x,q,cs,t,ugs,u1,u2,ug,r,v,dist,imes)
       parameter(eps1=0.001,mit=6,mdimn=3)
       parameter(eps2=0.001,eps3=0.001)
       real cs(*),q(ndimn,*),t(ndimn,*),x(ndimn)
       real h,r(mdimn),rd(mdimn),rdd(mdimn),rddd(mdimn)
       real xnn(mdimn),xnb(mdimn),v(mdimn)
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
        u   = ugs - is
       ug  = ugs
       ug0   = -10
c
c start Newton iterations
c
1000   continue
       xl = cs(is)
c
c get point and gradients
c
       call getp2(ndimn,n,q,t,cs,is,u,r,rd,rdd,rddd,s,0)
       h = 0
       dist = 0
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
       da = 0
c      do 30 nit = 1,mit
       xden =0.
       xpp  =0.
       do 40 id=1,ndimn
          xden = xden + rd(id)*rd(id)
          xpp  = xpp  + v(id)*rdd(id)
40     continue
       if((xpp+xden).gt.eps3*xl*xl)xden = xpp+xden
       da =-beta/xden
c
c dont let ug step over more than 1 segment
c
       if(abs(da).gt.1.)da = sign(1.,da)
       ug = ug +da*h
       ug = min(u2      ,max(u1,ug ))
       is = min(int(ug),n-1)
       u  = ug - is
       go to 1000
       return
       end
