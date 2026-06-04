       subroutine brasf(ndimn,x,n,q,t,cs,is1,is2,ns,index,itype)
c
c this routine tries to braket the segments of a cubic spline
c in which ( r - x)**2 may have a minimum . the search is limited to
c is1<is<is2 the segments at the extrema are checked for the presence
c of a global minimum only if they are the extrema of the whole spline
c itype=0 search only for the first segment, =1 get'em all
c
       real x(ndimn),q(ndimn,*),t(ndimn,*),cs(*)
       integer index(*)
       logical h1,h2
c
       ns = 0
       do 10 is = is1,is2
         i1 =is
         i2 =is+1
         r1 = 0.
         r2 = 0.
         hh = 0.
         cx = cs(is)
         do 15 id=1,ndimn
          r1 = r1+(x(id) - q(id,i1))*t(id,i1)
          r2 = r2+(x(id) - q(id,i2))*t(id,i2)
          r21= q(id,i1)-q(id,i2)
          t12= t(id,i1)+t(id,i2)
          hh = hh- (3*r21+cx*(t12+t(id,i1)))*
     1             (3*r21+cx*(t12+t(id,i2)))
15       continue
         hh=hh/(cx*cx)
         h1 = r1.ge.0
         h2 = r2.le.0
         if(is.eq.1.and..not.h1)then
           ns=ns+1
           index(ns)=is
         else if(is.eq.n-1.and..not.h2)then
           ns=ns+1
           index(ns)=is
         else if((h1.and.h2).or.hh.lt.-1.e-10)then
           ns=ns+1
           index(ns)=is
         endif
         if(itype.eq.0.and.ns.eq.1)return
10    continue
      return
      end
