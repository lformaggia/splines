C=NAME brauv
C=RCSfile $Source: /u/flosys3d/software/cvsroot/splines/brauv.f,v $
C=VERSION $Revision: 1.1.1.1 $ $Date: 1996/03/05 17:16:34 $
C=TYPE subroutine
C=PURPOSE Local_to xminli
C=KEYWORDS Parametric_polynomial_surfaces Line_search_Minimization
C=BLOCK ABSTRACT
c
c This routine is used by the real function xminli to bracket a minumum, i.e. to
c find an interval where a local minimum exist, in order to sucessively operate
c with the standard algorithm. The distance function
c indicated in the followin as F, and  is defined by:
c
c   F(x) = ||R(x) -xp|| = || r(ugs+x*xi(1),vgs+x*xi(2)) - xp||
c
C=END ABSTRACT
C=BLOCK USAGE
C
c call brauv(ax,bx,cx,fa,fb,fc,ugs,vgs,keypa,aps,n,m,xp,xi,bl1,bl2)
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
c xi          2             search direction
c ax,bx                     starting search interval
c bl1,bl2                   interval into which minimum is seeked, i.e. we look
c                           for the mimimum of the distance function in the
c                           interval [bl1,bl2]. The interval may be computed by
c                           using the routine limiuv
c
c   OUTPUT     DIMENSION       DESCRIPTION
c
c ax,bx,cx                  bracketing absissae, if F indicates the distance function we
c                           have F(ax)<=F(bx)<=F(cx) and the search intervall will be
c                           the one spanned by ax,bx,cx
c fa,fb,fc                  value of distance function at ax,bx,cx repectively
c
c
C
C=END USAGE
C
C=BLOCK SOURCE
C
      subroutine brauv(ax,bx,cx,fa,fb,fc,ug,vg,keypa,aps,n,m,
     1                 xp,xi,bl1,bl2)
      parameter (gold=1.618034, glimit=100., tiny=1.e-20)
      real aps(3,*),xi(2),rr(3),xp(3),ru(2),rv(2)
      real z(3),g(2)
      integer keypa(*)
      fa=distps(xp,ug,vg,xi,ax,keypa,aps,n,m,z,rr,der,g,
     1                     ru,rv,ugn,vgn,0)
      fb=distps(xp,ug,vg,xi,bx,keypa,aps,n,m,z,rr,der,g,
     1                     ru,rv,ugn,vgn,0)
c
77    if(fb.gt.fa)then
        dum=ax
        ax=bx
        bx=dum
        dum=fb
        fb=fa
        fa=dum
      endif
      cx=bx+gold*(bx-ax)
      if(cx.le.bl1.or.cx.ge.bl2)then
         cx = (bx+gold*ax)/(1.+gold)
         fc=distps(xp,ug,vg,xi,cx,keypa,aps,n,m,z,rr,der,g,
     1                     ru,rv,ugn,vgn,0)
c
         if(abs(ax-bx).le.tiny)goto 2
c
         if(fc.lt.fa)then
           dum = bx
           bx  = cx
           cx  = dum
           dum = fb
           fb  = fc
           fc  = dum
         else
           if(abs(ax-cx).le.tiny)goto 2
           ax = cx
           fa = fc
           goto 77
         endif
      else
        fc=distps(xp,ug,vg,xi,cx,keypa,aps,n,m,z,rr,der,g,
     1                     ru,rv,ugn,vgn,0)
      endif
1     if(abs(ax-cx).le.tiny)goto 2
c modification from numerical recipes routine : .gt. instead of .ge.
c      if(fb.gt.fc)then
      if(fb.ge.fc.and.(bx.ne.cx))then
c
        r=(bx-ax)*(fb-fc)
        q=(bx-cx)*(fb-fa)
        diff = q-r
        if(diff.ge.0)then
           diff=max(tiny,diff)
        else
           diff = min(diff,-tiny)
        endif
        u=bx-((bx-cx)*q-(bx-ax)*r)/(2.*sign(max(abs(q-r),tiny),diff))
        ulim=bx+glimit*(cx-bx)
        ulim= max(bl1,min(ulim,bl2))
        if((bx-u)*(u-cx).gt.0.)then
          fu=distps(xp,ug,vg,xi,u,keypa,aps,n,m,z,rr,der,g,
     1                     ru,rv,ugn,vgn,0)
          if(fu.lt.fc)then
            ax=bx
            fa=fb
            bx=u
            fb=fu
            go to 1
          else if(fu.gt.fb)then
            cx=u
            fc=fu
            go to 1
          endif
          u=cx+gold*(cx-bx)
          u = min(bl2,max(bl1,u))
          fu=distps(xp,ug,vg,xi,u,keypa,aps,n,m,z,rr,der,g,
     1                     ru,rv,ugn,vgn,0)
        else if((cx-u)*(u-ulim).gt.0.)then
          fu=distps(xp,ug,vg,xi,u,keypa,aps,n,m,z,rr,der,g,
     1                     ru,rv,ugn,vgn,0)
          if(fu.lt.fc)then
            bx=cx
            cx=u
            u=cx+gold*(cx-bx)
            u=min(bl2,max(bl1,u))
            fb=fc
            fc=fu
            fu=distps(xp,ug,vg,xi,u,keypa,aps,n,m,z,rr,der,g,
     1                     ru,rv,ugn,vgn,0)
          endif
        else if((u-ulim)*(ulim-cx).ge.0.)then
          u=ulim
          fu=distps(xp,ug,vg,xi,u,keypa,aps,n,m,z,rr,der,g,
     1                     ru,rv,ugn,vgn,0)
        else
          u=cx+gold*(cx-bx)
          u = min(bl2,max(bl1,u))
          fu=distps(xp,ug,vg,xi,u,keypa,aps,n,m,z,rr,der,g,
     1                     ru,rv,ugn,vgn,0)
        endif
        ax=bx
        bx=cx
        cx=u
        fa=fb
        fb=fc
        fc=fu
        go to 1
      endif
2     return
      end
C
C=END SOURCE
