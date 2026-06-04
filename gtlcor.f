C=NAME gtlcor
C=RCSfile $Source: /u/flosys3d/software/cvsroot/splines/gtlcor.f,v $
C=VERSION $Revision: 1.1.1.1 $ $Date: 1996/03/05 17:16:38 $
C=TYPE subroutine
C=PURPOSE general
C=KEYWORDS parametric_surfaces
C=BLOCK ABSTRACT
C
c It gives patch number and local coordinates from the global parametric coordinates
c and viceversa. It works also for curves
C
C
C=END ABSTRACT
C=BLOCK USAGE
C
C call gtlcor(ug,vg,u,v,ipa,i,j,n,m,inv,ierr)
C
C  INPUT       DIMENSION      DESCRIPTION
c
C  ug,vg                    global coordinates (if inv =1)
c  u,v                      local coordinates  (if inv=-1)
c  ipa                      patch number       (if inv=-1)
c  n,m                      n. of patches +1 along u and v respectively (if inv=-1)
c
C OUTPUT      DIMENSION     DESCRIPTION
c
c  u ,v                     local coordinates (if inv=1)
C  ipa                      patch number      (if inv=1)
c  n,m                      n. of patches +1 along u and v respectively (if inv=1)
C  inv                      switch:
c                            = 1 global to local
c                            =-1 local to global 
c  i,j                      position of patch in the lattice
c
c  ierr                     error indicator (if !=0 an error condition occurred:
c                            =1 (when inv  = 1) -> computed u or v less than 0
c                            =2 (when inv  = 1) -> computed u or v greater than 1.
c                            =1 (when inv  =-1) -> input patch number greater than 
c                                                  (n-1)*(m-1)
C
c NOTES
c
c - limits:  0.<=u<=1.  0.<=v<=1.     0 <ipa<=(n-1)*(m-1)   
c            1.<=ug<=n  1.<=vg<=m-1
c
c - ipa gives the absolute patch number. To evaluate patch position in the (i,j) lattice:
c
c     j = mod(ipa-1,m)+1   i = 1 + (ipa-j)/m (BEWARE! : it is an integer division)
c
c - the opposite operation:
c
c     ipa = j + (i-1)*m
c
c 
C
C=END USAGE
C
C=BLOCK SOURCE
C
        subroutine gtlcor(ug,vg,u,v,ipa,i,j,n,m,inv,ierr)
        implicit none
        integer*8 i,ierr,inv,ipa,j,m,n
        real*8 u,ug,v,vg
        ierr = 0
        if(inv.ne.-1)then
         i = max(1,min(int(ug),n-1))
         j = max(1,min(int(vg),m-1))
         u = ug - dble(i)
         v = vg - dble(j)
         if(u.lt.0.d0.or.v.lt.0.d0)ierr=1
         if(u.gt.1.d0.or.v.gt.1.d0)ierr=2
         u = max(0.d0,min(1.d0,u))
         v = max(0.d0,min(1.d0,v))
         ipa = j + (i-1)*m
        else
         j = mod(ipa-1,m) +1
         i = 1 + (ipa - j)/m
         if(i.gt.n-1)ierr=1
         i = min(i,n-1)
         ug = u+dble(i)
         vg = v+dble(j)
        endif
        return
        end
C
C=END SOURCE
