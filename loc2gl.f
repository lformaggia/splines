C=NAME loc2gl
C=RCSfile $Source: /u/flosys3d/software/cvsroot/splines/loc2gl.f,v $
C=VERSION $Revision: 1.1.1.1 $ $Date: 1996/03/05 17:16:38 $
C=TYPE subroutine
C=PURPOSE general
C=KEYWORDS parametric_surfaces
C=BLOCK ABSTRACT
C
c It gives arc number and local coordinate from the global parametric coordinate
c and viceversa. 
C
C
C=END ABSTRACT
C=BLOCK USAGE
C
C call loc2gl(ug,u,n,i,inv,ierr)
C
C  INPUT       DIMENSION      DESCRIPTION
c
C  ug                       global coordinate (if inv =1)
c  u                        local coordinate  (if inv=-1)
c  i                        arc number        (if inv=-1)
c
C OUTPUT      DIMENSION     DESCRIPTION
c
c  u                        local coordinate (if inv=1)
C  i                        arc number       (if inv=1)
C  inv                      switch:
c                            = 1 global to local
c                            =-1 local to global 
c  n                        number of nots (n.of arcs +1)
c  ierr                     error indicator (if !=0 an error condition occurred:
c                            =1 (when inv  = 1) -> computed u less than 0
c                            =2 (when inv  = 1) -> computed u greater than 1.
c
c NOTES
c
c - limits:  0.<=u<=1.    
c 
C
C=END USAGE
C
C=BLOCK SOURCE
C
        subroutine loc2gl(ug,u,n,i,inv,ierr)
        ierr = 0
        if(inv.eq.1)then
         i = max(1,min(int(ug),n-1))
         u = ug - float(i)
         if(u.lt.0.)ierr=1
         if(u.gt.1.)ierr=2
         u = max(0.,min(1.,u))
        else
         ug = u+float(i)
        endif
        return
        end
C
C=END SOURCE
