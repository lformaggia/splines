C=NAME getk2
C=RCSfile $Source: /u/flosys3d/software/cvsroot/splines/getk2.f,v $
C=VERSION $Revision: 1.1.1.1 $ $Date: 1996/03/05 17:16:37 $
C=TYPE subroutine
C=PURPOSE general
C=KEYWORDS parametric_curves
C=BLOCK ABSTRACT
C
C It computes the curvature,torsion, normal and binormal
c vectors at a point on a curve from the value of the derivatives
c at the point
C
C
C=END ABSTRACT
C=BLOCK USAGE
C
C call getk2(ndimn,s,xu,xuu,xuuu,xk,tau,xnd,xnb)
C
C
C  INPUT       DIMENSION      DESCRIPTION
c
c ndimn                     n. of dimensions (IT MUST BE equal to 2 or 3 !!)
c xu          ndimn         x,u 1st derivative with respect to parameter u at the
c                           point
c xuu         ndimn         x,uu 2nd derivative
c xuuu        ndimn         x,uuu 3rd derivative
c
C OUTPUT      DIMENSION     DESCRIPTION
c
c s                         ||x,u||
c xk                        curvature
c tau                       torsion
c xnd         ndimn         normal versor (||xnd||=1)
c xnb         ndimn         binormal      (||xnb||=1)
c
c NOTE:
c
c  ndimn must be equal either to 2 or 3.
c  The values of the various quantities are computed by using Frenet-Serret formulae.
c
C=END USAGE
C
C=BLOCK SOURCE
C
      subroutine getk2(ndimn,s,xd,xdd,xddd,xk,tau,xnd,xnb)
      implicit none
      integer*8 i,ndimn
      real*8 help(3),h2(3),rhe,s,s12,s2,tau,tau2,tiny,vbig,vsmall
      real*8 xd(*),xd2(3),xdd(*),xdd2(3),xddd(*),xddd2(3),xk,xk2
      real*8 xnb(*),xnb2(3),xnd(*),xnd2(3),xk12,xnx,xny,xnz
      parameter(tiny=1.d-10)
      parameter(VBIG=1.d+28,VSMALL=1.d-28)
c
c this routine gets all the other parameters that can be needed
c
c  xk = curvature   tau = torsion    xnd =normal direction
c  xnb = binormal direction
c
c  nb.: xd contains the derivative not the tangent unitary vector!!!!!
c
c Put into dble
c
      s2=0.d0
      do 67 i=1,ndimn
         xd2(i)   = dble(xd(i))       
         xdd2(i)  = dble(xdd(i))       
         xddd2(i) = dble(xddd(i))
         s2       = s2 + xd2(i)*xd2(i)
 67   continue
c
      s2 = sqrt(s2)
      s  = s2
      s12=1.d00/max(s2,VSMALL)
c
c evaluate  ru * ruu and store tangent vector 
c
      if(ndimn.eq.3)then
         help(1) = xd2(2)*xdd2(3)-xd2(3)*xdd2(2)
         help(2) = xd2(3)*xdd2(1)-xd2(1)*xdd2(3)
         help(3) = xd2(1)*xdd2(2)-xd2(2)*xdd2(1)
         h2  (1) = xd2(1)*s12
         h2  (2) = xd2(2)*s12
         h2  (3) = xd2(3)*s12
      else if(ndimn.eq.2)then
         help(1) = 0.d00
         help(2) = 0.d00
         help(3) = xd2(1)*xdd2(2)-xd2(2)*xdd2(1)
         h2  (1) = xd2(1)*s12
         h2  (2) = xd2(2)*s12
         h2  (3) = 0.d00
      endif
      xk2 = 0.d0
      do 10 i=1,3
         xk2 = xk2 + help(i)*help(i)
 10   continue
      xk2 = sqrt(xk2)
c
c  normalize ru*ruu  = b (binormal vector)
c
      if(ndimn.eq.2)then
         if(xk2.lt.VSMALL)then
            help(3)=1.0d00
         else
            help(3)=sign(1.d00,help(3))
         endif
      else
         if(xk2.lt.max(tiny*s2*s2,VSMALL))then
c
c Straight line: normal and binormal vectors are not uniquely
c defined. 
c
            rhe = sqrt(h2(3)*h2(3)+h2(2)*h2(2))
            if(xk2.gt.0.57d0)then
               help(1)= 0.d00
               help(2)=-h2(3)
               help(3)= h2(2)
            else
               rhe = sqrt(h2(1)*h2(1)+h2(2)*h2(2))
               help(1)=-h2(2)
               help(2)= h2(1)
               help(3)= 0.d00
            endif
         else
            rhe=xk2
         Endif
         rhe=1.d0/rhe
         do 15 i=1,3
            help(i) = help(i)*rhe
 15      continue
      endif
c
c   xk = |ru * ruu|/|ru|^3
c
      xk2 = xk2*s12*s12*s12
      xk  = xk2
c
c  n = b * t  (normal vector)
c
      xnx   = (help(2)*h2(3)-help(3)*h2(2))
      xny   = (help(3)*h2(1)-help(1)*h2(3)) 
      xnz   = (help(1)*h2(2)-help(2)*h2(1))
      h2(1) = xnx
      h2(2) = xny
      h2(3) = xnz
c
      do 50 i=1,ndimn
         xnb(i) = help(i)
         xnd(i) = h2(i)
 50   continue
c
c  torsion tau = (ru *ruu).ruuu/s**6*xk**2 = bn . ruuu/s**3*xk
c
c
      if(ndimn.eq.2)then
         tau=0.d0
      else
         tau2 = xddd2(1)*help(1)+help(2)*xddd2(2)+help(3)*xddd2(3)
         tau  = tau2*s12*s12*s12/max(VSMALL,xk2)
      endif
      return
      end
C
C=END SOURCE
