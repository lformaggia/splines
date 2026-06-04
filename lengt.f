C=DECK lengt
C=NAME lengt
C=RCSfile $Source: /u/flosys3d/software/cvsroot/splines/lengt.f,v $
C=VERSION $Revision: 1.1.1.1 $ $Date: 1996/03/05 17:16:38 $
C=TYPE subroutine
C=KEYWORDS Parametric_cubic_spline
C=PURPOSE general
C=BLOCK ABSTRACT
C
C It evaluates the length of the portion of a cubic spline arc between
C the parameter coordinates u1 and u2.It requires as input the value of
c the 5 coefficients of the polinomial expansion for (x,u)**2
c
c |x,u|**2 = a1 + a2*u + a3*u**2+a4*u**3+a5*u**4
c
c those coefficent may be computed by using the routine coeff
C
c It computes  int[u1,u2](||x,u||du)     (u1<=u2 )
c using Romberg's adaptive quadrature
C
C=END ABSTRACT
c=BLOCK USAGE
C
C  call lengt(in,a1,a2,a3,a4,a5,u1,u2,s,eps)
C
C  INPUT                    DESCRIPTION
c  a1...a5                  coefficents of the polinomial P = ||x,u||**2
c              P = a5*u**4+a4*u**3 + ... + a1     (see routine coeff)
C  u1,u2                    integration limits
C  eps                      tolerance for accuracy of the estimated length
C  in                       Obsolete parameter. The integral is always
c                           computed to fractional accuracy eps
c
c  OUTPUT
c  s                        computed length
c
c  DIAGNOSTIC
C
c  If the max. number of iterations are exceeded in the auxiliary routine
c  QROMBS, an error message is printed on the standard output
C
C=END USAGE
C
C=BLOCK SOURCE
C
      subroutine lengt(in,a1,a2,a3,a4,a5,u1,u2,s,eps)
c
c *** this sub. computes the length of a segment of a cubic.
c     using a romberg adaptive integration
c
      if(u1-u2.eq.0.)then
         s=0.
      else
         call qrombs(a1,a2,a3,a4,a5,u1,u2,s,eps)
      endif
      return
      end
C=END SOURCE
C
C
C=NAME lengt_aux
C=TYPE source_library
C=BLOCK ABSTRACT
C
C Auxiliary routines which perform Romberg integration.
c Modification of  some routines taken from Numerical Recipies (Ed.1)
c The name of the corresponding Numerical Recipes routines 
c is obtained by dropping  the last 's', i.e.  
c Routine name: polints -> Nume. Rec. name: polint  
c
c Refer to: Numerical Recipies, The Art of Scientific       
c           Computing (1st Ed.), by  W.H. Press et al.,
c           Cambridge Press, 1989.
c
c for details on the routines. 
C
C=END ABSTRACT
C=BLOCK SOURCE
      subroutine polints(xa,ya,n,x,y,dy)
      parameter (nmax=10) 
      dimension xa(n),ya(n),c(nmax),d(nmax)
      ns=1
      dif=abs(x-xa(1))
      do 11 i=1,n 
        dift=abs(x-xa(i))
        if (dift.lt.dif) then
          ns=i
          dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
11    continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if(den.eq.0.)then
             print*,'POLINTS: E001:  Error'
             return
          endif
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
12      continue
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
13    continue
      return
      end

      subroutine trapzds(a1,a2,a3,a4,a5,a,b,s,n)
      save it
      func(u1) = sqrt(a1+u1*(a2+u1*(a3+u1*(a4+u1*a5))))
c
      if (n.eq.1) then
        s=0.5*(b-a)*(func(a)+func(b))
        it=1
      else
        tnm=it
        del=(b-a)/tnm
        x=a+0.5*del
        sum=0.
        do 11 j=1,it
          sum=sum+func(x)
          x=x+del
11      continue
        s=0.5*(s+(b-a)*sum/tnm)
        it=2*it
      endif
      return
      end

      subroutine qrombs(a1,a2,a3,a4,a5,a,b,ss,eps)
      parameter (jmax=20, jmaxp=jmax+1, k=2, km=k-1)
      parameter (FLTMIN=1.e-30)
      dimension s(jmaxp),h(jmaxp)
      h(1)=1.
      do 11 j=1,jmax
        call trapzds(a1,a2,a3,a4,a5,a,b,s(j),j)
        if (j.ge.k) then
          call polints(h(j-km),s(j-km),k,0.,ss,dss)
          if (abs(dss).le.eps*abs(ss).or.
     $         (j.gt.2.and.abs(dss).lt.FLTMIN)) return
        endif
        s(j+1)=s(j)
        h(j+1)=0.25*h(j)
11    continue
      print*,  'QROMBS: too many steps.'
      end
C=END SOURCE
C=END DECK










