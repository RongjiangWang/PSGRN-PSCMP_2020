      subroutine swavelet(tau,df,fi,nf,wvf)
      implicit none
      integer*4 nf
      real*8 tau,df,fi
      complex*16 wvf(nf)
c
      integer*4 l,n
      real*8 f,omi,dt0
      complex*16 alfa,beta,gamma,eta,cx
c
      real*8 pi,pi2,eps
      complex*16 c1,ci
      data pi,pi2,eps/3.14159265358979d0,6.28318530717959d0,1.0d-06/
      data c1,ci/(1.d0,0.d0),(0.d0,1.d0)/
c
c     for wavelet: normalized square half-sinus
c
      do l=1,nf
        f=df*dble(l-1)
        cx=dcmplx(f*tau,fi*tau)
        if(cdabs(cx).eq.0.d0)then
          wvf(l)=(1.d0,0.d0)
        else if(cdabs(cx-c1).le.eps)then
          wvf(l)=-c1/cx/(c1+cx)
        else
          wvf(l)=ci/(dcmplx(pi2,0.d0)*cx*(c1+cx)*(c1-cx))
     &          *(cdexp(dcmplx(0.d0,-pi2)*cx)-c1)
        endif
      enddo
      return
      end
