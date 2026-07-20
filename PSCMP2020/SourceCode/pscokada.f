      subroutine pscokada(tmax)
      use pscalloc
      implicit none
      real*8 tmax
c
c     Last modified: Potsdam, Feb, 2019, by R. Wang
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c     from Okada's subroutine DC3D0:
c
      INTEGER IRET
      REAL*4 ALPHA,X,Y,Z,DEPTH,DIPS,
     &       UX,UY,UZ,UXX,UYX,UZX,UXY,UYY,UZY,UXZ,UYZ,UZZ
c
c     more from Okada's subroutine DC3D:
c
      REAL*4 AL1,AL2,AW1,AW2,DISL1,DISL2,DISL3
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c     LOCAL WORK SPACES
c     =================
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      integer*4 i,j,k,l,ieq,is,iptch,irec
      real*8 st,di,step_s,step_d,disn,dise
      real*8 csst,ssst,csra,ssra,csdi,ssdi
      real*8 cs2st,ss2st,eii,phi,azi,bazi,dis
      real*8 strain(6),sig(3,3),rot(3,3),swp(3,3)
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     PROCESSING:
c     ===========
c     coobs(1)=' Ux'
c     coobs(2)=' Uy'
c     coobs(3)=' Uz'
c     coobs(4)='Sxx'
c     coobs(5)='Syy'
c     coobs(6)='SZZ'
c     coobs(7)='Sxy'
c     coobs(8)='Syz'
c     coobs(9)='Szx'
c     coobs(10)=' Tx'
c     coobs(11)=' Ty'
c     coobs(12)='Rot'
c     coobs(13)=' Gd'
c     coobs(14)=' Gr'
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c     receiver and source independent variables
c
      ALPHA=sngl((larec+murec)/(larec+2.d0*murec))
      Z=-sngl(zrec)
c
      do is=1,ns
        if(tstart(is).gt.tmax)goto 100
        ieq=ieqno(is)
        st=strike(is)*DEG2RAD
        csst=dcos(st)
        ssst=dsin(st)
        cs2st=dcos(2.d0*st)
        ss2st=dsin(2.d0*st)
c
	  di=dip(is)*DEG2RAD
        csdi=dcos(di)
        ssdi=dsin(di)
c
        DEPTH=sngl(zref(is))
        DIPS=sngl(dip(is))
c
        step_s=length(is)/dble(nptch_s(is))
        step_d=width(is)/dble(nptch_d(is))
c
        do iptch=1,nptch_s(is)*nptch_d(is)
          DISL1=sngl( slip_s(is,iptch))
          DISL2=sngl(-slip_d(is,iptch))
          DISL3=sngl(opening(is,iptch))
c
c         for extended source
c
          AL1=sngl(ptch_s(is,iptch)-0.5d0*step_s)
          AL2=sngl(ptch_s(is,iptch)+0.5d0*step_s)
          AW1=sngl(-ptch_d(is,iptch)-0.5d0*step_d)
          AW2=sngl(-ptch_d(is,iptch)+0.5d0*step_d)
          do irec=1,nrec
c
c           transform from Aki's to Okada's system
c
            call disazi(REARTH,latref(is),lonref(is),
     &                  latrec(irec),lonrec(irec),disn,dise)
            X=sngl(disn*csst+dise*ssst)
            Y=sngl(disn*ssst-dise*csst)
c
            dis=dsqrt(disn**2+dise**2)
            if(dis.gt.ptch_s(is,iptch)+ptch_d(is,iptch))then
              azi=datan2(dise,disn)
              call disazi(REARTH,latrec(irec),lonrec(irec),
     &                    latref(is),lonref(is),disn,dise)
              bazi=datan2(dise,disn)
              phi=azi-(bazi-PI)
            else
              phi=0.d0
            endif
            rot(1,1)=dcos(phi)
            rot(1,2)=-dsin(phi)
            rot(1,3)=0.d0
            rot(2,1)=dsin(phi)
            rot(2,2)=dcos(phi)
            rot(2,3)=0.d0
            rot(3,1)=0.d0
            rot(3,2)=0.d0
            rot(3,3)=1.d0
c
            IRET=1
            call DC3D(ALPHA,X,Y,Z,DEPTH,DIPS,AL1,AL2,AW1,AW2,
     &            DISL1,DISL2,DISL3,UX,UY,UZ,
     &            UXX,UYX,UZX,UXY,UYY,UZY,UXZ,UYZ,UZZ,IRET)
c
c           transform from Okada's to Aki's system
c
            swp(1,1)=dble(UX)*csst+dble(UY)*ssst
            swp(2,1)=dble(UX)*ssst-dble(UY)*csst
            swp(3,1)=-dble(UZ)
c
            coobs(ieq,irec,1)=coobs(ieq,irec,1)
     &                       +rot(1,1)*swp(1,1)+rot(1,2)*swp(2,1)
            coobs(ieq,irec,2)=coobs(ieq,irec,2)
     &                       +rot(2,1)*swp(1,1)+rot(2,2)*swp(2,1)
            coobs(ieq,irec,3)=coobs(ieq,irec,3)+swp(3,1)
c
            swp(1,1)=dble(UXX)*csst*csst+dble(UYY)*ssst*ssst
     &               +0.5d0*dble(UXY+UYX)*ss2st
            swp(2,2)=dble(UXX)*ssst*ssst+dble(UYY)*csst*csst
     &               -0.5d0*dble(UXY+UYX)*ss2st
            swp(3,3)=dble(UZZ)
            swp(1,2)=0.5d0*dble(UXX-UYY)*ss2st
     &               -0.5d0*dble(UXY+UYX)*cs2st
            swp(2,1)=swp(1,2)
            swp(2,3)=-0.5d0*dble(UZX+UXZ)*ssst
     &               +0.5d0*dble(UYZ+UZY)*csst
            swp(3,2)=swp(2,3)
            swp(3,1)=-0.5d0*dble(UZX+UXZ)*csst
     &               -0.5d0*dble(UYZ+UZY)*ssst
            swp(1,3)=swp(3,1)
c
            do i=1,3
              do j=1,3
                sig(i,j)=0.d0
                do k=1,3
                  do l=1,3
                    sig(i,j)=sig(i,j)+rot(i,k)*swp(k,l)*rot(j,l)
                  enddo
                enddo
              enddo
            enddo
c
            strain(1)=sig(1,1)
            strain(2)=sig(2,2)
            strain(3)=sig(3,3)
            strain(4)=sig(1,2)
            strain(5)=sig(2,3)
            strain(6)=sig(3,1)
c
            eii=strain(1)+strain(2)+strain(3)
            coobs(ieq,irec,4)=coobs(ieq,irec,4)
     &                      +larec*eii+2.d0*murec*strain(1)
            coobs(ieq,irec,5)=coobs(ieq,irec,5)
     &                      +larec*eii+2.d0*murec*strain(2)
            coobs(ieq,irec,7)=coobs(ieq,irec,7)+2.d0*murec*strain(4)
            if(zrec.gt.0.d0)then
              coobs(ieq,irec,6)=coobs(ieq,irec,6)
     &                         +larec*eii+2.d0*murec*strain(3)
              coobs(ieq,irec,8)=coobs(ieq,irec,8)+2.d0*murec*strain(5)
              coobs(ieq,irec,9)=coobs(ieq,irec,9)+2.d0*murec*strain(6)
            endif
c
            swp(1,1)=-(dble(UXZ)*csst+dble(UYZ)*ssst)
            swp(2,1)=-(dble(UXZ)*ssst-dble(UYZ)*csst)
c
            coobs(ieq,irec,10)=coobs(ieq,irec,10)
     &                        +rot(1,1)*swp(1,1)+rot(1,2)*swp(2,1)
            coobs(ieq,irec,11)=coobs(ieq,irec,11)
     &                        +rot(2,1)*swp(1,1)+rot(2,2)*swp(2,1)
c
            coobs(ieq,irec,12)=coobs(ieq,irec,12)-0.5d0*dble(UYX-UXY)
            coobs(ieq,irec,14)=coobs(ieq,irec,14)
     &                        -(2.d0*G0/REARTH)*dble(UZ)
          enddo
        enddo
100     continue
      enddo
c
      return
      end
