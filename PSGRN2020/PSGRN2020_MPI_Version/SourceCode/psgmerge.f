      subroutine psgopenoutput(ierr,leninp,nlr)
      use psgalloc
      implicit none
c
      integer*4 ierr,leninp,nlr
      integer*4 i,istp,iunit
c
      ierr=0
      iunit=10
      do istp=1,4
        do i=1,14
          unit(i,istp)=0
          if(select(i,istp))then
            iunit=iunit+1
            unit(i,istp)=iunit
            open(unit(i,istp),file=green(i,istp),status='unknown')
            write(unit(i,istp),'(a)')'################################'
            write(unit(i,istp),'(a)')'# The input file used: '
     &                        //inputfile(1:leninp)
            write(unit(i,istp),'(a)')'################################'
            write(unit(i,istp),'(a)')'# Greens function component: '
     &                        //comptxt(i)
            write(unit(i,istp),'(a)')'#(Okada solutions subtracted)'
            write(unit(i,istp),'(a)')'# Source type: '//stype(istp)
            write(unit(i,istp),'(a)')'# Observation distance sampling:'
            write(unit(i,istp),'(a)')'#    nr        r1[m]        r2[m]'
     &                             //'  samp_ratio'
            write(unit(i,istp),'(i7,2E14.6,f10.4)')nr,r1,r2,sampratio
            write(unit(i,istp),'(a)')'# Uniform obs. site parameters:'
            write(unit(i,istp),'(a)')'#    depth[m]       la[Pa]       '
     &       //'mu[Pa]  rho[kg/m^3]    etk[Pa*s]    etm[Pa*s]     alpha'
            write(unit(i,istp),'(7E13.6)')zrec,la(nlr),
     &           mu(nlr),rho(nlr),etk(nlr),etm(nlr),alf(nlr)
            write(unit(i,istp),'(a)')'# Source depth sampling:'
            write(unit(i,istp),'(a)')'#   nzs       zs1[m]       zs2[m]'
            write(unit(i,istp),'(i7,2d14.6)')nzs,zs1,zs2
            write(unit(i,istp),'(a)')'# Time sampling:'
            write(unit(i,istp),'(a)')'#    nt        t-window[s]'
            write(unit(i,istp),'(i7,E24.16)')nt,twindow
            write(unit(i,istp),'(a)')'# Data in each source depth block'
            write(unit(i,istp),'(a)')'# ==============================='
            write(unit(i,istp),'(a)')'# Line 1: source layer parameters'
            write(unit(i,istp),'(a)')'#  s_depth, la, mu, rho,'
     &                              //' etk, etm, alpha'
            write(unit(i,istp),'(a)')'# Line 2: coseismic responses '
     &                             //'(f(ir,it=1),ir=1,nr)'
            write(unit(i,istp),'(a)')'# Line 3: (idec(ir),ir=1,nr)'
     &                //'(decimal exponents for postseismic responses)'
            write(unit(i,istp),'(a)')'# Line 4: (f(ir,it=2),ir=1,nr)'
            write(unit(i,istp),'(a)')'# Line 5: (f(ir,it=3),ir=1,nr)'
            write(unit(i,istp),'(a)')'#  ...'
          endif
        enddo
      enddo
      return
      end
c
      subroutine psgwriteblock(outunit,izs,nls,ierr)
      use psgalloc
      implicit none
c
      integer*4 outunit(14,4),izs,nls,ierr
      integer*4 i,ir,it,istp
      real*8 am
c
      ierr=0
      do istp=1,4
        do i=1,14
          if(.not.select(i,istp))goto 400
          write(outunit(i,istp),'(a)')'##############################'
     &                              //'###'
          write(outunit(i,istp),'(a,i2,a)')'# the ',izs,
     &                              '. source depth:'
          write(outunit(i,istp),'(a)')'##############################'
     &                              //'###'
          write(outunit(i,istp),'(7E13.6)')zs,la(nls),
     &      mu(nls),rho(nls),etk(nls),etm(nls),alf(nls)
          do ir=1,nr-1
            write(outunit(i,istp),'(E14.6,$)')dreal(du(1,ir,i,istp))
          enddo
          write(outunit(i,istp),'(E14.6)')dreal(du(1,nr,i,istp))
          do ir=1,nr
            du(1,ir,i,istp)=dcmplx(0.d0,dimag(du(1,ir,i,istp)))
            am=0.d0
            do it=1,(nt+1)/2
              am=dmax1(am,dabs(dreal(du(it,ir,i,istp))),
     &                    dabs(dimag(du(it,ir,i,istp))))
            enddo
            if(am.le.0.d0)then
              idec(ir)=0
            else
              idec(ir)=idint(dlog10(am))-ndigit
              do it=1,(nt+1)/2
                du(it,ir,i,istp)=du(it,ir,i,istp)
     &                *dcmplx(10.d0**(-idec(ir)),0.d0)
              enddo
            endif
          enddo
          call outint(outunit(i,istp),idec,nr)
          do it=1,(nt+1)/2
            if(it.gt.1)then
              do ir=1,nr
                nout(ir)=idnint(dreal(du(it,ir,i,istp)))
              enddo
              call outint(outunit(i,istp),nout,nr)
            endif
            if(it*2.le.nt)then
              do ir=1,nr
                nout(ir)=idnint(dimag(du(it,ir,i,istp)))
              enddo
              call outint(outunit(i,istp),nout,nr)
            endif
          enddo
400       continue
        enddo
      enddo
      return
      end
c
      subroutine psgtaskfile(izs,isp,tmpfile)
      use psgalloc
      implicit none
c
      integer*4 izs,isp,lend
      character*(*) tmpfile
      character*5 ztxt,ptxt
c
      lend=index(outdir,' ')-1
      if(lend.le.0)lend=len(outdir)
      write(ztxt,'(i5.5)')izs
      write(ptxt,'(i5.5)')isp
      tmpfile=' '
      tmpfile=outdir(1:lend)//'mpi_izs'//ztxt//'_isp'//ptxt//'.bin'
      return
      end
c
      subroutine psgwritetaskfile(izs,isp,nr1,nr2,ierr)
      use psgalloc
      implicit none
c
      integer*4 izs,isp,nr1,nr2,ierr,ntc,iunit
      character*255 tmpfile
c
      ierr=0
      ntc=(nt+1)/2
      iunit=80+mpi_rank
      call psgtaskfile(izs,isp,tmpfile)
      open(iunit,file=tmpfile,status='unknown',form='unformatted')
      write(iunit)du(1:ntc,nr1:nr2,1:14,1:4)
      close(iunit)
      return
      end
c
      subroutine psgreadtaskfile(izs,isp,nr1,nr2,ierr)
      use psgalloc
      implicit none
c
      integer*4 izs,isp,nr1,nr2,ierr,ntc,iunit
      character*255 tmpfile
      logical ex
c
      ierr=0
      ntc=(nt+1)/2
      iunit=95
      call psgtaskfile(izs,isp,tmpfile)
      inquire(file=tmpfile,exist=ex)
      if(.not.ex)then
        ierr=1
        return
      endif
      open(iunit,file=tmpfile,status='old',form='unformatted')
      read(iunit)du(1:ntc,nr1:nr2,1:14,1:4)
      close(iunit,status='delete')
      return
      end
c
      subroutine psgcloseoutput
      use psgalloc
      implicit none
c
      integer*4 i,istp
c
      do istp=1,4
        do i=1,14
          if(select(i,istp).and.unit(i,istp).gt.0)close(unit(i,istp))
        enddo
      enddo
      return
      end
