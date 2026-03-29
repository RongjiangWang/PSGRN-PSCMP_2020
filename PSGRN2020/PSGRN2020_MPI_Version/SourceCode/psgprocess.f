      subroutine psgprocess(ierr)
      use psgalloc
      implicit none
c
      integer*4 ierr
c
      integer*4 nlr,nls,leninp,izs,itask
      integer*4 task_tag,result_tag
      integer*4 task_ready,task_done,task_stop
      integer*4 nexttask,completed,worker,source
      integer*4 rizs,risp,rnr1,rnr2,rstatus
      integer*8 rid
c
      task_tag=101
      result_tag=102
      task_ready=1
      task_done=3
      task_stop=0
c
      call psgallocwork(ierr)
      zs=0.d0
      call psglayer(ierr)
      nlr=nno(lzrec)
      leninp=index(inputfile,' ')-1
      call psginitlabels
      call psgbsj(ierr)
      call psgbuildtasklist(ierr)
c
      if(mpi_size.le.1)then
        call psgopenoutput(ierr,leninp,nlr)
        do izs=1,nzs
          call psgprepareizs(izs,nls,ierr,.true.)
          do itask=1,ntasks
            if(task_izs(itask).eq.izs)then
              call psgspec(task_isp(itask),task_nr1(itask),
     &                     task_nr2(itask))
            endif
          enddo
          call psgwriteblock(unit,izs,nls,ierr)
        enddo
        call psgcloseoutput
        return
      endif
c
      if(mpi_rank.eq.0)then
        call psgopenoutput(ierr,leninp,nlr)
        nexttask=1
        completed=0
        do worker=1,mpi_size-1
          if(nexttask.le.ntasks)then
            call psgmpi_send_task(worker,task_tag,
     &          task_izs(nexttask),task_isp(nexttask),
     &          task_nr1(nexttask),task_nr2(nexttask),
     &          task_id(nexttask),task_ready,ierr)
            nexttask=nexttask+1
          else
            call psgmpi_send_task(worker,task_tag,0,0,0,0,0_8,
     &          task_stop,ierr)
          endif
        enddo
        do while(completed.lt.ntasks)
          call psgmpi_recv_result_any(source,result_tag,rizs,risp,rnr1,
     &          rnr2,rid,rstatus,ierr)
          completed=completed+1
          if(nexttask.le.ntasks)then
            call psgmpi_send_task(source,task_tag,
     &          task_izs(nexttask),task_isp(nexttask),
     &          task_nr1(nexttask),task_nr2(nexttask),
     &          task_id(nexttask),task_ready,ierr)
            nexttask=nexttask+1
          else
            call psgmpi_send_task(source,task_tag,0,0,0,0,0_8,
     &          task_stop,ierr)
          endif
        enddo
        do izs=1,nzs
          call psgprepareizs(izs,nls,ierr,.false.)
          do itask=1,ntasks
            if(task_izs(itask).eq.izs)then
              call psgreadtaskfile(task_izs(itask),task_isp(itask),
     &                             task_nr1(itask),task_nr2(itask),ierr)
            endif
          enddo
          call psgwriteblock(unit,izs,nls,ierr)
        enddo
        call psgcloseoutput
      else
100     call psgmpi_recv_task(0,task_tag,rizs,risp,rnr1,rnr2,rid,
     &                        rstatus,ierr)
        if(rstatus.eq.task_stop)return
        write(*,'(/,a,i4,a,i4,a,i4,a)')' Rank ',mpi_rank,
     &      ' task izs=',rizs,' isp=',risp,'.'
        call psgprepareizs(rizs,nls,ierr,.false.)
        call psgspec(risp,rnr1,rnr2)
        call psgwritetaskfile(rizs,risp,rnr1,rnr2,ierr)
        call psgmpi_send_result(0,result_tag,rizs,risp,rnr1,rnr2,rid,
     &                          task_done,ierr)
        goto 100
      endif
      return
      end
c
      subroutine psgallocwork(ierr)
      use psgalloc
      implicit none
c
      integer*4 ierr,ntmax
c
      if(taumin.le.0.d0)then
        nfmax=nfmin
      else
        nfmax=nfmin
150     nfmax=2*nfmax
        if(dble(2*nfmax-1)*0.1d0*taumin.lt.twindow.and.
     &     nfmax.lt.max0(nfmax2min*nfmin,nt/2))goto 150
      endif
      allocate(tgrn(2*nfmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in psgprocess: tgrn not allocated!'
      allocate(fgrn(2*nfmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in psgprocess: fgrn not allocated!'
      allocate(wvf(nfmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in psgprocess: wvf not allocated!'
      allocate(dswap(4*nfmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in psgprocess: dswap not allocated!'
      allocate(rs(nr),stat=ierr)
      if(ierr.ne.0)stop ' Error in psgprocess: rs not allocated!'
      allocate(geow(nr),stat=ierr)
      if(ierr.ne.0)stop ' Error in psgprocess: geow not allocated!'
      allocate(obs(nr,16,4),stat=ierr)
      if(ierr.ne.0)stop ' Error in psgprocess: obs not allocated!'
      allocate(obs0(nr,16,4),stat=ierr)
      if(ierr.ne.0)stop ' Error in psgprocess: obs0 not allocated!'
      allocate(idec(nr),stat=ierr)
      if(ierr.ne.0)stop ' Error in psgprocess: idec not allocated!'
      allocate(nout(nr),stat=ierr)
      if(ierr.ne.0)stop ' Error in psgprocess: nout not allocated!'
      ntmax=2*max0(nfmax,(1+nt)/2)
      allocate(dobs(ntmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in psgprocess: obs not allocated!'
      allocate(du(-1:ntmax/2,nr,14,4),stat=ierr)
      if(ierr.ne.0)stop ' Error in psgprocess: du not allocated!'
      lp=n0+2
      allocate(hp(lp),stat=ierr)
      if(ierr.ne.0)stop ' Error in psgprocess: hp not allocated!'
      allocate(nno(lp),stat=ierr)
      if(ierr.ne.0)stop ' Error in psgprocess: nno not allocated!'
      allocate(zp(lp),stat=ierr)
      if(ierr.ne.0)stop ' Error in psgprocess: zp not allocated!'
      allocate(hkup(2,2,lp),stat=ierr)
      if(ierr.ne.0)stop ' Error in psgprocess: hkup not allocated!'
      allocate(hklw(2,2,lp),stat=ierr)
      if(ierr.ne.0)stop ' Error in psgprocess: hklw not allocated!'
      allocate(maup(6,6,lp),stat=ierr)
      if(ierr.ne.0)stop ' Error in psgprocess: maup not allocated!'
      allocate(maiup(6,6,lp),stat=ierr)
      if(ierr.ne.0)stop ' Error in psgprocess: maiup not allocated!'
      allocate(malw(6,6,lp),stat=ierr)
      if(ierr.ne.0)stop ' Error in psgprocess: malw not allocated!'
      allocate(mailw(6,6,lp),stat=ierr)
      if(ierr.ne.0)stop ' Error in psgprocess: mailw not allocated!'
      return
      end
c
      subroutine psginitlabels
      use psgalloc
      implicit none
c
      stype(1)='explosion (M11=M22=M33=1*kappa)'
      stype(2)='strike-slip (M12=M21=1*mue)'
      stype(3)='dip-slip (M13=M31=1*mue)'
      stype(4)='clvd (M33=1*mue, M11=M22=-M33/2)'
      comptxt(1)='Uz (vertical displacement)'
      comptxt(2)='Ur (radial displacement)'
      comptxt(3)='Ut (tangential displacement)'
      comptxt(4)='Szz (linear stress)'
      comptxt(5)='Srr (linear stress)'
      comptxt(6)='Stt (linear stress)'
      comptxt(7)='Szr (shear stress)'
      comptxt(8)='Srt (shear stress)'
      comptxt(9)='Stz (shear stress)'
      comptxt(10)='Tr (tilt -dUr/dz)'
      comptxt(11)='Tt (tilt -dUt/dz)'
      comptxt(12)='Rot (rotation ar. z-axis)'
      comptxt(13)='Gd (geoid changes)'
      comptxt(14)='Gr (gravity changes)'
      return
      end
c
      subroutine psgbuildtasklist(ierr)
      use psgalloc
      implicit none
c
      integer*4 ierr,izs,isp,ir,nprf,maxtask,nr1,nr2
      real*8 zrs2w,dratio,swap
c
      ierr=0
      maxtask=max0(1,nzs*nr)
      if(.not.allocated(task_id))then
        allocate(task_id(maxtask),stat=ierr)
        if(ierr.ne.0)stop ' Error in psgprocess: task_id not allocated!'
        allocate(task_izs(maxtask),stat=ierr)
        if(ierr.ne.0)stop ' psgprocess: task_izs alloc error'
        allocate(task_isp(maxtask),stat=ierr)
        if(ierr.ne.0)stop ' psgprocess: task_isp alloc error'
        allocate(task_nr1(maxtask),stat=ierr)
        if(ierr.ne.0)stop ' psgprocess: task_nr1 alloc error'
        allocate(task_nr2(maxtask),stat=ierr)
        if(ierr.ne.0)stop ' psgprocess: task_nr2 alloc error'
        allocate(task_status(maxtask),stat=ierr)
        if(ierr.ne.0)stop ' psgprocess: task_status alloc error'
      endif
      ntasks=0
      do izs=1,nzs
        zrs2w=(zrec-(zs1+dble(izs-1)*dzs))**2
        do ir=1,nr
          rs(ir)=rsmin+0.01d0*dsqrt(zrs2w+r(ir)**2)
        enddo
        swap=dsqrt(zrs2w+(rs(nr)+r(nr))**2)
     &      /dsqrt(zrs2w+(rs(1)+r(1))**2)
        nprf=1+idnint(dlog(swap)/dlog(5.d0))
        if(nprf.gt.1)then
          dratio=swap**(1.d0/dble(nprf-1))
        else
          dratio=2.5d0
        endif
        isp=0
        nr2=0
200     isp=isp+1
        nr1=nr2+1
        nr2=nr1
        do ir=nr1+1,nr
          if(r(ir).le.dratio*dsqrt(zrs2w+(rs(nr1)+r(nr1))**2))nr2=ir
        enddo
        ntasks=ntasks+1
        task_id(ntasks)=int(ntasks,kind=8)
        task_izs(ntasks)=izs
        task_isp(ntasks)=isp
        task_nr1(ntasks)=nr1
        task_nr2(ntasks)=nr2
        task_status(ntasks)=1
        if(nr2.lt.nr)goto 200
      enddo
      return
      end
c
      subroutine psgprepareizs(izs,nls,ierr,verbose)
      use psgalloc
      implicit none
c
      logical verbose
      integer*4 izs,nls,ierr
      integer*4 i,l,ir,it,istp
      real*8 zrs2w
c
      zs=zs1+dble(izs-1)*dzs
      if(verbose)then
        write(*,'(/,a,i4,a,E13.4,a)')' Rank ',mpi_rank,
     &               ' processing source at depth:',zs,' m.'
      endif
      call psglayer(ierr)
      do l=1,lp
        zp(l)=0.d0
        do i=1,l-1
          if(nno(i).eq.nno(l))zp(l)=zp(l)+hp(i)
        enddo
      enddo
      nls=nno(ls)
      do istp=1,4
        do i=1,14
          do ir=1,nr
            do it=-1,ubound(du,1)
              du(it,ir,i,istp)=(0.d0,0.d0)
            enddo
          enddo
        enddo
      enddo
      zrs2w=(zrec-zs)**2
      zrs2=zrs2w
      do ir=1,nr
        rs(ir)=rsmin+0.01d0*dsqrt(zrs2w+r(ir)**2)
        geow(ir)=zrs2w+(rs(ir)+r(ir))**2
      enddo
      return
      end
