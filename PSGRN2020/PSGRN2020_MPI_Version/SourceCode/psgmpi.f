      subroutine psgmpi_init(myrank,nrank,ierr)
      implicit none
      include 'mpif.h'
c
      integer*4 myrank,nrank,ierr
c
      call MPI_INIT(ierr)
      if(ierr.ne.MPI_SUCCESS)then
        myrank=0
        nrank=1
        return
      endif
c
      call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
      if(ierr.ne.MPI_SUCCESS)return
      call MPI_COMM_SIZE(MPI_COMM_WORLD,nrank,ierr)
      return
      end
c
      subroutine psgmpi_finalize(ierr)
      implicit none
      include 'mpif.h'
c
      logical mpiflag
      integer*4 ierr
c
      call MPI_FINALIZED(mpiflag,ierr)
      if(ierr.ne.MPI_SUCCESS)return
      if(.not.mpiflag)then
        call MPI_FINALIZE(ierr)
      else
        ierr=0
      endif
      return
      end
c
      subroutine psgmpi_barrier(ierr)
      implicit none
      include 'mpif.h'
c
      integer*4 ierr
c
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      return
      end
c
      subroutine psgmpi_bcast_string(buffer,root,ierr)
      implicit none
      include 'mpif.h'
c
      character*(*) buffer
      integer*4 root,ierr
c
      call MPI_BCAST(buffer,len(buffer),MPI_CHARACTER,root,
     &               MPI_COMM_WORLD,ierr)
      return
      end
c
      subroutine psgmpi_bcast_i4(value,nvalue,root,ierr)
      implicit none
      include 'mpif.h'
c
      integer*4 nvalue,root,ierr
      integer*4 value(*)
c
      call MPI_BCAST(value,nvalue,MPI_INTEGER,root,
     &               MPI_COMM_WORLD,ierr)
      return
      end
c
      subroutine psgmpi_reduce_sum_i4(sendval,recvval,root,ierr)
      implicit none
      include 'mpif.h'
c
      integer*4 sendval,recvval,root,ierr
c
      call MPI_REDUCE(sendval,recvval,1,MPI_INTEGER,MPI_SUM,
     &                root,MPI_COMM_WORLD,ierr)
      return
      end
c
      subroutine psgmpi_send_task(dest,tag,izs,isp,nr1,nr2,
     &                            taskid,task_status,ierr)
      implicit none
      include 'mpif.h'
c
      integer*4 dest,tag,izs,isp,nr1,nr2,task_status,ierr
      integer*8 taskid
      integer*8 taskbuf(6)
c
      taskbuf(1)=int(izs,kind=8)
      taskbuf(2)=int(isp,kind=8)
      taskbuf(3)=int(nr1,kind=8)
      taskbuf(4)=int(nr2,kind=8)
      taskbuf(5)=taskid
      taskbuf(6)=int(task_status,kind=8)
      call MPI_SEND(taskbuf,6,MPI_INTEGER8,dest,tag,
     &              MPI_COMM_WORLD,ierr)
      return
      end
c
      subroutine psgmpi_recv_task(source,tag,izs,isp,nr1,nr2,
     &                            taskid,task_status,ierr)
      implicit none
      include 'mpif.h'
c
      integer*4 source,tag,izs,isp,nr1,nr2,task_status,ierr
      integer*4 mpistatus(MPI_STATUS_SIZE)
      integer*8 taskid
      integer*8 taskbuf(6)
c
      call MPI_RECV(taskbuf,6,MPI_INTEGER8,source,tag,
     &              MPI_COMM_WORLD,mpistatus,ierr)
      if(ierr.ne.MPI_SUCCESS)return
      izs=int(taskbuf(1),kind=4)
      isp=int(taskbuf(2),kind=4)
      nr1=int(taskbuf(3),kind=4)
      nr2=int(taskbuf(4),kind=4)
      taskid=taskbuf(5)
      task_status=int(taskbuf(6),kind=4)
      return
      end
c
      subroutine psgmpi_send_result(dest,tag,izs,isp,nr1,nr2,
     &                              resultid,result_status,ierr)
      implicit none
      include 'mpif.h'
c
      integer*4 dest,tag,izs,isp,nr1,nr2,result_status,ierr
      integer*8 resultid
      integer*8 resultbuf(6)
c
      resultbuf(1)=int(izs,kind=8)
      resultbuf(2)=int(isp,kind=8)
      resultbuf(3)=int(nr1,kind=8)
      resultbuf(4)=int(nr2,kind=8)
      resultbuf(5)=resultid
      resultbuf(6)=int(result_status,kind=8)
      call MPI_SEND(resultbuf,6,MPI_INTEGER8,dest,tag,
     &              MPI_COMM_WORLD,ierr)
      return
      end
c
      subroutine psgmpi_recv_result(source,tag,izs,isp,nr1,nr2,
     &                              resultid,result_status,ierr)
      implicit none
      include 'mpif.h'
c
      integer*4 source,tag,izs,isp,nr1,nr2,result_status,ierr
      integer*4 mpistatus(MPI_STATUS_SIZE)
      integer*8 resultid
      integer*8 resultbuf(6)
c
      call MPI_RECV(resultbuf,6,MPI_INTEGER8,source,tag,
     &              MPI_COMM_WORLD,mpistatus,ierr)
      if(ierr.ne.MPI_SUCCESS)return
      izs=int(resultbuf(1),kind=4)
      isp=int(resultbuf(2),kind=4)
      nr1=int(resultbuf(3),kind=4)
      nr2=int(resultbuf(4),kind=4)
      resultid=resultbuf(5)
      result_status=int(resultbuf(6),kind=4)
      return
      end
c
      subroutine psgmpi_recv_result_any(source,tag,izs,isp,nr1,nr2,
     &                                  resultid,result_status,ierr)
      implicit none
      include 'mpif.h'
c
      integer*4 source,tag,izs,isp,nr1,nr2,result_status,ierr
      integer*4 mpistatus(MPI_STATUS_SIZE)
      integer*8 resultid
      integer*8 resultbuf(6)
c
      call MPI_RECV(resultbuf,6,MPI_INTEGER8,MPI_ANY_SOURCE,tag,
     &              MPI_COMM_WORLD,mpistatus,ierr)
      if(ierr.ne.MPI_SUCCESS)return
      source=mpistatus(MPI_SOURCE)
      izs=int(resultbuf(1),kind=4)
      isp=int(resultbuf(2),kind=4)
      nr1=int(resultbuf(3),kind=4)
      nr2=int(resultbuf(4),kind=4)
      resultid=resultbuf(5)
      result_status=int(resultbuf(6),kind=4)
      return
      end
