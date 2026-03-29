      module psgtask
c===================================================================
c     MPI task/result metadata for (izs, isp) orchestration
c===================================================================
      implicit none
c
c     status flags
c     status values are intentionally simple integers so they can be
c     copied through MPI messages without extra packing logic.
c
      integer*4 psgtask_status_empty, psgtask_status_ready
      integer*4 psgtask_status_running, psgtask_status_done
      integer*4 psgtask_status_failed
      parameter(psgtask_status_empty=0,
     &          psgtask_status_ready=1,
     &          psgtask_status_running=2,
     &          psgtask_status_done=3,
     &          psgtask_status_failed=-1)
c
c     task metadata
c
      type psgtask_meta
        sequence
        integer*4 izs
        integer*4 isp
        integer*4 nr1
        integer*4 nr2
        integer*8 taskid
        integer*4 status
      end type psgtask_meta
c
c     result metadata
c
      type psgresult_meta
        sequence
        integer*4 izs
        integer*4 isp
        integer*4 nr1
        integer*4 nr2
        integer*8 resultid
        integer*4 status
      end type psgresult_meta
c
      contains
c===================================================================
c     helpers for task indexing and metadata initialization
c===================================================================
      integer*8 function psgtask_key(izs, isp)
      implicit none
      integer*4 izs, isp
c
c     pack (izs, isp) into a stable 64-bit identifier without relying
c     on a fixed span or 32-bit arithmetic.
c
      psgtask_key = ishft(int(izs,kind=8),32) + int(isp,kind=8)
      return
      end
c
      subroutine psgtask_init(task, izs, isp, nr1, nr2)
      implicit none
      type(psgtask_meta) task
      integer*4 izs, isp, nr1, nr2
      task%izs = izs
      task%isp = isp
      task%nr1 = nr1
      task%nr2 = nr2
      task%taskid = psgtask_key(izs, isp)
      task%status = psgtask_status_ready
      return
      end
c
      subroutine psgtask_reset(task)
      implicit none
      type(psgtask_meta) task
      task%izs = 0
      task%isp = 0
      task%nr1 = 0
      task%nr2 = 0
      task%taskid = 0
      task%status = psgtask_status_empty
      return
      end
c
      subroutine psgtask_set_status(task, status)
      implicit none
      type(psgtask_meta) task
      integer*4 status
      task%status = status
      return
      end
c
      logical function psgtask_is_complete(task)
      implicit none
      type(psgtask_meta) task
      psgtask_is_complete = (task%status.eq.psgtask_status_done)
      return
      end
c
      subroutine psgresult_init(result, izs, isp, nr1, nr2)
      implicit none
      type(psgresult_meta) result
      integer*4 izs, isp, nr1, nr2
      result%izs = izs
      result%isp = isp
      result%nr1 = nr1
      result%nr2 = nr2
      result%resultid = psgtask_key(izs, isp)
      result%status = psgtask_status_ready
      return
      end
c
      subroutine psgresult_from_task(result, task)
      implicit none
      type(psgresult_meta) result
      type(psgtask_meta) task
      result%izs = task%izs
      result%isp = task%isp
      result%nr1 = task%nr1
      result%nr2 = task%nr2
      result%resultid = task%taskid
      result%status = task%status
      return
      end
c
      subroutine psgresult_set_status(result, status)
      implicit none
      type(psgresult_meta) result
      integer*4 status
      result%status = status
      return
      end
c
      logical function psgresult_is_complete(result)
      implicit none
      type(psgresult_meta) result
      psgresult_is_complete = (result%status.eq.psgtask_status_done)
      return
      end
      end module
