      program psgmain
      use psgalloc
      implicit none
c
c     work space
c
      integer*4 ierr,mpiierr,myrank,nrank,nwarnsum
      integer*4 runtime
      integer*4 time
c
      call psgmpi_init(myrank,nrank,mpiierr)
      if(mpiierr.ne.0)stop ' Error in psgmain: MPI init failed!'
c
      inputfile=' '
      nwarn=0
      if(myrank.eq.0)then
        print *,'######################################################'
        print *,'#                                                    #'
        print *,'#                  Welcome to                        #'
        print *,'#                                                    #'
        print *,'#     PPPP     SSSS    GGGG   RRRR    N   N          #'
        print *,'#     P   P   S       G       R   R   NN  N          #'
        print *,'#     PPPP     SSS    G GGG   RRRR    N N N          #'
        print *,'#     P           S   G   G   R R     N  NN          #'
        print *,'#     P       SSSS     GGGG   R  R    N   N          #'
        print *,'#                                                    #'
        print *,'#                  Version 2020                      #'
        print *,'# (update of version 2008b -> use of dynamic memory) #'
        print *,'#                                                    #'
        print *,'#                      by                            #'
        print *,'#                 Rongjiang Wang                     #'
        print *,'#              (wang@gfz-potsdam.de)                 #'
        print *,'#                                                    #'
        print *,'#             Helmholtz Centre Potsdam               #'
        print *,'#    GFZ German Research Centre for Geosciences      #'
        print *,'#                                                    #'
        print *,'#                parallelized by                     #'
        print *,'#            Yong Zheng & Xuhao Zou                  #'
        print *,'#   China University of Geosciences, Wuhan           #'
        print *,'#             Last modified: Mar 2026               #'
        print *,'######################################################'
        print *,'                                                      '
        if(nrank.gt.1)then
          write(*,'(a,i6,a)')' MPI runtime enabled with ',nrank,
     &                        ' ranks.'
        endif
        write(*,'(a,$)')' Please type the file name of input data: '
        read(*,'(a)')inputfile
      endif
      call psgmpi_bcast_string(inputfile,0,mpiierr)
      if(mpiierr.ne.0)stop ' Error in psgmain: input broadcast failed!'
c
      mpi_rank=myrank
      mpi_size=nrank
c
      if(myrank.eq.0)runtime=time()
c
c     read input file
c
      call psggetinp(ierr)
c
c     construction of sublayers
c
      call psgsublay(ierr)
c
c     main computation procedure
c
      call psgprocess(ierr)
c
      call psgmpi_reduce_sum_i4(nwarn,nwarnsum,0,mpiierr)
      if(mpiierr.ne.0)stop ' Error in psgmain: MPI reduce failed!'
      if(myrank.eq.0)then
        runtime=time()-runtime
        write(*,'(a)')'################################################'
        write(*,'(a)')'#                                              #'
        write(*,'(a)')'#        End of computations with PSGRN        #'
        write(*,'(a)')'#                                              #'
        if(nwarnsum.eq.0)then
          write(*,'(a,i10,a)')'#        Run time: ',runtime,
     &                                             ' sec              #'
        else
          write(*,'(a,i10,a)')'#        Run time: ',runtime,
     &                                             ' sec              #'
          write(*,'(a,i10,a)')'#        Warnings: ',nwarnsum,
     &                     '                  #'
        endif
        write(*,'(a)')'################################################'
      endif
c
      call psgmpi_barrier(mpiierr)
      if(mpiierr.ne.0)stop ' Error in psgmain: MPI barrier failed!'
      call psgmpi_finalize(mpiierr)
      if(mpiierr.ne.0)stop ' Error in psgmain: MPI finalize failed!'
c
      stop
      end
