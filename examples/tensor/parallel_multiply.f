! This file is deprecated. tensorpotentialmultiply_generator.py can generate a more
! flexible version of this routine
!


module IndexTricks
contains
subroutine MapRowToColPacked(row, col, fullSize, bands, packedRow, packedCol)
      implicit none
      integer, intent(in) :: row, col, fullSize, bands
      integer, intent(out) :: packedRow, packedCol

      packedRow = col
      packedCol = bands - col + row
end subroutine

function GetGlobalStartIndex(fullSize, procCount, procId) result(globalStartIndex)
      implicit none
      integer, intent(in) :: fullSize, procCount, procId
      integer :: rest, paddedSize, globalStartIndex

      paddedSize  = fullSize
      rest = mod(fullSize, procCount)
      if (rest .ne. 0) then
          paddedSize = paddedSize + procCount - rest
      endif
      globalStartIndex = (paddedSize / procCount) * procId

      return
end function     

function GetDistributedSize(fullSize, procCount, procId) result(distributedSize)
      implicit none
      integer, intent(in) :: fullSize, procCount, procId
      integer :: rest, distributedSize, paddedDistribSize

      rest = mod(fullSize, procCount)
      if (rest.eq.0) then
          distributedSize = fullSize / procCount;
      else
          paddedDistribSize = (fullSize + procCount - rest) / procCount 
          distributedSize = fullSize - paddedDistribSize * procId
          distributedSize = max(distributedSize, 0)
          distributedSize = min(distributedSize, paddedDistribSize)
      endif

      return
end function 

function GetPaddedLocalSize(fullSize, procCount) result(localSize)
      implicit none
      integer, intent(in) :: fullSize, procCount
      integer :: rest, paddedSize, localSize

      paddedSize  = fullSize
      rest = mod(fullSize, procCount)
      if (rest .ne. 0) then
          paddedSize = paddedSize + procCount - rest
      endif
      localSize = paddedSize / procCount

      return
end function 

function GetOwnerProcId(fullSize, procCount, globalIndex) result(procId)
      implicit none
      integer, intent(in) :: fullSize, procCount, globalIndex
      integer :: rest, paddedSize, procId, localSize

      include "parameters.f"
        
      paddedSize = fullSize
      rest = mod(fullSize, procCount)
      if (rest .ne. 0) then
          paddedSize = paddedSize + procCount - rest
      endif
      localSize = paddedSize / procCount
      !Fortran does not 
      procId = floor( real(globalIndex, kind=dbl) / real(localSize, kind=dbl) )

      return
end function      
     
end module

subroutine BandedMatrixMultiply(potential, potentialExtent0, scaling, source, sourceExtent0, dest, destExtent0, globalSize0, bands0, procCount0, procId0, communicator0 )
      ! See parallel_multiply.py for a more thorough explanation of how this
      ! works
      use IndexTricks
      implicit none
      include "parameters.f"
      include "mpif.h"
      
      complex (kind=dbl), dimension(0:potentialExtent0-1), intent(in) :: potential
      integer, intent(in) :: potentialExtent0
      real (kind=dbl) :: scaling
      
      complex (kind=dbl), dimension(0:sourceExtent0-1), intent(in) :: source
      integer, intent(in) :: sourceExtent0
      
      complex (kind=dbl), dimension(0:destExtent0-1), intent(inout) :: dest
      integer, intent(in) :: destExtent0
      
      integer, intent(in) :: globalSize0, bands0, procCount0, procId0, communicator0

      integer :: N0, i0, row0, col0
     
      !temporary arrays
      complex (kind=dbl), dimension(0:0) :: tempValue, recvTemp0, sendTemp0

      !Row/Col indices for calculating row vector product
      integer :: globalPackedRow0, globalPackedCol0, localPackedRow0, localPackedCol0, localStartCol0, localEndCol0
      integer :: globalCol0, globalRow0 

      !Indices to map between local and global index ranges
      integer :: globalStartIndex0, paddedLocalSize0, localSize0

      !MPI variables
      integer :: deltaProc0, sourceProc0, destProc0
      integer :: error0, tag0, sendSize0, recvRequest0, sendRequest0, waitRecieve0, waitSend0

      !Indices for the recieved data
      integer :: sourceRow0, sourceGlobalStartIndex0, sourceGlobalRow0

      N0 = potentialExtent0
      localSize0 = sourceExtent0

      globalStartIndex0 = GetGlobalStartIndex(globalSize0, procCount0, procId0)
      paddedLocalSize0 = GetPaddedLocalSize(globalSize0, procCount0)

      waitRecieve0 = 0
      waitSend0 = 0
      sourceRow0 = -1

      sendSize0 = 1
      tag0 = 0

      col0 = 0
      row0 = 0
      i0 = 0

      !Iterate over all rows of the matrix for the columns stored on proc
      do row0 = -bands0, paddedLocalSize0+bands0-1

          !----------------------------------------------------------------------------
          !                         Computation
          !----------------------------------------------------------------------------

          ! Perform multiplication of the local part of the vector on 
          ! the part of the array that this proc has
          tempValue(0) = 0
          globalRow0 = row0 + globalStartIndex0
          ! If the local row-index actually maps to a row in the global matrix
          ! (otherwise we're on a proc near the boundary)
          if (0.le.globalRow0.and.globalRow0.lt.globalSize0) then
              
              !Loop over all cols that are stored locally
              localStartCol0 = max(0, row0 - bands0)
              localEndCol0 = min(row0 + bands0 + 1, localSize0)
              do col0 = localStartCol0, localEndCol0-1
                  globalCol0 = col0 + globalStartIndex0
                  
                  !PackedRows and -Cols are for C-style arrays
                  !i.e. row indexes the max strided rank
                  !     col indexes the min strided rank
                  !in fortran packed(col, row) gives the expected value
                  call MapRowToColPacked(globalRow0, globalCol0, globalSize0, bands0, globalPackedRow0, globalPackedCol0)

                  localPackedRow0 = globalPackedRow0 - globalStartIndex0
                  localPackedCol0 = globalPackedCol0
             
                  !index in local packed array
                  i0 = localPackedRow0 * (2 * bands0 + 1) + localPackedCol0

                  !write(*,*) globalPackedRow0, globalPackedCol0, i0

                  !Matrix Vector Multiplication
                  tempValue(0) = tempValue(0) + potential(i0) * scaling * source(col0)
              enddo
          endif

          !----------------------------------------------------------------------------
          !                         Communication
          !----------------------------------------------------------------------------

          ! Theese are the sends and recvs from the previous step. wait for them
          ! now instead of in the previous iteration in order to try to overlap computation 
          ! and communication
          if (waitSend0 .eq. 1) then
              call MPI_Wait(sendRequest0, MPI_STATUS_IGNORE, error0)
          endif
          if (waitRecieve0 .eq. 1) then
              call MPI_Wait(recvRequest0, MPI_STATUS_IGNORE, error0)
              dest(sourceRow0) = dest(sourceRow0) + recvTemp0(0)
          endif

          ! Find out which procs to send and recv from
          destProc0 = GetOwnerProcId(globalSize0, procCount0, globalRow0)
          deltaProc0 = destProc0 - procId0
          sourceProc0 = procId0 - deltaProc0

          !Check if we're to recieve some data this interation
          waitRecieve0 = 0
          sourceGlobalStartIndex0 = GetGlobalStartIndex(globalSize0, procCount0, sourceProc0)
          sourceGlobalRow0 = row0 + sourceGlobalStartIndex0
          if (0.le.sourceGlobalRow0.and.sourceGlobalRow0.lt.globalSize0.and.0.le.sourceProc0.and.sourceProc0.lt.procCount0.and.deltaProc0.ne.0) then
              !write(*,*) "Recieving ", procId0, " <- ", sourceProc0
              call MPI_Irecv(recvTemp0, sendSize0, MPI_DOUBLE_COMPLEX, sourceProc0, tag0, communicator0, recvRequest0, error0)
              waitRecieve0 = 1
              !Calculate where we want to store this value
              sourceRow0 =  sourceGlobalRow0 - globalStartIndex0
          endif

          !If we have computed a value this iteration we need to store it
          !locally or send it to another proc
          waitSend0 = 0
          if (0.le.globalRow0.and.globalRow0.lt.globalSize0) then
              if (deltaProc0.eq.0) then
                  dest(row0) = dest(row0) + tempValue(0)
              else
                  !Store the value we're sending in a temp value
                  !TODO: Eliminiate this copy by alternating between
                  !temp arrays
                  sendTemp0(0) = tempValue(0)
                  call MPI_ISend(sendTemp0, sendSize0, MPI_DOUBLE_COMPLEX, destProc0, tag0, communicator0, sendRequest0, error0)
                  waitSend0 = 1
              endif
          endif


      enddo

      ! Theese are the sends and recvs from the last step. 
      if (waitSend0 .eq. 1) then
          call MPI_Wait(sendRequest0, MPI_STATUS_IGNORE, error0)
      endif
      if (waitRecieve0 .eq. 1) then
          call MPI_Wait(recvRequest0, MPI_STATUS_IGNORE, error0)
          dest(sourceRow0) = dest(sourceRow0) + recvTemp0(0)
      endif


end subroutine 
