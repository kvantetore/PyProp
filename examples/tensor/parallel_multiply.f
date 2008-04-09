integer function MapRowColToPacked(row, col, fullSize, bands, packedRow, packedCol)
      implicit none
      integer, intent(in) :: row, col, fullSize, bands
      integer, intent(out) :: packedRow, packedCol
      packedRow = col
      packedCol = bands - col + row
      return packedRow, packedCol

integer function GetGlobalStartIndex(fullSize, procCount, procId):
      rest = fullSize % procCount
      if rest != 0:
          fullSize += procCount - rest
      return (fullSize / procCount) * procId

integer function GetDistributedShape(fullSize, procCount, procId):
      rest = fullSize % procCount
      if rest == 0:
          distrShape = fullSize / procCount;
      else:
          paddedDistrShape = (fullSize + procCount - rest) / procCount 
          shape = fullSize - paddedDistrShape * procId
          shape = max(shape, 0)
          shape = min(shape, paddedDistrShape)
          distrShape = shape
      return distrShape;

def GetOwnerProcId(fullSize, procCount, globalIndex):
      rest = fullSize % procCount
      if rest != 0:
          fullSize += procCount - rest
      
      return globalIndex / (fullSize / procCount)
      
      



subroutine BandedMatrixMultiply(potential, potentialExtent0, scaling, source, sourceExtent0, dest, destExtent0, procCount0, procId0, communicator0 )

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
      
      integer, intent(in) :: procCount0, procId0, communicator0

      integer :: N0, i0, row0, col0, bands0
      
      complex (kind=dbl), dimension(0:0) :: tempValue, recvTemp0
      integer :: deltaProc0, sourceProc0, destProc0, globalRow0, error0, globalPackedRow0, globalPackedCol0, localPackedRow0, localPackedCol0, localStartIndex0, localEndIndex0, globalSize0, sourceRow0, tag0, sendSize0, localSize0, globalCol0, recvRequest0, globalStartIndex0


      N0 = potentialExtent0
      localSize0 = sourceExtent0
      bands0 = ((N0 / localSize0) - 1) / 2

      globalStartIndex0 = procId0 * localSize0
      globalSize0 = localSize0 * procCount0

      sendSize0 = 1
      tag0 = 0

      col0 = 0
      row0 = 0
      i0 = 0

      do row0 = -bands0, localSize0+bands0-1

          tempValue(0) = 0
          globalRow0 = row0 + globalStartIndex0
          if (0.le.globalRow0.and.globalRow0.lt.globalSize0) then
              localStartIndex0 = max(0, row0 - bands0)
              localEndIndex0 = min(row0 + bands0 + 1, localSize0)
              do col0 = localStartIndex0, localEndIndex0-1
                  globalCol0 = col0 + globalStartIndex0
                  
                  !PackedRows and -Cols are for C-style arrays
                  !i.e. row indexes the max strided rank
                  !     col indexes the min strided rank
                  !in fortran packed(col, row) gives the expected value
                  globalPackedRow0 = globalCol0
                  globalPackedCol0 = bands0 - globalCol0 + globalRow0

                  localPackedRow0 = globalPackedRow0 - globalStartIndex0
                  localPackedCol0 = globalPackedCol0
             
                  !index in local packed array
                  i0 = localPackedRow0 * (2 * bands0 + 1) + localPackedCol0

                  !write(*,*) globalPackedRow0, globalPackedCol0, i0

                  !Matrix Vector Multiplication
                  tempValue(0) = tempValue(0) + potential(i0) * scaling * source(col0)
              enddo
          endif

          deltaProc0 = floor( real(row0, kind=dbl) / real(localSize0, kind=dbl) )
          destProc0 = procId0 + deltaProc0
          sourceProc0 = procId0 - deltaProc0

          if (deltaProc0.eq.0) then
              dest(row0) = dest(row0) + tempValue(0)
          else

              if (0.le.sourceProc0.and.sourceProc0.lt.procCount0) then
                  call MPI_Irecv(recvTemp0, sendSize0, MPI_DOUBLE_COMPLEX, sourceProc0, tag0, communicator0, recvRequest0, error0)
              endif
              if (0.le.destProc0.and.destProc0.lt.procCount0) then
                  call MPI_Send(tempValue, sendSize0, MPI_DOUBLE_COMPLEX, destProc0, tag0, communicator0, error0)
              endif
              if (0.le.sourceProc0.and.sourceProc0.lt.procCount0) then
                  call MPI_Wait(recvRequest0, MPI_STATUS_IGNORE, error0)
                    
                  sourceRow0 = row0 + localSize0 * (sourceProc0 - procId0) 
                  dest(sourceRow0) = dest(sourceRow0) + recvTemp0(0)
              endif
          endif
      enddo

end subroutine 
