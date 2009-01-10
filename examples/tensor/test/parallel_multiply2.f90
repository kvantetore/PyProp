module IndexTricks
contains
subroutine MapRowToColPacked(row, col, fullSize, bands, packedRow, packedCol)
    implicit none
    integer, intent(in) :: row, col, fullSize, bands
    integer, intent(out) :: packedRow, packedCol
    
    packedRow = col
    packedCol = bands - col + row
end subroutine

function GetLocalStartIndex(fullSize, procCount, procId) result(globalStartIndex)
    implicit none
    integer, intent(in) :: fullSize, procCount, procId
    integer :: rest, paddedSize, globalStartIndex, distribSize, firstSmallRank
    
    paddedSize    = fullSize
    rest = mod(fullSize, procCount)
    if (rest .ne. 0) then
        paddedSize = paddedSize + procCount - rest
    endif
    
    firstSmallRank = fullSize / paddedSize
    if (procId .le. firstSmallRank) then
        globalStartIndex = (paddedSize / procCount) * procId
        else
        globalStartIndex = fullSize - (procCount - procId)
    endif
    
    return
end function

function GetDistributedSize(fullSize, procCount, procId) result(distribSize)
    implicit none
    integer, intent(in) :: fullSize, procCount, procId
    integer :: rest, distribSize, paddedDistribSize
    
    rest = mod(fullSize, procCount)
    if (rest.eq.0) then
        distribSize = fullSize / procCount;
        else
        paddedDistribSize = (fullSize + procCount - rest) / procCount
        distribSize = fullSize - paddedDistribSize * procId
        
        if ( (distribSize .le. paddedDistribSize) .or. (distribSize .ge. 0) ) then
            distribSize = distribSize - (procCount - procId - 1)
        endif
        
        distribSize = max(distribSize, 1)
        distribSize = min(distribSize, paddedDistribSize)
    endif
    
    return
end function

function GetPaddedLocalSize(fullSize, procCount) result(localSize)
    implicit none
    integer, intent(in) :: fullSize, procCount
    integer :: rest, paddedSize, localSize
    
    paddedSize = fullSize
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
    integer :: rest, paddedSize, procId, localSize, firstSmallRank, distribPaddedSize
    
    include "parameters.f"
    
    paddedSize = fullSize
    rest = mod(fullSize, procCount)
    if (rest .ne. 0) then
        paddedSize = paddedSize + procCount - rest
    endif
    distribPaddedSize = paddedSize / procCount
    
    !round towards zero
    firstSmallRank = fullSize / paddedSize
    
    !round towards zero
    procId = globalIndex / distribPaddedSize
    
    !The procs after firstSmallRank have one datapoint each
    if (procId .gt. firstSmallRank) then
        procId = procCount - (fullSize - globalIndex)
    endif
    
    return
end function

end module



subroutine TensorPotentialMultiply_Distr2(potential, potentialExtent0, scaling, source, sourceExtent0, dest, destExtent0, &
                                globalSizeII, &
                                localMatrixIndexII, localMatrixIndexIIExtent0, &
                                globalRowII, globalRowIIExtent0, &
                                globalColII, globalColIIExtent0, &
                                sendProcII, sendProcIIExtent0, &
                                recvProcListII, recvProcListIIExtent0, recvProcListIIExtent1, &
                                recvLocalRowListII, recvLocalRowListIIExtent0, recvLocalRowListIIExtent1, &
                                recvCountII, recvCountIIExtent0)
    
    use IndexTricks
    implicit none
    include "mpif.h"
    include "parameters.f"
    
    integer, intent(in) :: potentialExtent0
    complex (kind=dbl), dimension(0:potentialExtent0-1), intent(in) :: potential
    real (kind=dbl) :: scaling
    
    integer, intent(in) :: sourceExtent0
    complex (kind=dbl), dimension(0:sourceExtent0-1), intent(in) :: source
    
    integer, intent(in) :: destExtent0
    complex (kind=dbl), dimension(0:destExtent0-1), intent(inout) :: dest
    
    integer, intent(in) :: globalSizeII

    integer, intent(in) :: localMatrixIndexIIExtent0
    integer, intent(in) :: globalRowIIExtent0
    integer, intent(in) :: globalColIIExtent0
    integer, intent(in) :: sendProcIIExtent0
    integer, intent(in) :: recvProcListIIExtent0, recvProcListIIExtent1
    integer, intent(in) :: recvLocalRowListIIExtent0, recvLocalRowListIIExtent1
    integer, intent(in) :: recvCountIIExtent0

    integer, dimension(0:localMatrixIndexIIExtent0), intent(in) :: localMatrixIndexII
    integer, dimension(0:globalRowIIExtent0), intent(in) :: globalRowII
    integer, dimension(0:globalColIIExtent0), intent(in) :: globalColII
    integer, dimension(0:sendProcIIExtent0), intent(in) :: sendProcII
    integer, dimension(0:recvProcListIIExtent1, 0:recvProcListIIExtent0), intent(in) :: recvProcListII
    integer, dimension(0:recvLocalRowListIIExtent1, 0:recvLocalRowListIIExtent0), intent(in) :: recvLocalRowListII
    integer, dimension(0:recvCountIIExtent0), intent(in) :: recvCountII

    integer :: iII, rowII, colII
    
    !temporary arrays TODO:FIX TEMPS
    complex (kind=dbl), dimension(0:255) :: recvTempII
    complex (kind=dbl), dimension(0:1) :: sendTempII
    integer, dimension(0:255) :: recvRequestII

    integer :: curSendII
    integer :: recvIdxII

    integer :: tempIndexII
    
    !Indices to map between local and global index ranges
    integer :: globalStartIndexII, paddedLocalSizeII
    
    !MPI variables
    integer :: errorII, tagII, sendSizeII
    integer :: sendRequestII
    integer :: waitRecieveII, waitSendII
    integer :: procIdII, procCountII, communicatorII
    
    !Indices for the recieved data
    integer :: sourceRowII
    
    call MPI_Comm_rank(MPI_COMM_WORLD, procIdII, errorII)
    call MPI_Comm_size(MPI_COMM_WORLD, procCountII, errorII)
    communicatorII = MPI_COMM_WORLD
    
    globalStartIndexII = GetLocalStartIndex(globalSizeII, procCountII, procIdII)
    paddedLocalSizeII = GetPaddedLocalSize(globalSizeII, procCountII)
    
    tagII = 0
    waitRecieveII = 0
    waitSendII = 0
    tempIndexII = 1
    curSendII = 0
    
    do iII = 0, localMatrixIndexIIExtent0-1
        !----------------------------------------------------------------------------
        !                         Calculation
        !----------------------------------------------------------------------------

         !Perform calculation for this step
        if (localMatrixIndexII(iII) .ne. -1) then
            colII = globalColII(iII) - globalStartIndexII
            rowII = globalRowII(iII) - globalStartIndexII

            if (sendProcII(iII) .eq. procIdII) then
                dest(rowII) = dest(rowII) + (potential(iII)) * scaling * source(colII)
                curSendII = 0
            else
                sendTempII(tempIndexII) = sendTempII(tempIndexII) + (potential(iII)) * scaling * source(colII)
                curSendII = 1
            endif
        endif

        !----------------------------------------------------------------------------
        !                         Communication
        !----------------------------------------------------------------------------
        
        ! Theese are the sends and recvs from the previous step. wait for them
        ! now instead of in the previous iteration in order to try to overlap computation
        ! and communication

        !Wait for previous recv
        do recvIdxII = 0, waitRecieveII-1
            call MPI_Wait(recvRequestII(recvIdxII), MPI_STATUS_IGNORE, errorII)
            sourceRowII = recvLocalRowListII(recvIdxII, iII-1)
            dest(sourceRowII) = dest(sourceRowII) + recvTempII(recvIdxII)
        enddo

        !Wait for previous send
        if (waitSendII .eq. 1) then
            call MPI_Wait(sendRequestII, MPI_STATUS_IGNORE, errorII)
        endif

        ! Theese are the sends and recvs for the current step. Post them now,
		! and wait for completion at the next iteration 

        !Post recvs for this step
        waitRecieveII = recvCountII(iII)
        do recvIdxII = 0, waitRecieveII-1
            call MPI_Irecv(recvTempII(recvIdxII), sendSizeII, MPI_DOUBLE_COMPLEX, recvProcListII(recvIdxII, iII), tagII, communicatorII, recvRequestII(recvIdxII), errorII)
        enddo 

        !Post send for this step if we have calculated anything
        waitSendII = 0
        if (curSendII .eq. 1) then
            call MPI_ISend(sendTempII(tempIndexII), sendSizeII, MPI_DOUBLE_COMPLEX, sendProcII(iII), tagII, communicatorII, sendRequestII, errorII)
            waitSendII = 1
        endif

        !Use the next tempIndex for the next computation
        tempIndexII = mod(tempIndexII+1, 1)
    enddo 

    !Wait for last recv
    do recvIdxII = 0, waitRecieveII-1
        call MPI_Wait(recvRequestII(recvIdxII), MPI_STATUS_IGNORE, errorII)
        sourceRowII = recvLocalRowListII(recvIdxII, iII-1)
        dest(sourceRowII) = dest(sourceRowII) + recvTempII(recvIdxII)
    enddo

    !Wait for last send
    if (waitSendII .eq. 1) then
        call MPI_Wait(sendRequestII, MPI_STATUS_IGNORE, errorII)
    endif

end subroutine TensorPotentialMultiply_Distr2


