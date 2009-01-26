
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
    integer :: rest, paddedSize, globalStartIndex, distribSize, firstSmallRank, distrPaddedSize
    
    paddedSize	= fullSize
    rest = mod(fullSize, procCount)
    if (rest .ne. 0) then
        paddedSize = paddedSize + procCount - rest
    endif
    distrPaddedSize = paddedSize/procCount
    
    firstSmallRank = fullSize / distrPaddedSize
    if (mod(fullSize, distrPaddedSize) .eq. 0) then
        firstSmallRank = firstSmallRank - 1
    endif
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
        
        if ( (distribSize .le. paddedDistribSize) .and. (distribSize .gt. 0) ) then
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
    integer :: rest, paddedSize, procId, firstSmallRank, distribPaddedSize
    
    include "parameters.f"
    
    paddedSize = fullSize
    rest = mod(fullSize, procCount)
    if (rest .ne. 0) then
        paddedSize = paddedSize + procCount - rest
    endif
    distribPaddedSize = paddedSize / procCount
    
    !round towards zero
    firstSmallRank = fullSize / distribPaddedSize
    
    !round towards zero
    procId = globalIndex / distribPaddedSize
    
    !The procs after firstSmallRank have one datapoint each
    if (procId .gt. firstSmallRank) then
        procId = procCount - (fullSize - globalIndex)
    endif
    
    return
end function

end module





subroutine TensorPotentialMultiply_SimpD_Simp_BandNH(potential, potentialExtent0, potentialExtent1, potentialExtent2, scaling, source, sourceExtent0, sourceExtent1, sourceExtent2, dest, destExtent0, destExtent1, destExtent2, globalSize0, localMatrixIndex0, localMatrixIndex0Extent0, globalRow0, globalRow0Extent0, globalCol0, globalCol0Extent0, sendProc0, sendProc0Extent0, recvProcList0, recvProcList0Extent0, recvProcList0Extent1, recvLocalRowList0, recvLocalRowList0Extent0, recvLocalRowList0Extent1, recvCount0, recvCount0Extent0, recvTemp0, recvTemp0Extent0, recvTemp0Extent1, recvTemp0Extent2, sendTemp0, sendTemp0Extent0, sendTemp0Extent1, sendTemp0Extent2, pair1, pair1Extent0, pair1Extent1)
    
    use IndexTricks
    implicit none
    include "mpif.h"
    include "parameters.f"
    
    integer, intent(in) :: potentialExtent2, potentialExtent1, potentialExtent0
    complex (kind=dbl), dimension(0:potentialExtent2-1, 0:potentialExtent1-1, 0:potentialExtent0-1), intent(in) :: potential
    real (kind=dbl) :: scaling
    
    integer, intent(in) :: sourceExtent2, sourceExtent1, sourceExtent0
    complex (kind=dbl), dimension(0:sourceExtent2-1, 0:sourceExtent1-1, 0:sourceExtent0-1), intent(in) :: source
    
    integer, intent(in) :: destExtent2, destExtent1, destExtent0
    complex (kind=dbl), dimension(0:destExtent2-1, 0:destExtent1-1, 0:destExtent0-1), intent(inout) :: dest
    
    integer, intent(in) :: localMatrixIndex0Extent0
    integer, dimension(0:localMatrixIndex0Extent0-1), intent(in) :: localMatrixIndex0
    
    integer, intent(in) :: globalRow0Extent0
    integer, dimension(0:globalRow0Extent0-1), intent(in) :: globalRow0
    
    integer, intent(in) :: globalCol0Extent0
    integer, dimension(0:globalCol0Extent0-1), intent(in) :: globalCol0
    
    integer, intent(in) :: sendProc0Extent0
    integer, dimension(0:sendProc0Extent0-1), intent(in) :: sendProc0
    
    integer, intent(in) :: recvProcList0Extent1, recvProcList0Extent0
    integer, dimension(0:recvProcList0Extent1-1, 0:recvProcList0Extent0-1), intent(in) :: recvProcList0
    
    integer, intent(in) :: recvLocalRowList0Extent1, recvLocalRowList0Extent0
    integer, dimension(0:recvLocalRowList0Extent1-1, 0:recvLocalRowList0Extent0-1), intent(in) :: recvLocalRowList0
    
    integer, intent(in) :: recvCount0Extent0
    integer, dimension(0:recvCount0Extent0-1), intent(in) :: recvCount0
    
    integer, intent(in) :: recvTemp0Extent2, recvTemp0Extent1, recvTemp0Extent0
    complex (kind=dbl), dimension(0:recvTemp0Extent2-1, 0:recvTemp0Extent1-1, 0:recvTemp0Extent0-1), intent(inout) :: recvTemp0
    
    integer, intent(in) :: sendTemp0Extent2, sendTemp0Extent1, sendTemp0Extent0
    complex (kind=dbl), dimension(0:sendTemp0Extent2-1, 0:sendTemp0Extent1-1, 0:sendTemp0Extent0-1), intent(inout) :: sendTemp0
    
    integer, intent(in) :: globalSize0
    integer :: i0, row0, col0
    
    !temporary arrays TODO:FIX TEMPS <= FIXED! :D
    !complex (kind=dbl), dimension(0:destExtent2-1, 0:destExtent1-1, 0:recvProcList0Extent1) :: recvTemp0
    !complex (kind=dbl), dimension(0:destExtent2-1, 0:destExtent1-1, 0:1) :: sendTemp0
    integer :: tempIndex0
    
    integer :: curSend0
    integer :: recvIdx0
    
    !Indices to map between local and global index ranges
    integer :: globalStartIndex0 !!!!, paddedLocalSize0, localSize0
    
    !MPI variables
    integer, dimension(0:recvProcList0Extent1) :: recvRequest0
    integer :: error0, tag0, sendSize0
    integer :: sendRequest0
    integer :: waitRecieve0, waitSend0
    integer :: procId0, procCount0, communicator0
    
    !Indices for the recieved data
    integer :: sourceRow0
    
    integer, intent(in) :: pair1Extent1, pair1Extent0
    integer, dimension(0:pair1Extent1-1, 0:pair1Extent0-1), intent(in) :: pair1
    
    integer :: N1, i1, row1, col1
    
    integer :: N2, i2, row2, col2, bsplineCount2, bandCount2
    integer :: subDiagonals2, sourceStride2, destStride2
    complex (kind=dbl) :: alpha2, beta2
    
    call MPI_Comm_rank(MPI_COMM_WORLD, procId0, error0)
    call MPI_Comm_size(MPI_COMM_WORLD, procCount0, error0)
    communicator0 = MPI_COMM_WORLD
    
    globalStartIndex0 = GetLocalStartIndex(globalSize0, procCount0, procId0)
    !paddedLocalSize0 = GetPaddedLocalSize(globalSize0, procCount0)
    
    sourceRow0 = -1
    sendSize0 = 1 * destExtent1 * destExtent2
    tag0 = 0
    tempIndex0 = 1
    curSend0 = 0
    
    N1 = potentialExtent1
    
    N2 = potentialExtent2
    bsplineCount2 = sourceExtent2
    bandCount2 = N2 / bsplineCount2
    
    subDiagonals2 = (bandCount2 - 1)/2
    sourceStride2 = 1
    destStride2 = 1
    alpha2 = scaling
    beta2 = 1.0d0
    col2 = 0
    row2 = 0
    i2 = 0
    
    
    waitRecieve0 = 0
    waitSend0 = 0
    tempIndex0 = 1
    
    !Iterate over all rows of the matrix for the columns stored on proc
    do i0 = 0, localMatrixIndex0Extent0-1	
		!call MPI_Barrier(MPI_COMM_WORLD, error0);
        write(*,*) "Procid ", procId0, " step ", i0
        !----------------------------------------------------------------------------
        !                         Calculation
        !----------------------------------------------------------------------------
        
        !Perform calculation for this step
        curSend0 = 0
        if (localMatrixIndex0(i0) .ne. -1) then
            col0 = globalCol0(i0) - globalStartIndex0
            row0 = globalRow0(i0) - globalStartIndex0
            
            if (sendProc0(i0) .eq. procId0) then
                
                do i1 = 0, N1-1
                    row1 = pair1(0, i1);
                    col1 = pair1(1, i1);
                    
                    i2 = 0
                    call zgbmv( &
                    "N", &
                    bsplineCount2, &
                    bsplineCount2, &
                    subDiagonals2, &
                    subDiagonals2, &
                    alpha2, &
                    potential(0, i1, i0), &
                    bandCount2, &
                    source(0, col1, col0), &
                    sourceStride2, &
                    beta2, &
                    dest(0, row1, row0), &
                    destStride2 &
                    )
                    
                enddo
                
            else
                sendTemp0(:, :, tempIndex0) = 0
                
                do i1 = 0, N1-1
                    row1 = pair1(0, i1);
                    col1 = pair1(1, i1);
                    
                    i2 = 0
                    call zgbmv( &
                    "N", &
                    bsplineCount2, &
                    bsplineCount2, &
                    subDiagonals2, &
                    subDiagonals2, &
                    alpha2, &
                    potential(0, i1, i0), &
                    bandCount2, &
                    source(0, col1, col0), &
                    sourceStride2, &
                    beta2, &
                    sendTemp0(0, row1, tempIndex0), &
                    destStride2 &
                    )
                    
                enddo
                
                curSend0 = 1
            endif
        endif
        
        !----------------------------------------------------------------------------
        !                         Communication
        !----------------------------------------------------------------------------
        
        ! Theese are the sends and recvs from the previous step. wait for them
        ! now instead of in the previous iteration in order to try to overlap computation
        ! and communication
        
        !Wait for previous recv
        do recvIdx0 = 0, waitRecieve0-1
            call MPI_Wait(recvRequest0(recvIdx0), MPI_STATUS_IGNORE, error0)
            sourceRow0 = recvLocalRowList0(recvIdx0, i0-1)
            dest(:, :, sourceRow0) = dest(:, :, sourceRow0) + recvTemp0(:, :, recvIdx0)
        enddo
        
        !Wait for previous send
        if (waitSend0 .eq. 1) then
            call MPI_Wait(sendRequest0, MPI_STATUS_IGNORE, error0)
        endif
        
        ! Theese are the sends and recvs for the current step. Post them now,
        ! and wait for completion at the next iteration
        
        !Post recvs for this step
        waitRecieve0 = recvCount0(i0)
        do recvIdx0 = 0, waitRecieve0-1
            call MPI_Irecv(recvTemp0(0, 0, recvIdx0), sendSize0, MPI_DOUBLE_COMPLEX, recvProcList0(recvIdx0, i0), tag0, communicator0, recvRequest0(recvIdx0), error0)
        enddo
        
        !Post send for this step if we have calculated anything
        waitSend0 = 0
        if (curSend0 .eq. 1) then
            call MPI_ISend(sendTemp0(0, 0, tempIndex0), sendSize0, MPI_DOUBLE_COMPLEX, sendProc0(i0), tag0, communicator0, sendRequest0, error0)
            waitSend0 = 1
        endif
        
        !Use the next tempIndex for the next computation
        tempIndex0 = mod(tempIndex0+1, 2)
    enddo
    
    !Wait for last recv
    do recvIdx0 = 0, waitRecieve0-1
        call MPI_Wait(recvRequest0(recvIdx0), MPI_STATUS_IGNORE, error0)
        sourceRow0 = recvLocalRowList0(recvIdx0, localMatrixIndex0Extent0-1)
        dest(:, :, sourceRow0) = dest(:, :, sourceRow0) + recvTemp0(:, :, recvIdx0)
    enddo
    
    !Wait for last send
    if (waitSend0 .eq. 1) then
        call MPI_Wait(sendRequest0, MPI_STATUS_IGNORE, error0)
    endif
    
    
end subroutine TensorPotentialMultiply_SimpD_Simp_BandNH


