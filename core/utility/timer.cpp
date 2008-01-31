#include "timer.h"
#include <mpi.h>

void Timer::Start()
{
	time -= MPI_Wtime();
}

void Timer::Stop()
{
	time += MPI_Wtime();
}


