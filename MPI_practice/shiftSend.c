#include <mpi.h>
#include <stdio.h>

/*
	В данном файле приведен код, представляющих реализацию
	межпроцессного взаимодействия, при котором данные передаются по кругу одновременно:
	каждый процесс отправляет и получает данные
*/

int main(int argc, char** argv) {
    MPI_Init(NULL, NULL);

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

	const int nextProcess = (world_rank + 1) % world_size;
	const int prevProcess = (world_rank - 1) % world_size;
	int sum = world_rank;

	MPI_Request requestMock;

	// процесс отправляет данные следующему процессу
	// получает данные от предыдущего процесса
	// это происходит n - 1 раз, а значит, данные проходят через все процессы
	for(int i = 0; i < world_size - 1; ++i) {
		MPI_Isend(&sum, 1, MPI_INT, nextProcess, 0, MPI_COMM_WORLD, &requestMock);
		MPI_Recv(&sum, 1, MPI_INT, prevProcess, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		sum += world_rank;
	}

	printf("Process %d: total sum is %d\n", world_rank, sum);

    MPI_Finalize();
}
