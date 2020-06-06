#include <mpi.h>
#include <stdio.h>

/*
	В данном файле приведен код, представляющих реализацию
	межпроцессного взаимодействия, при котором данные передаются по кругу последовательно:
	один процесс получает данные от предыдущего и отправляет данные следующему
*/

int main(int argc, char** argv) {
    MPI_Init(NULL, NULL);

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

	const int nextProcess = (world_rank + 1) % world_size;
	const int prevProcess = (world_rank - 1) % world_size;
	int sum = 0;
	
	// все процессы, кроме самого первого сначала ожидают получения данных
	if (world_rank != 0) {
		MPI_Recv(&sum, 1, MPI_INT, prevProcess, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}

	// затем все процессы отправляют данные
	sum += world_rank;
	MPI_Send(&sum, 1, MPI_INT, nextProcess, 0, MPI_COMM_WORLD);

	// самый первый процесс получает данные в самом конце, чтобы замкнуть передачу эстафеты
	if (0 == world_rank) {
		MPI_Recv( &sum, 1, MPI_INT, prevProcess, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		printf("Sum = %d\n", sum);
	}

    MPI_Finalize();
}
