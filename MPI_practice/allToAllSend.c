#include <mpi.h>
#include <stdio.h>

/*
	В данном файле приведен код, представляющих реализацию
	передачи данных от каждого процесса каждому с помощью связи
	тип точка-точка
*/

int main(int argc, char** argv) {
    MPI_Init(NULL, NULL);

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

	int sum = world_rank;
	MPI_Request requestMock;

	// неблокирующая отправка всем процессам
	for(int i = 0; i < world_size; ++i) {
		if (i == world_rank) continue;
		MPI_Isend(&world_rank, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &requestMock);
	}

	// блокирующее получение результата от всех процессов
	for(int i = 0; i < world_size; ++i) {
		if (i == world_rank) continue;
		int receivedData;
		MPI_Recv(&receivedData, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		sum += receivedData;
	}

	printf("Process %d: total sum is %d\n", world_rank, sum);

    MPI_Finalize();
}
