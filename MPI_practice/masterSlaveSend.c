#include <mpi.h>
#include <stdio.h>

/*
	В данном файле приведен код, представляющих реализацию
	межпроцессного взаимодействия типа master-slave
*/

int main(int argc, char** argv) {
    MPI_Init(NULL, NULL);

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

	if (0 == world_rank) {
		// код для главного процесса
		// главный процесс ждем данных от всех процессов
		int sum = 0;
		for(int i = 1; i < world_size; ++i) {
			int recievedData;
			MPI_Recv(&recievedData, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			sum += recievedData;
		}
		printf("Sum = %d\n", sum);
	}
	else {
		// отправка данных главномц процессу
		MPI_Send(&world_rank, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
	}

    MPI_Finalize();
}
