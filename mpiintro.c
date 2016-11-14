#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#include <unistd.h> 

#define N 8
#define M 8
#define DEPTH 10
#define TAG_PART_READY 0
#define TAG_PREV_READY 1000
#define TAG_NEXT_READY 2000
#define TAG_STEPS 500

int min(int x, int y) {
	if (x < y) {
		return x;
	}
	return y;
}

void step_cell(int x, int y, int level, int mp[][N][M]) {
	int i, j;
	int curx, cury;
	int neighbours = 0;
	for (int i = -1; i <= 1; i++) {
		for (int j = -1; j <= 1; j++) {
			if (i == 0 && j == 0) {
				continue;
			}
			curx = x + i;
			cury = y + j;
			if (curx < 0) {
				curx = N - 1;
			} else if (curx == N) {
				curx = 0;
			}
			if (cury < 0) {
				cury = M - 1;
			} else if (cury == M) {
				cury = 0;
			}
			if (mp[level - 1][curx][cury] == 1) {
				neighbours++;
			}
		}
	}
	if (mp[level - 1][x][y] == 0 && neighbours == 3) {
		mp[level][x][y] = 1;
	} else if (mp[level - 1][x][y] == 1 && (neighbours == 2 || neighbours == 3)) {
		mp[level][x][y] = 1;
	} else {
		mp[level][x][y] = 0;
	}
}

void calculate_one_frame(int frame, int level, int rank,
						 int size, int mp[][N][M]) {
	int min_y = rank * N / size + frame;
	int max_y = min(N - 1, (rank + 1) * N / size - 1) - frame;
	int i, j;
	if (min_y > max_y) {
		return;
	}
	// printf("frame: %d, level: %d, rank: %d, min_y: %d, max_y: %d\n", frame, level, rank, min_y, max_y);
	if (level == DEPTH) {
		for (i = min_y; i <= max_y; i++) {
			for (j = 0; j < M; j++) {
				step_cell(i, j, level, mp);
			}
		}
		return;
	}
	for (j = frame; j < M - frame; j++) {
		step_cell(min_y, j, level, mp);
		step_cell(max_y, j, level, mp);
	}
	for (i = min_y + 1; i < max_y; i++) {
		step_cell(i, frame, level, mp);
		step_cell(i, M - 1 - frame, level, mp);
	}
}

void calculate_pyramid(int rank, int size, int mp[][N][M]) { // ERROR!!!
	int i, j;
	for (i = 0; i < DEPTH - 1; i++) {
		for (j = i; j < (M + 1) / 2; j++) {
			calculate_one_frame(j, i + 1, rank, size, mp);
		}
	}
	calculate_one_frame(DEPTH - 1, DEPTH, rank, size, mp);
}

void add_layer(int level, int rank, int size, int mp[][N][M]) {
	int i;
	for (i = level; i < DEPTH; i++) {
		calculate_one_frame(i - level, i + 1, rank, size, mp);
	}
}

int main (int argc, char **argv) {
	int size, prev_rk, next_rk;
	MPI_Status status;
	int i, j, k;
	int mp[DEPTH + 1][N][M];
	int steps;
	int lol[N];
	int zero_ready = 0;
	int one_ready = 0;
	int rank = 0;
	int min_y, max_y;

	MPI_Init(&argc, &argv);
	
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	prev_rk = (rank + size - 1) % size;
	next_rk = (rank + 1) % size;

	min_y = rank * N / size;
	max_y = min(N - 1, (rank + 1) * N / size - 1);

	if (rank == 0) {

		// INITIALIZING

		for (i = 0; i < N; i++) {
			for (j = 0; j < M; j++) {
				mp[0][i][j] = 0;
			}
		}

		FILE *fi = fopen("input.txt", "r");

		fscanf(fi, "%d", &steps);

		for (i = 0; i < N; i++) {
			for (j = 0; j < M; j++) {
				fscanf(fi, "%d", &mp[0][i][j]);
			}
		}

		printf("-begin-\n");
		for (i = 0; i < N; i++) {
			for (j = 0; j < M; j++) {
				printf("%d ", mp[0][i][j]);
			}
			printf("\n");
		}

		for (i = 1; i < size - 1; i++) {
			MPI_Send(&mp[0][i * N / size - 1][0], N / size * M + 2 * M, MPI_INT, 
				i, TAG_PART_READY + i, MPI_COMM_WORLD);
		}

		MPI_Send(&mp[0][i * N / size - 1][0], N * M - N / size * M * i + M, MPI_INT, 
			i, TAG_PART_READY + i, MPI_COMM_WORLD);

		MPI_Send(&mp[0][0][0], M, MPI_INT,
			i, TAG_PART_READY + i, MPI_COMM_WORLD);

		for (i = 1; i < size; i++) {
			MPI_Send(&steps, 1, MPI_INT, i, TAG_STEPS + i, MPI_COMM_WORLD);
		}

	} else {


		if (rank != size - 1) {
			MPI_Recv(&mp[0][rank * N / size - 1][0], N / size * M + 2 * M, MPI_INT, 
				0, TAG_PART_READY + rank, MPI_COMM_WORLD, NULL);
		} else {
			MPI_Recv(&mp[0][rank * N / size - 1][0], N * M - N / size * M * (size - 1) + M, MPI_INT, 
				0, TAG_PART_READY + rank, MPI_COMM_WORLD, NULL);
			MPI_Recv(&mp[0][0][0], M, MPI_INT,
				0, TAG_PART_READY + rank, MPI_COMM_WORLD, NULL);
		}


		MPI_Recv(&steps, 1, MPI_INT, 0, TAG_STEPS + rank, MPI_COMM_WORLD, NULL);

	}

	// if (rank == 0) {
	// 	sleep(2);
	// }
	
	// printf("my rank: %d\n", rank);

	// for (i = 0; i < N; i++) {
	// 	for (j = 0; j < M; j++) {
	// 		printf("%d ", mp[0][i][j]);
	// 	}
	// 	printf("\n");
	// }

	do {
		// printf("%d\n", steps);

		calculate_pyramid(rank, size, mp);


		// for (i = min_y; i <= max_y; i++) {
		// 	for (j = 0; j < M; j++) {
		// 		printf("rank %d: %d ", rank, mp[1][i][j]);
		// 	}

		// 	printf("\n");
		// }

		MPI_Send(&mp[1][rank * N / size][0], M, MPI_INT, 
			prev_rk, TAG_NEXT_READY, MPI_COMM_WORLD);

		MPI_Send(&mp[1][min(N - 1, (rank + 1) * N / size - 1)][0], M, MPI_INT, 
			next_rk, TAG_PREV_READY, MPI_COMM_WORLD);

		for (i = 1; i < DEPTH; i++) {

			MPI_Recv(&mp[i][min(N - 1, (prev_rk + 1) * N / size - 1)][0], M, MPI_INT, 
				prev_rk, TAG_PREV_READY + i - 1, MPI_COMM_WORLD, NULL);

			MPI_Recv(&mp[i][next_rk * N / size][0], M, MPI_INT, 
				next_rk, TAG_NEXT_READY + i - 1, MPI_COMM_WORLD, NULL);

			add_layer(i, rank, size, mp);

			if (i + 1 == DEPTH) {
				break;
			}

			MPI_Send(&mp[i + 1][rank * N / size][0], M, MPI_INT, 
				prev_rk, TAG_NEXT_READY + i, MPI_COMM_WORLD);

			MPI_Send(&mp[i + 1][min(N - 1, (rank + 1) * N / size - 1)][0], M, MPI_INT, 
				next_rk, TAG_PREV_READY + i, MPI_COMM_WORLD);
		}

		for (j = 0; j < M; j++) {
			mp[0][min_y][j] = mp[DEPTH][min_y][j];
			mp[0][max_y][j] = mp[DEPTH][max_y][j];
		}

		MPI_Send(&mp[0][min(N - 1, (rank + 1) * N / size - 1)][0], M, MPI_INT,
			next_rk, TAG_PREV_READY, MPI_COMM_WORLD);

		MPI_Send(&mp[0][rank * N / size][0], M, MPI_INT,
			prev_rk, TAG_NEXT_READY, MPI_COMM_WORLD);

		for (i = min_y + 1; i < max_y; i++) {
			for (j = 0; j < M; j++) {
				mp[0][i][j] = mp[DEPTH][i][j];
			}
		}

		MPI_Recv(&mp[0][min(N - 1, (prev_rk + 1) * N / size - 1)][0], M, MPI_INT, 
			prev_rk, TAG_PREV_READY, MPI_COMM_WORLD, NULL);

		MPI_Recv(&mp[0][next_rk * N / size][0], M, MPI_INT, 
			next_rk, TAG_NEXT_READY, MPI_COMM_WORLD, NULL);

		steps -= DEPTH;
	} while (steps > 0);

	steps += DEPTH;

	// for (i = min_y; i <= max_y; i++) {
	// 	for (j = 0; j < M; j++) {
	// 		printf("rank %d: %d ", rank, mp[2][i][j]);
	// 	}
	// 	printf("\n");
	// }

	if (rank != 0) {
		if (rank != size - 1) {
			MPI_Send(&mp[steps][rank * N / size][0], N / size * M, MPI_INT, 
				0, TAG_PART_READY + rank, MPI_COMM_WORLD);
		} else {
			MPI_Send(&mp[steps][rank * N / size][0], N * M - N / size * M * (size - 1), MPI_INT, 
				0, TAG_PART_READY + rank, MPI_COMM_WORLD);
		}
		return 0;
	} else {
		for (i = 1; i < size - 1; i++) {
			MPI_Recv(&mp[steps][i * N / size][0], N / size * M, MPI_INT, 
				i, TAG_PART_READY + i, MPI_COMM_WORLD, NULL);
		}
		MPI_Recv(&mp[steps][i * N / size][0], N * M - N / size * M * (size - 1), MPI_INT, 
			i, TAG_PART_READY + i, MPI_COMM_WORLD, NULL);
		printf("\n--end--\n");
		for (i = 0; i < N; i++) {
			for (j = 0; j < M; j++) {
				printf("%d ", mp[steps][i][j]);
			}
			printf("\n");
		}

		FILE *fo = fopen("output.txt", "w");
		for (i = 0; i < N; i++) {
			for (j = 0; j < M; j++) {
				fprintf(fo, "%d ", mp[steps][i][j]);
			}
			fprintf(fo, "\n");
		}
	}

	MPI_Finalize();

	return 0;
}