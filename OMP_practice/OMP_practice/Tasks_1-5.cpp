#include <iostream>
#include "omp.h"
#include <iomanip>
#include <stdlib.h>
#include <cmath>
#include <ctime>
#include <list>

#define __RAND_MAX (unsigned int)(1 << 31)

using namespace std;

typedef double** mat;
typedef double* vec;

void MxV_byRows(mat A, vec H, int h, int w, vec result) {
#pragma omp parallel for schedule(static) shared(result)
	for (int i = 0; i < h; ++i) {
		double localRes = 0.0;

		for (int j = 0; j < w; ++j) {
			localRes += A[i][j] * H[j];
		}
		result[i] = localRes;
	}
}

void MxV_byCols(const mat A, const vec H, int h, int w, vec result) {
	for (int i = 0; i < h; ++i)
		result[i] = 0.0;

#pragma omp parallel for schedule(static) shared(result)
	for (int j = 0; j < w; ++j) {
		double Hj = H[j];
		for (int i = 0; i < h; ++i) {
#pragma omp atomic
			result[i] += A[i][j] * Hj;
		}
	}
}

void MxV_byBlocks(const mat A, const vec H, int h, int w, vec result) {
	for (int i = 0; i < h; ++i) {
		result[i] = 0.0;
	}

	int s, q = 1, tn, k, l;
#pragma omp parallel shared(s, q, tn, k, l)
	{
#pragma omp single
		tn = omp_get_num_threads();

#pragma omp barrier

#pragma omp for reduction(max: q)
		for (int i = 2; i <= (int)sqrt(tn); ++i) {
			if (tn % i == 0)
				q = i;
		}

#pragma omp single
		{
			if (w < h)
				s = tn / q;
			else {
				s = q;
				q = tn / s;
			}
			k = (h + h % s) / s;
			l = (w + w % q) / q;
		}
#pragma omp barrier

		int num = omp_get_thread_num();
		int bi = num / q;
		int bj = num % q;

		int bh = k;
		if (bi == s - 1)
			bh = h - bi * k;

		int bw = l;
		if (bj == q - 1)
			bw = w - bj * l;

		for (int i = bi * k; i < bi * k + bh; ++i) {
			double localRes = 0.0;

			for (int j = bj * l; j < bj * l + bw; ++j) {
				localRes += A[i][j] * H[j];
			}

#pragma omp atomic
			result[i] += localRes;
		}
	}
}

void MxM_Tape(mat A, mat B, int ha, int wa, int wb, mat result) {
#pragma omp parallel for schedule(static) shared(result)
	for (int n = 0; n < ha * wb; ++n) {
		int i = n / wb;
		int j = n % wb;
		double localRes = 0.0;
		for (int k = 0; k < wa; ++k) {
			localRes += A[i][k] * B[k][j];
		}
		result[i][j] = localRes;
	}
}

void MxM_Block(mat A, mat B, int ha, int wa, int wb, mat result) {
	int s, q = 1, tn, k, l;
#pragma omp parallel shared(s, q, tn, k, l)
	{
#pragma omp single
		tn = omp_get_num_threads();

#pragma omp barrier

#pragma omp for reduction(max: q)
		for (int i = 2; i <= (int)sqrt(tn); ++i) {
			if (tn % i == 0)
				q = i;
		}

#pragma omp single
		{
			if (wb < ha)
				s = tn / q;
			else {
				s = q;
				q = tn / s;
			}
			k = (ha + ha % s) / s;
			l = (wb + wb % q) / q;
		}
#pragma omp barrier

		int num = omp_get_thread_num();
		int bi = num / q;
		int bj = num % q;

		int bh = k;
		if (bi == s - 1)
			bh = ha - bi * k;

		int bw = l;
		if (bj == q - 1)
			bw = wb - bj * l;

		for (int i = bi * k; i < bi * k + bh; ++i) {
			for (int j = bj * l; j < bj * l + bw; ++j) {
				double localRes = 0.0;
				for (int k = 0; k < wa; ++k) {
					localRes += A[i][k] * B[k][j];
				}
				result[i][j] = localRes;
			}
		}
	}
}


void printMatrix(mat A, int h, int w) {
	for (int i = 0; i < h; ++i) {
		for (int j = 0; j < w; ++j) {
			cout << A[i][j] << " ";
		}
		cout << endl;
	}
}

void printVec(const vec V, int h) {
	for (int i = 0; i < h; ++i) {
		cout << V[i] << " ";
	}
	cout << endl;
}

mat generateMat(int h, int w, double mn, double mx) {
	mat A = new vec[h];
	for (int i = 0; i < h; ++i) {
		A[i] = new double[w];
		for (int j = 0; j < w; ++j) {
			A[i][j] = mn + ((double)rand() / (__RAND_MAX)) * (mx - mn);
		}
	}
	return A;
}

vec generateVec(int h, double mn, double mx) {
	vec newVec = new double[h];

	for (int i = 0; i < h; ++i) {
		newVec[i] = mn + ((double)rand() / (__RAND_MAX)) * (mx - mn);
	}
	return newVec;
}

void findAllEntries(const char* text, int length, const char* substr, int subLength, list<int>& positions) {
	int l, tn;
#pragma omp parallel shared(text, length, substr, subLength, positions)
	{
#pragma omp single
		{
			tn = omp_get_num_threads();
			l = (length - subLength + 1) / tn;
			if (0 == l) l = 1;
		}
#pragma omp barrier
		list<int> localResult;
		int num = omp_get_thread_num();
		for (int i = l * num; i < l * (num + 1) && i < length - subLength + 1; ++i) {
			int j = 0;
			while (
				j < subLength &&
				i < l * (num + 1) &&
				i < length - subLength + 1 &&
				text[i + j] == substr[j]
				) {
				++j;
			}

			if (j == subLength) {
				localResult.push_back(i);
			}
			else if (j != 0) {
				i += j - 1;
			}
		}
#pragma omp critical
		{
			positions.insert(positions.end(), localResult.begin(), localResult.end());
		}
	}
}

int MaxOfMinInRows(mat A, int h, int w) {
	int i, j;
	int max_e = -int(1e9);
	int min_e;
#pragma omp parallel
	{
#pragma omp for reduction(max: max_e) private(min_e, j)
		for (i = 0; i < h; ++i) {
			min_e = int(1e9);
			for (j = 0; j < w; ++j) {
				if (A[i][j] < min_e)
					min_e = A[i][j];
			}
			if (min_e > max_e)
				max_e = min_e;
		}
	}

	return max_e;
}

int main(int argc, char** argv) {
	srand(0);
	
	int m = atoi(argv[1]);
	int n = m;

	mat A = generateMat(m, n, -5, 5);

	const int rep = 10;
	double totalTime = 0.0;

	for (int i = 0; i < rep; ++i) {
		list<int> result;
		double start = omp_get_wtime();
		MaxOfMinInRows(A, m, n);
		double end = omp_get_wtime();
		totalTime += end - start;
	}
	totalTime /= rep;

	double ms = 1000 * totalTime;

	cout << ms;
	return 0;
}