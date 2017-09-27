#include<omp.h>
#include<iostream>
#include<time.h>

using namespace std;

const int REQUESTED_THREADS = 4;

double slow_sum(double* arr, int num) {
		double* dASum = new double[REQUESTED_THREADS];

		double dSum = 0;

		for (int j = 0; j < REQUESTED_THREADS; j++) {
			dASum[j] = 0.0;
		}

		omp_set_num_threads(REQUESTED_THREADS);

		#pragma omp parallel
		{
			int i, numThreads, threadID;

			threadID = omp_get_thread_num();

			numThreads = omp_get_num_threads();

			for (i = threadID; i < num; i += numThreads) {

				dASum[threadID] += arr[i];

			}
		}

		for (int j = 0; j < REQUESTED_THREADS; j++) {
			dSum += dASum[j];
		}

		delete[] dASum;

		//cout << numberOfThreads << endl;

		return dSum;
}

double fast_sum(double* A, int num) {

	double sum = 0;

	double *B = new double[num * (log2(num) + 1)];
	double *C = new double[num * (log2(num) + 1)];

#pragma omp parallel for
	for (int j = 0; j < num; j++) {
		B[j] = A[j];
	}

	int logNum = log2(num);

	for (int h = 0; h <= logNum; h++) {
#pragma omp parallel for
		for (int j = 0; j <= (num / 2 ^ h); j++) {
			B[num*h + j] = B[num*(h - 1) + (2 * j - 1)] + B[num*(h - 1) + (2 * j)];
		}
	}

	for (int h = logNum; h >= 0; h--) {
#pragma omp parallel for
		for (int j = 0; j <= (num / 2 ^ h); j++) {
			if (j == 0) {
				C[num*h] = B[num*h];
			}
			switch (j % 2) {
			case 0:
				C[num*h + j] = C[num*(h + 1) + (j / 2)];
				break;
			case 1:
				C[num*h + j] = C[num*(h + 1) + ((j - 1) / 2)];
				break;
			}
		}
	}

	/*
#pragma omp parallel for
	for (int i = 0; i < num; i++) {
		sum += arr[i];
	}
	*/
	return C[num];
}


int main() {
	const int i_N = 1000000;
	const int i_R = 1000;
	double* inputArr = new double[i_N+1];
	double sum = 0.0;
	clock_t t;
	double runTime, mflops;

	for (int i = 0; i < i_N+1; i++) {
		inputArr[i] = i;
	}

	t = clock();
	for (int i = 0; i <= i_R; i++) {
		sum = slow_sum(inputArr, i_N + 1);
	}

	runTime = 1.0*(clock() - t) / CLOCKS_PER_SEC;

	mflops = (i_R*i_N) / (runTime*1e6);

	cout << "Slow Sum MFLOPS: " << mflops << endl;// ':' << runTime << endl;

	cout << sum << endl;




	t = clock();
	for (int i = 0; i <= i_R; i++) {
		sum = fast_sum(inputArr, i_N + 1);
	}

	runTime = 1.0*(clock() - t) / CLOCKS_PER_SEC;

	mflops = (i_R*i_N) / (runTime*1e6);

	cout << "Fast Sum MFLOPS: " << mflops << endl;// ':' << runTime << endl;

	cout << sum << endl;


	delete[] inputArr;
	return 0;
}