#include<omp.h>
#include<iostream>
#include<time.h>
#include<sstream>
#include<fstream>
#include<algorithm>
#include<iomanip>

using namespace std;

int REQUESTED_THREADS = 4;
//ofstream outArr;

double slow_sum(double* arr, int num, int depth = 0) {

	int dASumSize = (num + 1) / 2; // equivalent of ceil (num/2) without floating point division
	int maxThreads = num / 2;

	double* dASum = new double[dASumSize];

	double sum = 0.0;

	omp_set_num_threads(min(maxThreads, REQUESTED_THREADS));

	#pragma omp parallel
	{
		int i, iEnd, numThreads, threadID;
		threadID = omp_get_thread_num();

		numThreads = omp_get_num_threads();

		int numDivNumThreads = ((num + 1) & ~1) / numThreads;

		i = numDivNumThreads * threadID;
		iEnd = numDivNumThreads * (threadID + 1);

		i = (i + 1) & ~1; // Round i to even that is greater than itself

		if (threadID == (numThreads - 1)) {
			iEnd = num - 1;
			//cout << numThreads << endl;
		}

		for (; i < iEnd; i += 2) {
			dASum[i/2] = arr[i] + arr[i+1];
		}
	}

	if ((num % 2) == 1)
		dASum[dASumSize - 1] = arr[num - 1];

	/*
	for (int i = 0; i < dASumSize; i++) {
		outArr << setw(5) << dASum[i];
	}
	outArr << endl;
	*/

	if (depth != 0) {
		delete[] arr;
	}

	if (dASumSize == 1) {
		sum = dASum[0];
		delete[] dASum;
		return sum;
	}

	depth++;
	return slow_sum(dASum, dASumSize, depth);
}

double fast_sum(double* A, int num) {
	double sum = 0.0;

	omp_set_num_threads(REQUESTED_THREADS);

#pragma omp parallel for reduction(+:sum)
	for (int i = 0; i < num; i++) {
		sum += A[i];
	}
	return sum;
}


int main(int argc, char *argv[]) {
	
	int i_N = 1024;
	int i_R = 1;
	
	ofstream fileFast("outputFast.txt", ios_base::app);
	ofstream fileSlow("outputSlow.txt", ios_base::app);// file2;

	//fileSlow.open("outputSlow.txt", ios_base::app);
	//fileFast.open("outputFast.txt", ios_base::app);
	
	//outArr.open("outArr.txt");

	//Assign R, N, and Threads from arguments
	if (argc == 4) {
		stringstream sstrm;

		sstrm << argv[1];
		sstrm >> i_N;
		sstrm.clear();

		sstrm << argv[2];
		sstrm >> i_R;
		sstrm.clear();

		sstrm << argv[3];
		sstrm >> REQUESTED_THREADS;
	}

	for (; REQUESTED_THREADS >= 1; REQUESTED_THREADS--) {

		//Create input array with dynamic allocation
		double* inputArr = new double[i_N];
		double sum = 0.0;

		//Use these to calculate MFLOPS
		clock_t t;
		double runTime;

		for (int i = 0; i < i_N + 1; i++) {
			inputArr[i] = i;
		}

		t = clock();
		for (int i = 0; i < i_R; i++) {
			sum = slow_sum(inputArr, i_N);
		}

		runTime = 1.0*(clock() - t) / CLOCKS_PER_SEC;

		fileSlow << runTime << "\r\n";

		t = clock();
		for (int i = 0; i < i_R; i++) {
			sum = fast_sum(inputArr, i_N);
		}

		runTime = 1.0*(clock() - t) / CLOCKS_PER_SEC;

		fileFast << runTime << "\r\n";
	}
	fileSlow << "--" << "\r\n";
	fileFast << "--" << "\r\n";

	fileSlow.close();
	fileFast.close();

	//cout.clear();
	//cout << sum << endl;
	
	return 0;
}