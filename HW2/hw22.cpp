#include <iostream>
#include <omp.h>
using namespace std;

const int N = 10000;

int main(int argc, char** argv) {
    float* A = new float[N * N];
    float* B = new float[N * N];

    //This part is not modifiable 
    for(int i = 0; i < N; i++) {
        for(int j = 0; j < N; j++) {
            B[i * N + j] = 0;
            A[i * N + j] = (i + j) / N;
        } 
    }

    //TODO: parallelize B = A + Atranspose  -----------------------------------------------------------------------
    //Report timings with 1,2,4,8 threads - Copy paste your code to the report and explain the details 
    #pragma omp parallel for collapse(2) schedule(static)
    for(int i = 0; i < N; i++) {
        for(int j = 0; j < N; j += 4) {
            B[i * N + j] = A[i * N + j] + A[j * N + i];
            if (j + 1 < N) B[i * N + j + 1] = A[i * N + j + 1] + A[j * N + i + 1];
            if (j + 2 < N) B[i * N + j + 2] = A[i * N + j + 2] + A[j * N + i + 2];
            if (j + 3 < N) B[i * N + j + 3] = A[i * N + j + 3] + A[j * N + i + 3];
        }
    }
    
    //This is a basic checksum for the output
    float sum1 = 0, sum2 = 0;
    for(int i = 0; i < N; i++) {
        for(int j = 0; j < N; j++) {
            sum1 += A[i * N + j];
            sum2 += B[i * N + j];
        }
    }
    cout << sum1 << " " << sum2 << endl;

    //memory leak?
    delete[] A;
    delete[] B;

    return 0;
}


