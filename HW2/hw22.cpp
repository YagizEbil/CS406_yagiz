#include <iostream>
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
    for(int i = 0; i < N; i++) {
        for(int j = 0; j < N; j++) {
            B[i * N + j] = A[i * N + j] + A[j * N + i];
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

    return 0;
}


