#include <iostream>
using namespace std;

const int N = 2000;

int main(int argc, char** argv) {
    float* A = new float[N * N];
    float* B = new float[N * N];

    //This part is not modifiable 
    for(int i = 0; i < N; i++) {
        for(int j = 0; j < N; j++) {
            B[i * N + j] = 0;
            A[i * N + j] = float(i + j) / N;
        } 
    }

    //TODO: parallelize B = A * A  -----------------------------------------------------------------------
    //  Copy paste your code to the report and explain the details for your parallelization 
    //  Report timings with 1,2,4,8 threads  
    for(int i = 0; i < N; i++) {
        for(int j = 0; j < N; j++) {
            for(int k = 0; k < N; k++) {
                B[i * N + j] += A[i * N + k] + A[k * N + j];
            }
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


