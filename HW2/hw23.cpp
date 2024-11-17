#include <iostream>
#include <cstdlib>
#include <immintrin.h> // For SIMD
#include <omp.h>
using namespace std;

const int N = 2000;

int main(int argc, char** argv) {
    float* A;
    float* B;
    posix_memalign((void**)&A, 64, N * N * sizeof(float));
    posix_memalign((void**)&B, 64, N * N * sizeof(float));

    //This part is not modifiable 
    for(int i = 0; i < N; i++) {
        for(int j = 0; j < N; j++) {
            B[i * N + j] = 0;
            A[i * N + j] = float(i + j) / N;
        } 
    }

    // Parallelize B = A * A with SIMD
    #pragma omp parallel for collapse(2) schedule(guided)
    for(int i = 0; i < N; i++) {
        for(int j = 0; j < N; j++) {
            __m256 sum = _mm256_setzero_ps();
            for(int k = 0; k < N; k += 8) {
                __m256 a = _mm256_load_ps(&A[i * N + k]);
                __m256 b = _mm256_load_ps(&A[j * N + k]);
                sum = _mm256_fmadd_ps(a, b, sum);
            }
            float temp[8];
            _mm256_store_ps(temp, sum);
            B[i * N + j] = temp[0] + temp[1] + temp[2] + temp[3] + temp[4] + temp[5] + temp[6] + temp[7];
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

    free(A);
    free(B);

    return 0;
}


