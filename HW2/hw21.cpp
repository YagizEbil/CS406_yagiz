#include <iostream>
using namespace std;
const int N = 40;
int main(int argc, char** argv) {
    float cube[N][N][N];
    float cube2[N][N][N];
    float matrix[N][N];
    float matrix2[N][N];

    //This part is not modifiable 
    for(int i = 0; i < N; i++) {
        for(int j = 0; j < N; j++) {
            matrix[i][j] = matrix2[i][j] = 0;
            for(int k = 0; k < N; k++) {
                cube[i][j][k] = (i + j) / N;
                cube2[i][j][k] = (j + k) / N;
            } 
        } 
    }

    //TODO: parallelize this  -----------------------------------------------------------------------
    //Report timings with 1,2,4,8 threads - Copy paste your code to the report and explain the details for your parallelization 
    for(int i = 0; i < N; i++) {
        for(int j = 0; j < N; j++) {
            for(int k = 0; k < N; k++) {
                for(int x = 0; x < N; x++) {
                    for(int y = 0; y < N; y++) {
                        for(int z = 0; z < N; z++) {
                            matrix[i][x] += cube[i][x][k] * cube2[z][y][j];
                            matrix2[y][z] += cube[z][k][x] * cube2[y][j][i];
                        }
                    }
                }
            }
        }
    }
    
    //This is a basic checksum for the output
    float sum1 = 0, sum2 = 0;
    for(int i = 0; i < N; i++) {
        for(int j = 0; j < N; j++) {
            sum1 += matrix[i][j];
            sum2 += matrix2[i][j];
        }
    }
    cout << sum1 << " " << sum2 << endl;

    return 0;
}


