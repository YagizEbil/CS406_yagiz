#include <iostream>
#include <fstream>
#include <chrono>
using namespace std;

int main() {
    //Reading matrix from a binary file --------------------------------------------
    int n = -1; //number of rows/columns
    int nnz = -1; //number of nonzeros in a matrix;

    ifstream infile("/data/sparse_matrix.bin", ios::binary);

    // Read size n + 1 for row_ptrs
    infile.read(reinterpret_cast<char*>(&n), sizeof(int));
    int* row_ptrs = new int[n + 1];
    infile.read(reinterpret_cast<char*>(row_ptrs), (n + 1) * sizeof(int));

    // Read size nnz for col_ids and values
    infile.read(reinterpret_cast<char*>(&nnz), sizeof(int));
    int* col_ids = new int[nnz];
    double* values = new double[nnz];
    infile.read(reinterpret_cast<char*>(col_ids), nnz * sizeof(int));
    infile.read(reinterpret_cast<char*>(values), nnz * sizeof(double));

    infile.close();
    cout << "Matrix is read from binary file:" << endl;
    cout << "\tNumber of rows/columns: " << n << endl;
    cout << "\tNumber of nonzeros: " << nnz << endl << endl;

    //Initial SpMV iterations -------------------------------------------------------
    int num_iterations = 10;

    double* y = new double[n]; // Output vector
    double* x = new double[n]; // First input vector

    for (int i = 0; i < n; i++) { // Do not modify the initialization
        x[i] = 0.00001f;  
        y[i] = 0;
    }
  
    auto start_time = chrono::high_resolution_clock::now();
    for (int iter = 0; iter < num_iterations; iter++) {
        for (int i = 0; i < n; i++) {
            double dotproduct = 0;
            for (int idx = row_ptrs[i]; idx < row_ptrs[i + 1]; idx++) {
                dotproduct += values[idx] * x[col_ids[idx]];
            }
            y[i] = dotproduct;
        }
        // Use y as the next x
        double* temp = x;
        x = y;
        y = temp;
    }
    auto end_time = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = end_time - start_time;

    cout << "Time taken for the original SpMV over " << num_iterations << " iterations: " << elapsed.count() << " seconds" << endl;

    double avg = x[0] / n, maxv = x[0], minv = x[0];
    for(int i = 1; i < n; i++) {
      avg += x[i] / n;
      maxv = max(maxv, x[i]);
      minv = min(minv, x[i]);
    }

    //Only for correctness - will be used to check correctness
    cout << "Original statistics: " << avg << " " << maxv << " " << minv << endl << endl;

    //your cannot modify before this line 
    // Do your optimizations here -----------------------------------------------------------------------------
    //......// 
    //......//

    //---------------------------------------------------------------------------------------------------------
    //you cannot modify after this line

    // Repeat the SpMV - this is the same code above (do not change) ------------------------------------------
    for (int i = 0; i < n; i++) { // Do not modify the initialization
        x[i] = 0.00001f;  
        y[i] = 0;
    }
  
    start_time = chrono::high_resolution_clock::now();
    for (int iter = 0; iter < num_iterations; iter++) {
        for (int i = 0; i < n; i++) {
            double dotproduct = 0;
            for (int idx = row_ptrs[i]; idx < row_ptrs[i + 1]; idx++) {
                dotproduct += values[idx] * x[col_ids[idx]];
            }
            y[i] = dotproduct;
        }
        // Use y as the next x
        double* temp = x;
        x = y;
        y = temp;
    }
    end_time = chrono::high_resolution_clock::now();
    elapsed = end_time - start_time;
    cout << "Time taken for the modified SpMV over " << num_iterations << " iterations: " << elapsed.count() << " seconds" << endl;

    avg = x[0] / n, maxv = x[0], minv = x[0];
    for(int i = 1; i < n; i++) {
      avg += x[i] / n;
      maxv = max(maxv, x[i]);
      minv = min(minv, x[i]);
    }

    //These must be the same with the above
    cout << "With optimization statistics: " << avg << " " << maxv << " " << minv << endl;

    return 0;
}