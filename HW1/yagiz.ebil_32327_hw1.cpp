#include <iostream>
#include <fstream>
#include <chrono>
#include <limits.h>
#include <vector>
#include <queue>
#include <unordered_map>
#include <algorithm>
#include <set>
#include <functional>
#include <numeric>

using namespace std;

// Kadir Yagiz Ebil - 32327

void nested_dissection(int left, int right, vector<int>& ordering) {
    if (left >= right) return;

    int mid = (left + right) / 2;

    nested_dissection(left, mid, ordering);
    nested_dissection(mid + 1, right, ordering);

    ordering.push_back(mid);
}

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
#define RCM

#ifdef RCM
    unordered_map<int, vector<int>> graph;

    for (int i = 0; i < n; ++i) {
        for (int idx = row_ptrs[i]; idx < row_ptrs[i + 1]; ++idx) {
            int col = col_ids[idx];
            graph[i].push_back(col);  
            graph[col].push_back(i);  // Symmetric for undirected graph
        }
    }

    // Function to find the node with the minimum degree
    int min_degree = INT_MAX;
    int start_node = -1;
    for (const auto& entry : graph) {
        if (entry.second.size() < min_degree) {
            min_degree = entry.second.size();
            start_node = entry.first;
        }
    }

    vector<bool> visited(n, false);
    vector<int> ordering;
    queue<int> q;
    if (start_node != -1) {
        visited[start_node] = true;
    } else {
        cerr << "Error: No valid start node found." << endl;
        return -1;
    }
    q.push(start_node);

    while (!q.empty()) {
        int node = q.front();
        q.pop();
        ordering.push_back(node);

        vector<int> neighbors = graph[node];
        sort(neighbors.begin(), neighbors.end(), [&](int a, int b) {
            return graph[a].size() < graph[b].size();
        });

        for (int neighbor : neighbors) {
            if (!visited[neighbor]) {
                visited[neighbor] = true;
                q.push(neighbor);
            }
        }
    }

    reverse(ordering.begin(), ordering.end());

    unordered_map<int, int> old_to_new;
    for (int i = 0; i < n; ++i) {
        old_to_new[ordering[i]] = i;
    }

    int* new_row_ptrs_rcm = new int[n + 1];
    vector<int> new_col_ids_rcm(nnz);
    vector<double> new_values_rcm(nnz);

    new_row_ptrs_rcm[0] = 0;
    int new_nnz_rcm = 0;
    for (int i = 0; i < n; ++i) {
        int old_row = ordering[i];
        for (int idx = row_ptrs[old_row]; idx < row_ptrs[old_row + 1]; ++idx) {
            int old_col = col_ids[idx];
            int new_col = old_to_new[old_col];
            new_col_ids_rcm[new_nnz_rcm] = new_col;
            new_values_rcm[new_nnz_rcm] = values[idx];
            ++new_nnz_rcm;
        }
        new_row_ptrs_rcm[i + 1] = new_nnz_rcm;
    }

    std::copy(new_row_ptrs_rcm, new_row_ptrs_rcm + n + 1, row_ptrs);
    std::copy(new_col_ids_rcm.begin(), new_col_ids_rcm.end(), col_ids);
    std::copy(new_values_rcm.begin(), new_values_rcm.end(), values);

    delete[] new_row_ptrs_rcm;
    cout << "Applied RCM ordering to the matrix." << endl;
#endif

#ifdef RCM2
    unordered_map<int, vector<int>> graph;

    for (int i = 0; i < n; ++i) {
        for (int idx = row_ptrs[i]; idx < row_ptrs[i + 1]; ++idx) {
            int col = col_ids[idx];
            graph[i].push_back(col);
            graph[col].push_back(i); 
        }
    }

    int min_degree = INT_MAX;
    int start_node = -1;
    for (const auto& entry : graph) {
        if (entry.second.size() < min_degree) {
            min_degree = entry.second.size();
            start_node = entry.first;
        }
    }

    // Priority Queue for degree-based ordering (min-heap)
    vector<bool> visited(n, false);
    vector<int> ordering;
    priority_queue<pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>>> pq;

    if (start_node != -1) {
        pq.push({min_degree, start_node});
    } else {
        cerr << "Error: No valid start node found." << endl;
        return -1;
    }

    while (!pq.empty()) {
        int node = pq.top().second;
        pq.pop();

        if (visited[node]) continue; 
        visited[node] = true;
        ordering.push_back(node);

        // Push unvisited neighbors to the priority queue
        for (int neighbor : graph[node]) {
            if (!visited[neighbor]) {
                pq.push({graph[neighbor].size(), neighbor});  // Neighbor degree as priority
            }
        }
    }

    reverse(ordering.begin(), ordering.end());

    unordered_map<int, int> old_to_new;
    for (int i = 0; i < n; ++i) {
        old_to_new[ordering[i]] = i;
    }

    int* new_row_ptrs_rcm = new int[n + 1];
    vector<int> new_col_ids_rcm(nnz);
    vector<double> new_values_rcm(nnz);

    new_row_ptrs_rcm[0] = 0;
    int new_nnz_rcm = 0;
    for (int i = 0; i < n; ++i) {
        int old_row = ordering[i];
        for (int idx = row_ptrs[old_row]; idx < row_ptrs[old_row + 1]; ++idx) {
            int old_col = col_ids[idx];
            int new_col = old_to_new[old_col];
            new_col_ids_rcm[new_nnz_rcm] = new_col;
            new_values_rcm[new_nnz_rcm] = values[idx];
            ++new_nnz_rcm;
        }
        new_row_ptrs_rcm[i + 1] = new_nnz_rcm;
    }

    std::copy(new_row_ptrs_rcm, new_row_ptrs_rcm + n + 1, row_ptrs);
    std::copy(new_col_ids_rcm.begin(), new_col_ids_rcm.end(), col_ids);
    std::copy(new_values_rcm.begin(), new_values_rcm.end(), values);

    delete[] new_row_ptrs_rcm;
    cout << "Applied RCM-2 ordering with priority queue to the matrix." << endl;
#endif

#ifdef BLOCK
    int block_size = 64; 

    vector<int> block_indices(n / block_size + (n % block_size != 0));
    for (int i = 0; i < block_indices.size(); ++i) {
        block_indices[i] = i * block_size;
    }

    vector<int> new_order(n);
    int new_index = 0;
    for (int block_start : block_indices) {
        for (int i = 0; i < block_size; ++i) {
            int row = block_start + i;
            if (row < n) {
                new_order[row] = new_index++;
            }
        }
    }

    int* new_row_ptrs = new int[n + 1];
    vector<int> new_col_ids(row_ptrs[n]); 
    vector<double> new_values(row_ptrs[n]); 

    new_row_ptrs[0] = 0;
    int new_nnz = 0;
    for (int i = 0; i < n; ++i) {
        for (int idx = row_ptrs[i]; idx < row_ptrs[i + 1]; ++idx) {
            int old_col = col_ids[idx];
            int new_row = new_order[i];
            int new_col = new_order[old_col]; 
            new_col_ids[new_nnz] = new_col;
            new_values[new_nnz] = values[idx];
            ++new_nnz;
        }
        new_row_ptrs[new_order[i] + 1] = new_nnz;
    }

    for (int i = 1; i < n + 1; ++i) {
        if (new_row_ptrs[i] == 0) {
            new_row_ptrs[i] = new_row_ptrs[i - 1];
        }
    }

    std::copy(new_row_ptrs, new_row_ptrs + n + 1, row_ptrs);
    std::copy(new_col_ids.begin(), new_col_ids.end(), col_ids);
    std::copy(new_values.begin(), new_values.end(), values);

    delete[] new_row_ptrs;
    cout << "Applied Block ordering to the matrix." << endl;
#endif

#ifdef ND
    //Not complete 
    vector<int> nd_ordering(n);
    iota(nd_ordering.begin(), nd_ordering.end(), 0);

    nested_dissection(0, n - 1, nd_ordering);

    int* new_row_ptrs_nd = new int[n + 1];
    vector<int> new_col_ids_nd(nnz);
    vector<double> new_values_nd(nnz);

    new_row_ptrs_nd[0] = 0;
    int new_nnz_nd = 0;
    for (int i = 0; i < n; ++i) {
        int old_row = nd_ordering[i];
        for (int idx = row_ptrs[old_row]; idx < row_ptrs[old_row + 1]; ++idx) {
            int old_col = col_ids[idx];
            int new_col = nd_ordering[old_col];
            new_col_ids_nd[new_nnz_nd] = new_col;
            new_values_nd[new_nnz_nd] = values[idx];
            ++new_nnz_nd;
        }
        new_row_ptrs_nd[i + 1] = new_nnz_nd;
    }

    delete[] row_ptrs;
    row_ptrs = new_row_ptrs_nd;
    std::copy(new_col_ids_nd.begin(), new_col_ids_nd.end(), col_ids);
    std::copy(new_values_nd.begin(), new_values_nd.end(), values);

    cout << "Applied Nested Dissection (ND) ordering to the matrix." << endl;
#endif

#ifdef RABBIT
    //Not complete
    unordered_map<int, vector<int>> graph;

    for (int i = 0; i < n; ++i) {
        for (int idx = row_ptrs[i]; idx < row_ptrs[i + 1]; ++idx) {
            int col = col_ids[idx];
            graph[i].push_back(col);  
        }
    }

    vector<pair<int, int>> degrees(n);  
    for (int i = 0; i < n; ++i) {
        degrees[i] = {graph[i].size(), i};
    }

    sort(degrees.begin(), degrees.end());

    unordered_map<int, int> new_index_map;
    for (int new_index = 0; new_index < n; ++new_index) {
        new_index_map[degrees[new_index].second] = new_index;
    }

    int* new_row_ptrs = new int[n + 1];
    vector<int> new_col_ids(nnz);
    vector<double> new_values(nnz);

    new_row_ptrs[0] = 0;
    int new_nnz = 0;

    for (int old_row = 0; old_row < n; ++old_row) {
        for (int idx = row_ptrs[old_row]; idx < row_ptrs[old_row + 1]; ++idx) {
            int old_col = col_ids[idx];
            int new_col = new_index_map[old_col]; 
            new_col_ids[new_nnz] = new_col;
            new_values[new_nnz] = values[idx];
            ++new_nnz;
        }
        new_row_ptrs[new_index_map[old_row] + 1] = new_nnz; 
    }

    for (int i = 1; i < n + 1; ++i) {
        if (new_row_ptrs[i] == 0) {
            new_row_ptrs[i] = new_row_ptrs[i - 1];
        }
    }

    std::copy(new_row_ptrs, new_row_ptrs + n + 1, row_ptrs);
    std::copy(new_col_ids.begin(), new_col_ids.end(), col_ids);
    std::copy(new_values.begin(), new_values.end(), values);

    delete[] new_row_ptrs;
    cout << "Applied Rabbit ordering to the matrix." << endl;
#endif

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

//Memory leak maybe? If modification is not allowed please delete.
    delete[] row_ptrs;
    delete[] col_ids;
    delete[] values;
    delete[] x;
    delete[] y;

    return 0;
}