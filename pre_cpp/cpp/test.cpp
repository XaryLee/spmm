#include <iostream>
#include "csr.cpp"
#include <vector>
#include <algorithm>
#include <numeric>
using namespace std;

int main(){
    vector<double> nnz{1, 1, 1};
    vector<int> colidx{0, 1, 2};
    vector<int> indptr{0, 1, 2, 3};
    vector<vector<int>> vv;
    vv.push_back(colidx);
    vv.push_back(indptr);
    cout << vv[0][1] << endl;
    // int shape[2]{3, 4};
    // cout << SpM(nnz, colidx, indptr, shape).shape[0] << endl;
    return 0;
}