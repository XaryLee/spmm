#include <iostream>
#include "csr.cpp"
#include <vector>
using namespace std;

int main(){
    vector<double> nnz{1, 1, 1};
    vector<int> colidx{0, 1, 2};
    vector<int> indptr{0, 1, 2, 3};
    int shape[2]{3, 3};
    SpM test(nnz, colidx, indptr, shape);
    // test = eye(3);
    vector<SpM> mtx_list{test};
    return 0;
}