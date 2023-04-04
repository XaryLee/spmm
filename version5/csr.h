#ifndef _CSR_H_
#define _CSR_H_

#include <iostream>
#include <cstring>
using namespace std;

class SpM{
public:
    double* data;
    int* indices;
    int* indptr;
    int* shape;
    SpM(double*, int*, int*, int*);
    SpM(int*);
    SpM();
    ~SpM();
    void operator=(const SpM &);
    void check(bool);
};

SpM::SpM(double* data, int* indices, int* indptr, int* shape){
    // shape:[nrow ncol nnz]
    // this->data = data;
    // this->indices = indices;
    // this->indptr = indptr;
    // this->shape = shape;
    this->data = new double[shape[2]];
    this->indices = new int[shape[2]];
    this->indptr = new int[shape[0]+1];
    this->shape = new int[3];

    int nnz = shape[2];
    int row = shape[0] + 1;

    #pragma omp parallel for num_threads(CORENUM)
    for(int i = 0; i < nnz; i++){
        this->data[i] = data[i];
        this->indices[i] = indices[i];
        if(i<row) this->indptr[i] = indptr[i];
    }
    // #pragma omp parallel for num_threads(CORENUM)
    // for(int i = 0; i < shape[0]+1; i++)
    //     this->indptr[i] = indptr[i];
    for(int i = 0; i < 3; i++)
        this->shape[i] = shape[i];
    // memcpy(this->data, data, sizeof(double)*shape[2]);
    // memcpy(this->indices, indices, sizeof(int)*shape[2]);
    // memcpy(this->indptr, indptr, sizeof(int)*(shape[0]+1));
    // memcpy(this->shape, shape, sizeof(int)*3);
}

SpM::SpM(){
    data = NULL;
    indices = NULL;
    indptr = NULL;
    shape = NULL;
}

SpM::SpM(int* shape){
    this->shape = new int[3];
    memcpy(this->shape, shape, sizeof(int)*3);
    data = new double[shape[2]]();
    indices = new int[shape[2]]();
    indptr = new int[shape[0]+1]();
}

void SpM::operator=(const SpM &mr){
    // cout << "=" << endl;
    delete[] data, indices, indptr, shape;
    // cout << "delete done" << endl;
    data = new double[mr.shape[2]];
    indices = new int[mr.shape[2]];
    indptr = new int[mr.shape[0]+1];
    shape = new int[3];
    // cout << "new done" << endl;
    memcpy(shape, mr.shape, sizeof(int)*3);
    memcpy(data, mr.data, sizeof(double)*shape[2]);
    // cout << "0" << endl;
    memcpy(indices, mr.indices, sizeof(int)*shape[2]);
    // cout << "0" << endl;
    memcpy(indptr, mr.indptr, sizeof(int)*(shape[0]+1));

    // cout << "done" << endl;
}

SpM::~SpM(){
    delete[] data; data = NULL;
    delete[] indices; indices = NULL;
    delete[] indptr; indptr = NULL;
    delete[] shape; shape = NULL;
}

void SpM::check(bool show_indptr=true){
    // cout << "data_size: " << data.size() << endl;
    // cout << "indices_size: " << indices.size() << endl;
    // cout << "indptr_size: " << indptr.size() << endl;
    // if(show_indptr){
    //     cout << "indptr:\n";
    //     for(auto each:indptr)
    //         cout << each << ' ';
    //     cout << endl;
    // }
    cout << "shape: ";
    for(int i = 0; i < 3; i++)
        cout << shape[i] << ' ';
    cout << endl;
    if(show_indptr){
        cout << "indptr:\n";
        for(int i = 0; i < min(shape[0]+1, 100); i++){
            cout << indptr[i] << ' ';
        }
        cout << endl;
    }
    cout << "------end------" << endl;
}

#endif
