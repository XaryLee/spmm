#ifndef _CSR_H_
#define _CSR_H_

#include <iostream>
using namespace std;

class SpM{
public:
    double* data;
    int* indices;
    int* indptr;
    int* shape;
    SpM(double*, int*, int*, int*);
    SpM(){}
    ~SpM();
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
    // for(int i = 0; i < shape[2]; i++)
    //     this->data[i] = data[i];
    // for(int i = 0; i < shape[2]; i++)
    //     this->indices[i] = indices[i];
    // for(int i = 0; i < shape[0]+1; i++)
    //     this->indptr[i] = indptr[i];
    // for(int i = 0; i < 3; i++)
    //     this->shape[i] = shape[i];
    memcpy(this->data, data, sizeof(*data)*shape[2]);
    memcpy(this->indices, indices, sizeof(*indices)*shape[2]);
    memcpy(this->indptr, indptr, sizeof(*indptr)*(shape[0]+1));
    memcpy(this->shape, shape, sizeof(*shape)*3);
}

SpM::~SpM(){
    delete[] data;
    delete[] indices;
    delete[] indptr;
    delete[] shape;
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
        for(int i = 0; i < shape[0]+1; i++)
            cout << indptr[i] << ' ';
        cout << endl;
    }
    cout << "------end------" << endl;
}

#endif