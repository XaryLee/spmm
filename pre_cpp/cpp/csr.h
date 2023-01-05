#ifndef _CSR_H_
#define _CSR_H_

#include <iostream>
#include <vector>
using namespace std;

class SpM{
public:
    vector<double> data;
    vector<int> indices;
    vector<int> indptr;
    int* shape;
    SpM(vector <double>, vector<int>, vector<int>, int*);
};

SpM::SpM(vector <double> data, vector<int> indices, vector<int> indptr, int* shape){
    this->data = data;
    this->indices = indices;
    this->indptr = indptr;
    this->shape = shape;
}

#endif
