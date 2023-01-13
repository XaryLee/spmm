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
    vector<int> shape;
    SpM(vector <double>, vector<int>, vector<int>, vector<int>);
    SpM(){}
    void check(bool);
};

SpM::SpM(vector <double> data, vector<int> indices, vector<int> indptr, vector<int> shape){
    this->data = data;
    this->indices = indices;
    this->indptr = indptr;
    this->shape = shape;
}

void SpM::check(bool show_indptr=true){
    cout << "data_size: " << data.size() << endl;
    cout << "indices_size: " << indices.size() << endl;
    cout << "indptr_size: " << indptr.size() << endl;
    if(show_indptr){
        cout << "indptr:\n";
        for(auto each:indptr)
            cout << each << ' ';
        cout << endl;
    }
    cout << "------end------" << endl;
}

#endif
