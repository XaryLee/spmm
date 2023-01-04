#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include "csr.cpp"
using namespace std;

void gen_mtx_txt(string path, int nnzcount, int shape[], vector<double> nnz, vector<int> colidx,vector<int> rowptr,int xsize){
    string command = "mkdir " + path;
    system(command.c_str());
    ofstream fout;
    
    fout.open(path + "/info.txt");
    fout << to_string(nnzcount) + "\n";
    fout << to_string(shape[0]) + "\n";
    fout << to_string(shape[1]) + "\n";
    fout.close();
    
    fout.open(path + "/col.txt");
    for(auto each:nnz)
        fout << to_string(each) + " ";
    fout.close();

    fout.open(path + "/row.txt");
    for(auto each:nnz)
        fout << to_string(each) + " ";
    fout.close();

    fout.open(path + "/nnz.txt");
    for(auto each:nnz)
        fout << to_string(each) + " ";
    fout.close();

    fout.open(path + "/x.txt");
    for( int i = 0; i < xsize; i++)
        fout << "1 ";
    fout.close();
}

SpM reorder_row(SpM mtx, vector<int> seq){
    vector<double> nnz;
    vector<int> colidx;
    vector<int> rowptr{0};
    for(int s:seq){
        if(mtx.indptr[s] == mtx.indptr[s+1])
            rowptr.push_back(*rowptr.end());
        else{
            rowptr.push_back(*rowptr.end() + (mtx.indptr[s+1] - mtx.indptr[s]));
            for( i = mtx.indptr[s]; i < mtx.indptr[s+1]; i++ ){
                nnz.push_back(mtx.data[i]);
                colidx.push_back(mtx.indices[i])
            }
        }
    }
    return SpM(nnz, colidx, rowptr, mtx.shape);
}

void gen_new_panels(SpM mtx, vector<SpM> &plist, vector<int> &psize_list, int &bnum){}

int main(){
    return 0;
}