#include <iostream>
#include <vector>
#include <fstream>
#include <string>
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

int main(){
    return 0;
}