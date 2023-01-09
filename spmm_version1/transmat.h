#ifndef _TRANSMAT_H_
#define _TRANSMAT_H_

#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include "csr.h"

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
            for( int i = mtx.indptr[s]; i < mtx.indptr[s+1]; i++ ){
                nnz.push_back(mtx.data[i]);
                colidx.push_back(mtx.indices[i]);
            }
        }
    }
    return SpM(nnz, colidx, rowptr, mtx.shape);
}

void gen_new_panels(SpM mtx, vector<SpM> &plist, vector<int> &psize_list, int &bnum){
    // plist = []
    // psize_list = [0]
    //!!!!!
    // for(auto i:mtx.indptr) cout<<i<<" ";
    int threshold = 512 * 1024 / 8;
    // int threshold = 4;
    int counter = 0;
    vector<bool> element_array(mtx.shape[1],0);
    vector<int>::const_iterator beg, end;
    vector<double>::const_iterator dbeg, dend;
    for(int index = 0; index < (int)(mtx.indptr.size() - 1); index++ ){
        beg = mtx.indices.begin() + mtx.indptr[index];
        end = mtx.indices.begin() + mtx.indptr[index + 1];
        cout<<mtx.indptr[index]<< " "<<mtx.indptr[index+1]<<" ";
        vector<int> row_indices(beg, end);
        for(int value:row_indices){
            cout<<value<<" ";
            if(element_array[value-1] == 0)
                counter++;
            element_array[value-1] = 1;
            cout<<value<<" ";
        }
        if(counter >= threshold){
            element_array.resize(mtx.shape[1],0);
            counter = 0;
            psize_list.push_back(index + 1);
            // p_indptr = mtx.indptr[psize_list[-2]:psize_list[-1]+1]
            beg = mtx.indptr.begin() + *(psize_list.end() - 2);
            end = mtx.indptr.begin() + *(psize_list.end() - 1);
            vector<int> p_indptr(beg, end);
            // p_nnz = mtx.data[p_indptr[0]:p_indptr[-1]]
            dbeg = mtx.data.begin() + *(p_indptr.begin());
            dend = mtx.data.begin() + *(p_indptr.end() - 1);
            vector<double> p_nnz(dbeg, dend);
            // p_indices = mtx.indices[p_indptr[0]:p_indptr[-1]]
            beg = mtx.indices.begin() + *(p_indptr.begin());
            end = mtx.indices.begin() + *(p_indptr.end() - 1);
            vector<int> p_indices(beg, end); 
            int offset = p_indptr[0];
            // p_indptr = p_indptr - offset
            for(int each:p_indptr)
                each -= offset;
            // pm = scipy.sparse.csr_matrix((p_nnz, p_indices, p_indptr), shape=(psize_list[-1]-psize_list[-2], mtx.shape[1]))
            vector<int> pm_shape(3);
            pm_shape[0] = *(psize_list.end() - 1) - *(psize_list.end() - 2);
            pm_shape[1] = mtx.shape[1];
            pm_shape[2] = mtx.shape[2];
            SpM pm(p_nnz, p_indices, p_indptr, pm_shape);
            plist.push_back(pm);
        }
    }
    psize_list.push_back(mtx.indptr.size() - 1);
    // p_indptr = mtx.indptr[psize_list[-2]:psize_list[-1]+1]
    beg = mtx.indptr.begin() + *(psize_list.end() - 2);
    end = mtx.indptr.begin() + *(psize_list.end() - 1);
    vector<int> p_indptr(beg, end);
    // p_nnz = mtx.data[p_indptr[0]:p_indptr[-1]]
    dbeg = mtx.data.begin() + *(p_indptr.begin());
    dend = mtx.data.begin() + *(p_indptr.end() - 1);
    vector<double> p_nnz(dbeg, dend);
    // p_indices = mtx.indices[p_indptr[0]:p_indptr[-1]]
    beg = mtx.indices.begin() + *(p_indptr.begin());
    end = mtx.indices.begin() + *(p_indptr.end() - 1);
    vector<int> p_indices(beg, end); 
    int offset = p_indptr[0];
    // p_indptr = p_indptr - offset
    for(int each:p_indptr)
        each -= offset;
    // pm = scipy.sparse.csr_matrix((p_nnz, p_indices, p_indptr), shape=(psize_list[-1]-psize_list[-2], mtx.shape[1]))
    vector<int> pm_shape(3);
    pm_shape[0] = *(psize_list.end() - 1) - *(psize_list.end() - 2);
    pm_shape[1] = mtx.shape[1];
    pm_shape[2] = mtx.shape[2];
    SpM pm(p_nnz, p_indices, p_indptr, pm_shape);
    plist.push_back(pm);
    bnum = psize_list.size() - 1;
    cout << "the number of blocks is" << bnum << endl;

}

// int main(){
//     return 0;
// }
#endif