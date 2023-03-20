#ifndef _TRANSMAT_H_
#define _TRANSMAT_H_
#include <iostream>
#include <vector>
#include "csr.h"

using namespace std;

SpM reorder_row(SpM &mtx, int* seq){
    // len(seq) = mtx.shape[0]
    double* nnz = new double[mtx.shape[2]];
    int* colidx = new int[mtx.shape[2]];
    int* rowptr = new int[mtx.shape[0]+1]();
    int tail = 0;
    for(int i = 0; i < mtx.shape[0]; i++){
        int s = seq[i];
        if(mtx.indptr[s] == mtx.indptr[s+1])
            rowptr[i+1] = rowptr[i];
        else{
            rowptr[i+1] = (rowptr[i] + (mtx.indptr[s+1] - mtx.indptr[s]));
            for( int j = mtx.indptr[s]; j < mtx.indptr[s+1]; j++ ){
                nnz[tail] = (mtx.data[j]);
                colidx[tail] = (mtx.indices[j]);
                tail++;
                // 这里可以用memcpy加速
            }
        }
    }
    SpM mr(nnz, colidx, rowptr, mtx.shape);
    delete[] nnz, colidx, rowptr;
    return mr;
}

void gen_new_panels(SpM &mtx, vector<SpM> &plist, vector<int> &psize_list, int &bnum){
    // plist = []
    // psize_list = [0]
    //!!!!!
    // for(auto i:mtx.indptr) cout<<i<<" ";
    int threshold = 512 * 1024 / 8;
    // int threshold = 4;
    int counter = 0;
    vector<bool> element_array(mtx.shape[1],0);
    // vector<int>::const_iterator beg, end;
    int *beg, *end;
    double *dbeg, *dend;
    for(int index = 0; index < (mtx.shape[0]); index++ ){
        beg = mtx.indices + mtx.indptr[index];
        end = mtx.indices + mtx.indptr[index + 1];
        // for(auto i : mtx.indptr) cout<<i<<" ";
        // cout<<"index"<<index<<endl;
        // cout<<"indptr"<<mtx.indptr[index]<< " "<<mtx.indptr[index+1]<<" "<<endl;
        vector<int> row_indices(beg, end);
        for(int value:row_indices){
            // cout<<value<<" ";
            if(element_array[value] == 0)
                counter++;
            element_array[value] = 1;
            // cout<<value<<" ";
        }
        if(counter >= threshold){
            element_array.resize(mtx.shape[1],0);
            counter = 0;
            psize_list.push_back(index + 1);
            // p_indptr = mtx.indptr[psize_list[-2]:psize_list[-1]+1]
            beg = mtx.indptr + *(psize_list.end() - 2);
            end = mtx.indptr + *(psize_list.end() - 1) + 1;
            auto p_indptr_len = end - beg;
            int* p_indptr = new int[p_indptr_len];
            memcpy(p_indptr, beg, sizeof(int)*(end-beg));
            // vector<int> p_indptr(beg, end);
            // p_nnz = mtx.data[p_indptr[0]:p_indptr[-1]]
            dbeg = mtx.data + p_indptr[0];
            dend = mtx.data + p_indptr[p_indptr_len-1];
            double* p_nnz = new double[dbeg - dend];
            memcpy(p_nnz, dbeg, sizeof(double)*(dbeg-dend));
            // vector<double> p_nnz(dbeg, dend);
            // p_indices = mtx.indices[p_indptr[0]:p_indptr[-1]]
            beg = mtx.indices + p_indptr[0];
            end = mtx.indices + p_indptr[p_indptr_len-1];
            int* p_indices = new int[end-beg];
            memcpy(p_indices, beg, sizeof(int)*(end-beg));
            // vector<int> p_indices(beg, end); 
            int offset = p_indptr[0];
            // p_indptr = p_indptr - offset
            for(int i = 0; i < p_indptr_len; i++)
                p_indptr[i] -= offset;
            // pm = scipy.sparse.csr_matrix((p_nnz, p_indices, p_indptr), shape=(psize_list[-1]-psize_list[-2], mtx.shape[1]))
            int pm_shape[3];
            pm_shape[0] = *(psize_list.end() - 1) - *(psize_list.end() - 2);
            pm_shape[1] = mtx.shape[1];
            pm_shape[2] = mtx.shape[2];
            SpM pm(p_nnz, p_indices, p_indptr, pm_shape);
            delete[] p_nnz, p_indices, p_indptr;
            plist.push_back(pm);
        }
    }
    psize_list.push_back(mtx.shape[0]);
    // p_indptr = mtx.indptr[psize_list[-2]:psize_list[-1]+1]
    beg = mtx.indptr + *(psize_list.end() - 2);
    end = mtx.indptr + *(psize_list.end() - 1) + 1;
    int p_indptr_len = end - beg;
    int* p_indptr = new int[p_indptr_len];
    memcpy(p_indptr, beg, sizeof(int)*(end-beg));
    // vector<int> p_indptr(beg, end);
    // p_nnz = mtx.data[p_indptr[0]:p_indptr[-1]]
    dbeg = mtx.data + p_indptr[0];
    dend = mtx.data + p_indptr[p_indptr_len-1];
    double* p_nnz = new double[dbeg - dend];
    memcpy(p_nnz, dbeg, sizeof(double)*(dbeg-dend));
    // vector<double> p_nnz(dbeg, dend);
    // p_indices = mtx.indices[p_indptr[0]:p_indptr[-1]]
    beg = mtx.indices + p_indptr[0];
    end = mtx.indices + p_indptr[p_indptr_len-1];
    int* p_indices = new int[end-beg];
    memcpy(p_indices, beg, sizeof(int)*(end-beg));
    // vector<int> p_indices(beg, end); 
    int offset = p_indptr[0];
    // p_indptr = p_indptr - offset
    for(int i = 0; i < p_indptr_len; i++)
        p_indptr[i] -= offset;
    // pm = scipy.sparse.csr_matrix((p_nnz, p_indices, p_indptr), shape=(psize_list[-1]-psize_list[-2], mtx.shape[1]))
    // vector<int> pm_shape(3);
    int pm_shape[3];
    pm_shape[0] = *(psize_list.end() - 1) - *(psize_list.end() - 2);
    pm_shape[1] = mtx.shape[1];
    pm_shape[2] = mtx.shape[2];
    SpM pm(p_nnz, p_indices, p_indptr, pm_shape);
    delete[] p_nnz, p_indices, p_indptr;
    plist.push_back(pm);
    bnum = psize_list.size() - 1;
    cout << "the number of blocks is" << bnum << endl;

}

#endif
