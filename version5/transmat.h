#ifndef _TRANSMAT_H_
#define _TRANSMAT_H_
#include <iostream>
#include <vector>
#include <cstring>
#include <ctime>
#include "csr.h"

using namespace std;

// SpM reorder_row(SpM &mtx, int* seq){
    
//     SpM mr(mtx.shape);
//     // mr.check(false);
//     int shape = mtx.shape[0];
//     int tail = 0;
//     for(int i = 0; i < shape; i++){
//             int s = seq[i];
//             mr.indptr[i+1] = (mr.indptr[i] + (mtx.indptr[s+1] - mtx.indptr[s]));
//             for( int j = mtx.indptr[s]; j < mtx.indptr[s+1]; j++ ){
//                 mr.data[tail] = mtx.data[j];
//                 mr.indices[tail] = mtx.indices[j];
//                 tail++;
//             }
//         }
//     return mr;
// }

SpM reorder_row(SpM &mtx, int* seq){
    
    SpM mr(mtx.shape);
    int shape = mtx.shape[0];
    vector<int> tail(CORENUM,0);
    int MIN_CHUNK_SIZE = 20000;
    if (shape > MIN_CHUNK_SIZE){
        for(int c=1;c<CORENUM;c++){
            tail[c] += tail[c-1];
            int sup = floor(shape * (c-1) / CORENUM);
            int up = floor(shape * (c) / CORENUM);
            for(int i = sup ; i < up; i++) {
                tail[c]+=mtx.indptr[seq[i]+1] - mtx.indptr[seq[i]];
            }
            mr.indptr[up] = tail[c];
        }
        #pragma omp parallel for num_threads(CORENUM)
        for(int c =0;c<CORENUM;c++)
        {
            for(int i = floor(shape * (c) / CORENUM) ; i < floor(shape * (c+1) / CORENUM); i++){
                int s = seq[i];
                //竞争关系
                mr.indptr[i+1] = (mr.indptr[i] + (mtx.indptr[s+1] - mtx.indptr[s]));
                for( int j = mtx.indptr[s]; j < mtx.indptr[s+1]; j++ ){
                    mr.data[tail[c]] = mtx.data[j];
                    mr.indices[tail[c]] = mtx.indices[j];
                    tail[c]++;
                }
            }
        }
    }
    else{
        int tail = 0;
        for(int i = 0; i < shape; i++){
            int s = seq[i];
            mr.indptr[i+1] = (mr.indptr[i] + (mtx.indptr[s+1] - mtx.indptr[s]));
            for( int j = mtx.indptr[s]; j < mtx.indptr[s+1]; j++ ){
                mr.data[tail] = mtx.data[j];
                mr.indices[tail] = mtx.indices[j];
                tail++;
            }
        }
    }
    return mr;
}

int gen_new_panels(SpM &mtx, SpM* &plist, int* &psize_list,int &end_psize_list, int &bnum){

    // 返回值是plist的长度

    // plist = []
    // psize_list = [0]
    //!!!!!
    // for(auto i:mtx.indptr) cout<<i<<" ";
    // cout << "enter gen_new_panels" << endl;
    int tail = 0;
    int threshold = 512 * 1024 / 8; //65536
    // int threshold = 4;
    int counter = 0;
    bool* element_array = new bool[mtx.shape[1]]();
    // vector<int>::const_iterator beg, end;
    int *beg, *end;
    double *dbeg, *dend;
    int* row_indices = NULL;

    // cout << mtx.shape[0] << ' '<< mtx.shape[1] << ' ' << mtx.shape[2] << endl;

    for(int index = 0; index < mtx.shape[0]; index++ ){

        // cout << "enter for" << endl;
        // if(index % 10000 == 0)
        //     cout << index << endl;
        // if(counter % 10000 == 0)
        //     cout << counter << endl;

        beg = mtx.indices + mtx.indptr[index];
        end = mtx.indices + mtx.indptr[index + 1];
        // for(auto i : mtx.indptr) cout<<i<<" ";
        // cout<<"index"<<index<<endl;
        // cout<<"indptr"<<mtx.indptr[index]<< " "<<mtx.indptr[index+1]<<" "<<endl;
        delete[] row_indices; row_indices = NULL;
        row_indices = new int[end-beg];
        memcpy(row_indices, beg, sizeof(int)*(end-beg));
        for(int i = 0; i < end-beg; i++){
            // cout<<value<<" ";
            auto value = row_indices[i];
            if(element_array[value] == 0)
                counter++;
            element_array[value] = 1;
            // cout<<value<<" ";
        }
        if(counter >= threshold){

            // cout << "if" << endl;

            // element_array.resize(mtx.shape[1],0);
            delete[] element_array; element_array = NULL;
            element_array = new bool[mtx.shape[1]]();

            // cout << "000" << endl;

            counter = 0;
            psize_list[end_psize_list++]=index + 1;
            // p_indptr = mtx.indptr[psize_list[-2]:psize_list[-1]+1]
            beg = mtx.indptr + psize_list[end_psize_list - 2];
            end = mtx.indptr + psize_list[end_psize_list - 1] + 1;
            auto p_indptr_len = end - beg;

            // cout << p_indptr_len << endl;

            int* p_indptr = new int[p_indptr_len];
            memcpy(p_indptr, beg, sizeof(int)*(end-beg));
            // vector<int> p_indptr(beg, end);
            // p_nnz = mtx.data[p_indptr[0]:p_indptr[-1]]

            // cout << "111" << endl;

            dbeg = mtx.data + p_indptr[0];
            dend = mtx.data + p_indptr[p_indptr_len-1];
            // cout << dend - dbeg << endl;
            double* p_nnz = new double[dend - dbeg];
            // cout << "111" << endl;
            memcpy(p_nnz, dbeg, sizeof(double)*(dend - dbeg));
            // vector<double> p_nnz(dbeg, dend);
            // p_indices = mtx.indices[p_indptr[0]:p_indptr[-1]]

            // cout << "222" << endl;

            beg = mtx.indices + p_indptr[0];
            end = mtx.indices + p_indptr[p_indptr_len-1];
            // cout << "length: " << end-beg << ' ' << dend-dbeg << ' ' << p_indptr_len << endl;
            int* p_indices = new int[end-beg];
            memcpy(p_indices, beg, sizeof(int)*(end-beg));
            // vector<int> p_indices(beg, end); 
            int offset = p_indptr[0];
            // p_indptr = p_indptr - offset
            for(int i = 0; i < p_indptr_len; i++)
                p_indptr[i] -= offset;
            // pm = scipy.sparse.csr_matrix((p_nnz, p_indices, p_indptr), shape=(psize_list[-1]-psize_list[-2], mtx.shape[1]))

            // cout << "333" << endl;

            int pm_shape[3];
            pm_shape[0] = psize_list[end_psize_list-1] - psize_list[end_psize_list- 2];
            // cout << "333" << endl;
            pm_shape[1] = mtx.shape[1];
            pm_shape[2] = dend - dbeg;
            // cout << pm_shape[0] << ' ' << pm_shape[1] << ' ' << pm_shape[2] << endl;
            SpM pm(p_nnz, p_indices, p_indptr, pm_shape);
            // pm.check(false);
            delete[] p_nnz, p_indices, p_indptr;
            // cout << "444" << endl;
            // plist.push_back(pm);
            plist[tail++] = pm;
            // cout << "end if" << endl;
        }
    }
    delete[] row_indices, element_array;
    row_indices = NULL; element_array = NULL;

    // cout << "end for" << endl;

    psize_list[end_psize_list++]=mtx.shape[0];
    // p_indptr = mtx.indptr[psize_list[-2]:psize_list[-1]+1]
    // cout << "000" << endl;
    beg = mtx.indptr + psize_list[end_psize_list - 2];
    end = mtx.indptr + psize_list[end_psize_list - 1] + 1;
    int p_indptr_len = end - beg;
    int* p_indptr = new int[p_indptr_len];
    memcpy(p_indptr, beg, sizeof(int)*(end-beg));
    // vector<int> p_indptr(beg, end);
    // p_nnz = mtx.data[p_indptr[0]:p_indptr[-1]]
    dbeg = mtx.data + p_indptr[0];
    dend = mtx.data + p_indptr[p_indptr_len-1];
    double* p_nnz = new double[dend-dbeg];
    memcpy(p_nnz, dbeg, sizeof(double)*(dend-dbeg));
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
    pm_shape[0] = psize_list[end_psize_list- 1] - psize_list[end_psize_list - 2];
    pm_shape[1] = mtx.shape[1];
    pm_shape[2] = dend - dbeg;
    SpM pm(p_nnz, p_indices, p_indptr, pm_shape);
    delete[] p_nnz, p_indices, p_indptr;
    // plist.push_back(pm);
    // cout << "000" << endl;
    plist[tail++] = pm;
    bnum = end_psize_list - 1;
    // cout << "the number of blocks is " << bnum << endl;
    return tail;
}

#endif