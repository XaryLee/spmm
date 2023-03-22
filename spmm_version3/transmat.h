#ifndef _TRANSMAT_H_
#define _TRANSMAT_H_
#include <iostream>
#include <vector>
#include "csr.h"

using namespace std;

SpM reorder_row(SpM &mtx, int* seq){
    // len(seq) = mtx.shape[0]
    // mtx.check();

    // cout << "seq: " << endl;
    // for(int i = 0; i < min(100, mtx.shape[0]); i++)
    //     cout << seq[i] << ' ';
    // cout << endl;
    cout << mtx.shape[2] << endl;
    double* nnz = new double[mtx.shape[2]];
    int* colidx = new int[mtx.shape[2]];
    int* rowptr = new int[mtx.shape[0]+1]();
    int tail = 0;
    int rowptr_tail = 0;
    rowptr[rowptr_tail++] = 0;
    cout << "for" << endl;
    for(int i = 0; i < mtx.shape[0]; i++){
        if(i % 10000 == 0)
            cout << i << endl;
        int s = seq[i];
        // cout<<"(i,s)"<<i<<" "<<s<<endl;
        if(mtx.indptr[s] == mtx.indptr[s+1]){
            // cout << "if" << endl;
            rowptr[rowptr_tail] = rowptr[rowptr_tail-1];
            rowptr_tail++;
        }
        else{
            // cout << "else" <<" ";
            rowptr[rowptr_tail] = (rowptr[rowptr_tail-1] + (mtx.indptr[s+1] - mtx.indptr[s]));
            rowptr_tail++;
            // cout<< "c1"<<" "<<mtx.indptr[s]<<" "<<mtx.indptr[s+1]<<" ";
            for( int j = mtx.indptr[s]; j < mtx.indptr[s+1]; j++ ){
                // if(tail<258535 && tail > 200000) cout<<tail<<" ";
                if(tail == 228856){
                    cout<<"enter228856"<<endl;
                    cout<<mtx.data[j]<<" "<<mtx.indices[j]<<endl;
                    cout<<"000"<<endl;
                    colidx[tail] = mtx.indices[j];
                    cout<<"111"<<endl;
                    // for(int k = 0;k<tail;k++) if(nnz[k]!=1) cout<<tail<<endl;
                    nnz[tail] = mtx.data[j];
                    cout<<"222"<<endl;
                    tail++;
                    cout<<"333"<<endl;
                }
                else{
                    nnz[tail] = mtx.data[j];
                    colidx[tail] = mtx.indices[j];
                    tail++;
                }

                // cout<<tail<<" "<<mtx.shape[2]<<endl;
                // 这里可以用memcpy加速
            }
            // cout<<"c2"<<endl;
        }
    }
    cout << "end for" << endl;
    SpM mr(nnz, colidx, rowptr, mtx.shape);
    // cout << mr.shape[0] << ' ' << mr.shape[1] << ' ' << mr.shape[2] << endl;
    delete[] nnz, colidx, rowptr;
    // cout << mr.shape[0] << ' ' << mr.shape[1] << ' ' << mr.shape[2] << endl;
    return mr;
}

SpM reorder_row_test(SpM &mtx, int* seq){
    // len(seq) = mtx.shape[0]
    // mtx.check();

    // cout << "seq: " << endl;
    // for(int i = 0; i < min(100, mtx.shape[0]); i++)
    //     cout << seq[i] << ' ';
    // cout << endl;
    
    double* nnz = new double[mtx.shape[2]];
    int* colidx = new int[mtx.shape[2]];
    int* rowptr = new int[mtx.shape[0]+1]();
    int tail = 0;
    int rowptr_tail = 0;
    rowptr[rowptr_tail++] = 0;
    cout << "for" << endl;
        int i = 14806;
        int s = seq[i];
        cout<<"(i,s)"<<i<<" "<<s<<endl;
        if(mtx.indptr[s] == mtx.indptr[s+1]){
            cout << "if" << endl;
            rowptr[rowptr_tail] = rowptr[rowptr_tail-1];
            rowptr_tail++;
        }
        else{
            cout << "else" << endl;
            rowptr[rowptr_tail] = (rowptr[rowptr_tail-1] + (mtx.indptr[s+1] - mtx.indptr[s]));
            rowptr_tail++;
            for( int j = mtx.indptr[s]; j < mtx.indptr[s+1]; j++ ){
                nnz[tail] = (mtx.data[j]);
                colidx[tail] = (mtx.indices[j]);
                tail++;
                // 这里可以用memcpy加速
            }
        }
    
    cout << "end for" << endl;
    SpM mr(nnz, colidx, rowptr, mtx.shape);
    // cout << mr.shape[0] << ' ' << mr.shape[1] << ' ' << mr.shape[2] << endl;
    delete[] nnz, colidx, rowptr;
    // cout << mr.shape[0] << ' ' << mr.shape[1] << ' ' << mr.shape[2] << endl;
    return mr;
}


int gen_new_panels(SpM &mtx, SpM* &plist, vector<int> &psize_list, int &bnum){

    // 返回值是plist的长度

    // plist = []
    // psize_list = [0]
    //!!!!!
    // for(auto i:mtx.indptr) cout<<i<<" ";
    cout << "enter gen_new_panels" << endl;
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
            psize_list.push_back(index + 1);
            // p_indptr = mtx.indptr[psize_list[-2]:psize_list[-1]+1]
            beg = mtx.indptr + *(psize_list.end() - 2);
            end = mtx.indptr + *(psize_list.end() - 1) + 1;
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
            pm_shape[0] = *(psize_list.end() - 1) - *(psize_list.end() - 2);
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

    psize_list.push_back(mtx.shape[0]);
    // p_indptr = mtx.indptr[psize_list[-2]:psize_list[-1]+1]
    // cout << "000" << endl;
    beg = mtx.indptr + *(psize_list.end() - 2);
    end = mtx.indptr + *(psize_list.end() - 1) + 1;
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
    pm_shape[0] = *(psize_list.end() - 1) - *(psize_list.end() - 2);
    pm_shape[1] = mtx.shape[1];
    pm_shape[2] = dend - dbeg;
    SpM pm(p_nnz, p_indices, p_indptr, pm_shape);
    delete[] p_nnz, p_indices, p_indptr;
    // plist.push_back(pm);
    // cout << "000" << endl;
    plist[tail++] = pm;
    bnum = psize_list.size() - 1;
    cout << "the number of blocks is " << bnum << endl;
    return tail;
}

#endif