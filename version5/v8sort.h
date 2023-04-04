#ifndef _V8SORT_H_
#define _V8SORT_H_

#include <iostream>
#include <vector>
#include "csr.h"
#include "transmat.h"
#include "bitmap.h"
#include<chrono>
using namespace chrono;
using namespace std;

void get_row_lens(SpM &mr,int* &seq){
    //vector<int> rlen;
    for(int r = 0; r < mr.shape[0]; r++)
        seq[r]= (mr.indptr[r+1]-mr.indptr[r]);
}

void get_same_row_lens(int* rlen,int end_rlen, int* &ret){//把长度为33以下的不同长度的行的个数进行了统计
    //vector<int> ret(33,0);
    for(int i=0;i<end_rlen;i++)
        if(rlen[i] < 33)
            ret[rlen[i]]++;
    //return ret;
}

void gen_panel_list(SpM &mr, int* &panel_list ,int &end_panel_list,int panel_num=0){
    // mr.check();
    int nnz = mr.indptr[mr.shape[0]];
    int row = mr.shape[0];
    if(panel_num == 0)
        panel_num = (int)(row / 2048) + 1;
    if(panel_num == 0)
        panel_num = 1;
    int panel_size = (int)(nnz / panel_num) + 1;
    int boundary = panel_size;
    panel_list[end_panel_list++]=0;
    for(int i = 0; i < row; i += 8)
        if(mr.indptr[i] > boundary){
            boundary += panel_size;
            panel_list[end_panel_list++]=i;
        }
    panel_list[end_panel_list++]=row;
    // cout << "the correct number of panels is " << panel_num << endl;
    // cout << "the true number of panels is " << panel_list.size() - 1 << endl;
    //return panel_list;
}

void gen_panel_list(SpM &mr, int row_beg, int row_end, int* &panel_list ,int &end_panel_list,int panel_num=0){
    // mr.check();
    int row = row_end - row_beg;
    // int nnz = mr.indptr[row];
    auto indptr_beg = mr.indptr + row_beg;
    auto indptr_end = mr.indptr + row_end + 1;
    int offset = indptr_beg[0];
    int nnz = indptr_beg[row] - offset;
    if(panel_num == 0)
        panel_num = (int)(row / 2048) + 1;
    if(panel_num == 0)
        panel_num = 1;
    int panel_size = (int)(nnz / panel_num) + 1;
    int boundary = panel_size;
    panel_list[end_panel_list++]=0;
    for(int i = 0; i < row; i += 8)
        if(indptr_beg[i]-offset > boundary){
            boundary += panel_size;
            panel_list[end_panel_list++]=i;
        }
    panel_list[end_panel_list++]=row;
    // cout << "the correct number of panels is " << panel_num << endl;
    // cout << "the true number of panels is " << panel_list.size() - 1 << endl;
    //return panel_list;
}

void panel_sort_nnz(int* &out_seq, int &end_out_seq,vector<int>& spv8_list, SpM &mr, int* panelsize_list, int end_panelsize_list,int panelsize = 2 * 1024){
    int* seq = new int[mr.indptr[mr.shape[0]]];
    //int end_seq=0;
    get_row_lens(mr,seq);
    int* order=new int[mr.shape[0]];
    int end_order = 0;
    int* remain = new int[mr.shape[0]];
    int end_remain=0;
    //vector<int> same_len_rows;
    //vector<int> pseq;
    int count, spv8_len;
    double time_sort=0;//for test
    int* same_len_rows=new int[33];
    int* pseq=new int[mr.shape[0]];
    int end_pseq=0;
    for (int i=0;i<33;i++) same_len_rows[i]=0;
    for(int j = 0; j < end_panelsize_list - 1; j++){
        //vector<int>().swap(remain);
        end_pseq = panelsize_list[j+1]-panelsize_list[j];
        end_order=0;
        end_remain=0;
        //vector<int>().swap(same_len_rows);
        count = 0;
        spv8_len = 0;
        auto begin_sort= high_resolution_clock::now();
        argsort(&seq[0] + panelsize_list[j], end_pseq,pseq);
        //cout<<end_pseq<<endl;
        auto end_sort= high_resolution_clock::now();
        auto duration_sort = duration_cast<microseconds>(end_sort - begin_sort);
        time_sort +=  double(duration_sort.count());
        
        for(int i=0;i<end_pseq ;i++)
            pseq[i] += panelsize_list[j];
            get_same_row_lens(&seq[0] + panelsize_list[j], end_pseq,same_len_rows);
        for(int i = 0; i < 33; i++){
            int begin = count;
            int end = count + same_len_rows[i];
            int boundary = end - same_len_rows[i] % 8;
            spv8_len += boundary - begin;
            
            // cout << "boundary = " << boundary << endl;
            // cout << "begin = " << begin << endl;
            // cout << "spv8_len = " << spv8_len << endl;

            for(int i = begin; i < boundary; i++)
                order[end_order++]=pseq[i];
            for(int i = boundary; i < end; i++)
                remain[end_remain++]=pseq[i];
            count += same_len_rows[i];
            same_len_rows[i]=0;//added
        }
        spv8_list.push_back(spv8_len);
        for(int i = count; i < (int)(end_pseq); i++)
            remain[end_remain++]=pseq[i];

        auto add_seq = order;
        for (int i=0;i<end_order;i++)
        {
            out_seq[end_out_seq++]=order[i];
        }
        for (int i=0;i<end_remain;i++)
        {
            out_seq[end_out_seq++]=remain[i];
        }
        // for(auto each:remain)
        //     add_seq.push_back(each);
        // for(auto each:add_seq)
        //     out_seq[end_out_seq++]=each;
    }
    delete[] pseq; pseq=NULL;
    delete[] seq;seq=NULL;
    delete[] same_len_rows;same_len_rows=NULL;
    delete[] order;order = NULL;
    delete[] remain;remain=NULL;
    // cout<<"----time for sort:"<<time_sort<<"----"<<endl;
}

void panel_sort_nnz(int* &out_seq, int &end_out_seq,vector<int>& spv8_list, SpM &mr, int row_beg, int row_end, int* panelsize_list, int end_panelsize_list,int panelsize = 2 * 1024){
    int nrows = row_end - row_beg;
    auto indptr_beg = mr.indptr + row_beg;
    auto indptr_end = mr.indptr + row_end + 1;
    int offset = indptr_beg[0];
    int* seq = new int[indptr_beg[nrows]-offset];
    //int end_seq=0;
    for(int r = 0; r < nrows; r++)
        seq[r] = indptr_beg[r+1] - indptr_beg[r];
    int* order=new int[nrows];
    int end_order = 0;
    int* remain = new int[nrows];
    int end_remain=0;
    //vector<int> same_len_rows;
    //vector<int> pseq;
    int count, spv8_len;
    double time_sort=0;//for test
    int* same_len_rows=new int[33];
    int* pseq=new int[nrows];
    int end_pseq=0;
    for (int i=0;i<33;i++) same_len_rows[i]=0;
    for(int j = 0; j < end_panelsize_list - 1; j++){
        //vector<int>().swap(remain);
        end_pseq = panelsize_list[j+1]-panelsize_list[j];
        end_order=0;
        end_remain=0;
        //vector<int>().swap(same_len_rows);
        count = 0;
        spv8_len = 0;
        auto begin_sort= high_resolution_clock::now();
        argsort(&seq[0] + panelsize_list[j], end_pseq,pseq);
        //cout<<end_pseq<<endl;
        auto end_sort= high_resolution_clock::now();
        auto duration_sort = duration_cast<microseconds>(end_sort - begin_sort);
        time_sort +=  double(duration_sort.count());
        
        for(int i=0;i<end_pseq ;i++)
            pseq[i] += panelsize_list[j];
            get_same_row_lens(&seq[0] + panelsize_list[j], end_pseq,same_len_rows);
        for(int i = 0; i < 33; i++){
            int begin = count;
            int end = count + same_len_rows[i];
            int boundary = end - same_len_rows[i] % 8;
            spv8_len += boundary - begin;
            
            // cout << "boundary = " << boundary << endl;
            // cout << "begin = " << begin << endl;
            // cout << "spv8_len = " << spv8_len << endl;

            for(int i = begin; i < boundary; i++)
                order[end_order++]=pseq[i];
            for(int i = boundary; i < end; i++)
                remain[end_remain++]=pseq[i];
            count += same_len_rows[i];
            same_len_rows[i]=0;//added
        }
        spv8_list.push_back(spv8_len);
        for(int i = count; i < (int)(end_pseq); i++)
            remain[end_remain++]=pseq[i];

        auto add_seq = order;
        for (int i=0;i<end_order;i++)
        {
            out_seq[end_out_seq++]=order[i];
        }
        for (int i=0;i<end_remain;i++)
        {
            out_seq[end_out_seq++]=remain[i];
        }
        // for(auto each:remain)
        //     add_seq.push_back(each);
        // for(auto each:add_seq)
        //     out_seq[end_out_seq++]=each;
    }
    delete[] pseq; pseq=NULL;
    delete[] seq;seq=NULL;
    delete[] same_len_rows;same_len_rows=NULL;
    delete[] order;order = NULL;
    delete[] remain;remain=NULL;
    // cout<<"----time for sort:"<<time_sort<<"----"<<endl;
}

SpM transpose_v8(SpM &mr, vector<int> spv8_list, int panelsize=2*1024){
    auto rowptr = mr.indptr;
    auto colidx = mr.indices;
    auto data = mr.data;
    auto new_rowptr = mr.indptr;
    int* new_colidx = new int[mr.shape[2]];
    int new_colidx_tail = 0;
    double* new_data = new double[mr.shape[2]];
    int new_data_tail = 0;
    int base = 0;
    for(int i = 0; i < (int)(spv8_list.size()); i++){
        int count = spv8_list[i];
        for(int row_index = base; row_index < base + count; row_index += 8){
            int rowlen = rowptr[row_index + 1] - rowptr[row_index];
            int nnz_index = rowptr[row_index];
            int *add_data = new int[rowlen*8]();
            int *add_colidx = new int[rowlen*8]();
            for(int row = 0; row < rowlen; row++)
                for(int col = 0; col < 8; col++){
                    add_data[8 * row + col] = data[nnz_index + col * rowlen + row];
                    add_colidx[8 * row + col] = colidx[nnz_index + col * rowlen + row];
                }
            for(int j = 0; j < rowlen*8; j++){
                new_data[new_data_tail] = add_data[j];
                new_data_tail++;
            }
            for(int j = 0; j < rowlen*8; j++){
                new_colidx[new_colidx_tail] = add_colidx[j];
                new_colidx_tail++;
            }
            delete[] add_data;
            delete[] add_colidx;
        }
        if(i == (int)(spv8_list.size()-1)){
            int row_begin = base + count;
            int nnz_begin = rowptr[row_begin];
            for(int k = nnz_begin; k < (int)(mr.shape[2]); k++){
                new_data[new_data_tail] = data[k];
                new_data_tail++;
            }
            for(int k = nnz_begin; k < (int)(mr.shape[2]); k++){
                new_colidx[new_colidx_tail] = colidx[k];
                new_colidx_tail++;
            }
        }
        else{
            int row_begin = base + count;
            int row_end = base + panelsize;
            int nnz_begin = rowptr[row_begin];
            int nnz_end = rowptr[row_end];
            for(int k = nnz_begin; k < nnz_end; k++){
                new_data[new_data_tail] = data[k];
                new_data_tail++;
            }
            for(int k = nnz_begin; k < nnz_end; k++){
                new_colidx[new_colidx_tail] = colidx[k];
                new_colidx_tail++;
            }
        }
        base += panelsize;
    }
    SpM new_mr(new_data, new_colidx, new_rowptr, mr.shape);
    delete[] new_data, new_colidx;
    return new_mr;
}

SpM transpose_spv8_nnz(SpM &mr, vector<int> spv8_list, int* panelsize_list){
    // cout << "enter transpose_spv8_nnz" << endl;
    auto rowptr = mr.indptr;

    // cout << "rowptr:\n";
    // for(auto each:rowptr) cout << each << ' '; cout << endl;

    auto colidx = mr.indices;
    auto data = mr.data;
    auto new_rowptr = mr.indptr;
    int* new_colidx = new int[mr.shape[2]];
    int new_colidx_tail = 0;
    double* new_data = new double[mr.shape[2]];
    int new_data_tail = 0;
    int base = 0;
    int times = 0;
    int times1 = 0;
    int times2 = 0;
    // cout << spv8_list.size() << endl;
    for(int i = 0; i < (int)(spv8_list.size()); i++){
        // cout << i << endl;
        int count = spv8_list[i];
        int panelsize = panelsize_list[i+1] - panelsize_list[i];
        for(int row_index = base; row_index < base + count; row_index += 8){
            int rowlen = rowptr[row_index + 1] - rowptr[row_index];

            // cout << "rowlen: " << rowlen << endl;

            int nnz_index = rowptr[row_index];
            double *add_data = new double[rowlen*8]();
            int *add_colidx = new int[rowlen*8]();
            // vector<double> add_data(rowlen*8);
            // vector<int> add_colidx(rowlen*8);

            for(int row = 0; row < rowlen; row++)
                for(int col = 0; col < 8; col++){
                    add_data[8 * row + col] = data[nnz_index + col * rowlen + row];
                    add_colidx[8 * row + col] = colidx[nnz_index + col * rowlen + row];
                }
            for(int j = 0; j < rowlen*8; j++){
                new_data[new_data_tail++] = add_data[j];
                // new_data.push_back(add_data[j]);
            //for(int j = 0; j < rowlen*8; j++){
                new_colidx[new_colidx_tail++] = add_colidx[j];
                // new_colidx.push_back(add_colidx[j]);
                times++;
            }

            // cout << "type 1, " << new_colidx.size() << endl;

            delete[] add_data;
            delete[] add_colidx;
        }
        // cout << "end for" << endl;
        if(i == (int)(spv8_list.size()-1)){
            // cout << "if" << endl;
            int row_begin = base + count;
            int nnz_begin = rowptr[row_begin];

            // cout << "row_begin = " << row_begin << endl;
            // cout << "nnz_begin = " << nnz_begin << endl;
            // cout << "data_size = " << data.size() << endl;
            for(int k = nnz_begin; k < mr.shape[2]; k++){
                new_data[new_data_tail] = data[k];
                new_data_tail++;
                times1++;
                new_colidx[new_colidx_tail++] = colidx[k];
                // new_colidx.push_back(colidx[k]);
            }

            // cout << "type 2, " << new_colidx.size() << endl;
            // cout << "end if" << endl;
        }
        else{
            int row_begin = base + count;
            int row_end = base + panelsize;
            int nnz_begin = rowptr[row_begin];
            int nnz_end = rowptr[row_end];
            for(int k = nnz_begin; k < nnz_end; k++){
                new_data[new_data_tail++] = data[k];
                
                times2++;
                new_colidx[new_colidx_tail++] = colidx[k];
                // new_colidx.push_back(colidx[k]);
            }

            // cout << "type 3, " << new_colidx.size() << endl;

        }
        base += panelsize;
    }
    // cout << "end for" << endl;
    SpM new_mr(new_data, new_colidx, new_rowptr, mr.shape);
    delete[] new_data, new_colidx;
    return new_mr;
}

#endif
