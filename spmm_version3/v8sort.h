#ifndef _V8SORT_H_
#define _V8SORT_H_

#include <iostream>
#include <vector>
#include "csr.h"
#include "transmat.h"
#include "bitmap.h"
using namespace std;

vector<int> get_row_lens(SpM &mr){
    vector<int> rlen;
    for(int r = 0; r < mr.shape[0]; r++)
        rlen.push_back(mr.indptr[r+1]-mr.indptr[r]);
    return rlen;
}

vector<int> get_same_row_lens(vector<int> rlen){
    vector<int> ret(33,0);
    for(auto &x:rlen)
        if(x < 33)
            ret[x]++;
    return ret;
}

vector<int> gen_panel_list(SpM &mr, int panel_num=0){
    // mr.check();
    int nnz = mr.indptr[mr.shape[0]];
    int row = mr.shape[0];
    if(panel_num == 0)
        panel_num = (int)(row / 2048) + 1;
    if(panel_num == 0)
        panel_num = 1;
    int panel_size = (int)(nnz / panel_num) + 1;
    int boundary = panel_size;
    vector<int> panel_list{0};
    for(int i = 0; i < row; i += 8)
        if(mr.indptr[i] > boundary){
            boundary += panel_size;
            panel_list.push_back(i);
        }
    panel_list.push_back(row);
    cout << "the correct number of panels is " << panel_num << endl;
    cout << "the true number of panels is " << panel_list.size() - 1 << endl;
    return panel_list;
}

void panel_sort(vector<int>& out_seq, vector<int>& spv8_list, SpM &mr, int panelsize = 2 * 1024){
    cout << "Sort Panel Size " << panelsize << endl;
    auto seq = get_row_lens(mr);
    int pnum = (int)(seq.size() / panelsize);
    int iterations = 0;
    int count, spv8_len;
    vector<int> order;
    vector<int> remain;
    vector<int> same_len_rows;
    vector<int> pseq;
    if(pnum == (float)(seq.size()) / (float)(panelsize))
        iterations = pnum;
    else
        iterations = pnum + 1;
    for(int p = 0; p < iterations; p++){
        count = 0;   
        spv8_len = 0;
        if(p == pnum){
            pseq = argsort(vector<int>(seq.begin() + p * panelsize, seq.end()));
            for(auto& each:pseq)
                each += p * panelsize;
            same_len_rows = get_same_row_lens(vector<int>(seq.begin() + p * panelsize, seq.end()));
        }
        else{
            auto pseq = argsort(vector<int>(seq.begin() + p * panelsize, seq.begin() + (p+1) * panelsize));
            for(auto& each:pseq)
                each += p * panelsize;
            same_len_rows = get_same_row_lens(vector<int>(seq.begin() + p * panelsize, seq.begin() + (p+1) * panelsize));
        }
        for(int i = 0; i < 33; i++){
            int begin = count;
            int end = count + same_len_rows[i];
            int boundary = end - same_len_rows[i] % 8;
            spv8_len += boundary - begin;
            for(int i = begin; i < boundary; i++)
                order.push_back(pseq[i]);
            for(int i = boundary; i < end; i++)
                remain.push_back(pseq[i]);
            count += same_len_rows[i];
        }
        spv8_list.push_back(spv8_len);
        for(int i = count; i < pseq.size(); i++)
            remain.push_back(pseq[i]);
        auto add_seq = order;
        for(auto each:remain)
            add_seq.push_back(each);
        for(auto each:add_seq)
            out_seq.push_back(each);
    }
}

void panel_sort_nnz(vector<int>& out_seq, vector<int>& spv8_list, SpM &mr, vector<int> panelsize_list, int panelsize = 2 * 1024){
    auto seq = get_row_lens(mr);
    vector<int> order;
    vector<int> remain;
    vector<int> same_len_rows;
    vector<int> pseq;
    int count, spv8_len;
    for(int j = 0; j < panelsize_list.size() - 1; j++){
        count = 0;
        spv8_len = 0;
        pseq = argsort(vector<int>(seq.begin() + panelsize_list[j], seq.begin() + panelsize_list[j+1]));
        for(auto &each:pseq)
            each += panelsize_list[j];
        same_len_rows = get_same_row_lens(vector<int>(seq.begin() + panelsize_list[j], seq.begin() + panelsize_list[j+1]));
        for(int i = 0; i < 33; i++){
            int begin = count;
            int end = count + same_len_rows[i];
            int boundary = end - same_len_rows[i] % 8;
            spv8_len += boundary - begin;
            
            // cout << "boundary = " << boundary << endl;
            // cout << "begin = " << begin << endl;
            // cout << "spv8_len = " << spv8_len << endl;

            for(int i = begin; i < boundary; i++)
                order.push_back(pseq[i]);
            for(int i = boundary; i < end; i++)
                remain.push_back(pseq[i]);
            count += same_len_rows[i];
        }
        spv8_list.push_back(spv8_len);
        for(int i = count; i < (int)(pseq.size()); i++)
            remain.push_back(pseq[i]);
        auto add_seq = order;
        for(auto each:remain)
            add_seq.push_back(each);
        for(auto each:add_seq)
            out_seq.push_back(each);
    }
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

SpM transpose_spv8_nnz(SpM &mr, vector<int> spv8_list, vector<int> panelsize_list){
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
            for(int j = 0; j < rowlen*8; j++)
                new_data[new_data_tail++] = add_data[j];
                // new_data.push_back(add_data[j]);
            for(int j = 0; j < rowlen*8; j++){
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
                // new_data.push_back(data[k]);
            }
            // cout << "000" << endl;
            for(int k = nnz_begin; k < (mr.shape[2]); k++){
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
            for(int k = nnz_begin; k < nnz_end; k++)
                new_data[new_data_tail++] = data[k];
                // new_data.push_back(data[k]);
            for(int k = nnz_begin; k < nnz_end; k++){
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