#include <iostream>
#include <vector>

#ifndef _CSR_H_
#define _CSR_H_
#include "csr.h"
#endif

#include "transmat.h"
#include "bitmap.h"
using namespace std;

vector<int> get_row_lens(SpM mr){
    vector<int> rlen;
    for(int r = 0; r < mr.shape[0]; r++)
        rlen.push_back(mr.indptr[r+1]-mr.indptr[r]);
    return rlen;
}

vector<int> get_same_row_lens(vector<int> rlen){
    vector<int> ret(0, 33);
    for(auto &x:ret)
        if(x < 33)
            ret[x]++;
    return ret;
}

vector<int> gen_panel_list(SpM mr, int panel_num=0){
    int nnz = *(mr.indptr.end() - 1);
    int row = mr.indptr.size() - 1;
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
    cout << "the correct number of panels is" << panel_num << endl;
    cout << "the true number of panels is" << panel_list.size() - 1 << endl;
    return panel_list;
}

void panel_sort(vector<int>& out_seq, vector<int>& spv8_list, SpM mr, int panelsize = 2 * 1024){
    cout << "Sort Panel Size " << panelsize << endl;
    auto seq = get_row_lens(mr);
    int pnum = (int)(seq.size() / panelsize);
    int iterations = 0;
    if(pnum == (float)(seq.size()) / (float)(panelsize))
        iterations = pnum;
    else
        iterations = pnum + 1;
    for(int p = 0; p < iterations; p++){
        int count = 0;
        vector<int> order;
        vector<int> remain;
        vector<int> same_len_rows;
        vector<int> pseq;
        int spv8_len = 0;
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

void panel_sort_nnz(vector<int>& out_seq, vector<int>& spv8_list, SpM mr, vector<int> panelsize_list, int panelsize = 2 * 1024){
    auto seq = get_row_lens(mr);
    for(int j = 0; j < panelsize_list.size() - 1; j++){
        int count = 0;
        vector<int> order;
        vector<int> remain;
        vector<int> same_len_rows;
        vector<int> pseq;
        int spv8_len = 0;
        pseq = argsort(vector<int>(seq.begin() + panelsize_list[j], seq.begin() + panelsize_list[j+1]));
        for(auto &each:pseq)
            each += panelsize_list[j];
        same_len_rows = get_same_row_lens(vector<int>(seq.begin() + panelsize_list[j], seq.begin() + panelsize_list[j+1]));
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
        for(int i = count; i < (int)(pseq.size()); i++)
            remain.push_back(pseq[i]);
        auto add_seq = order;
        for(auto each:remain)
            add_seq.push_back(each);
        for(auto each:add_seq)
            out_seq.push_back(each);
    }
}

int main(){
    return 0;
}
