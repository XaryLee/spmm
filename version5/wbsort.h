#ifndef _WBSORT_H_
#define _WBSORT_H_

#include <iostream>
#include <vector>
#include <map>
#include <unordered_map>
#include <chrono>
#include <omp.h>
#include "csr.h"

using namespace chrono;

using namespace std;

int* SeqReverse(int* seq, int len){
    // 调用后记得delete[]
    int* reseq = new int[len];
    int MIN_CHUNK_SIZE = 20000;
    if (len < MIN_CHUNK_SIZE || CORENUM == 1) {
        //不并行
        for (int i = 0; i < len; i++) {
            reseq[seq[i]] = i;
        }
    } else {
        //并行
        #pragma omp parallel for num_threads(CORENUM)
           for(int c =0;c<CORENUM;c++)
            {
                for(int i = floor(len * (c) / CORENUM) ; i < floor(len * (c+1) / CORENUM); i++) reseq[seq[i]] = i;
            }
    }
    return reseq;
}

int* gen_rseq(int* seq_bitmap, vector<int> seq_v8){
    // 记得delete[]
    int* seq_row = new int[seq_v8.size()];
    int len = seq_v8.size();
    // #pragma omp parallel for num_threads(CORENUM)
    int MIN_CHUNK_SIZE = 20000;
    if (len < MIN_CHUNK_SIZE || CORENUM == 1) {
        //不并行
        for (int i = 0; i < len; i++) {
            seq_row[i] = seq_bitmap[seq_v8[i]];
        }
    } else {
        //并行
        #pragma omp parallel for num_threads(CORENUM)
            for(int c =0;c<CORENUM;c++)
            {
                for(int i = floor(len * (c) / CORENUM) ; i < floor(len * (c+1) / CORENUM); i++) seq_row[i] = seq_bitmap[seq_v8[i]];
            }
    }
    return seq_row;
}

int* gen_rseq(int* seq_bitmap, int bitmap_len, int* seq_v8, int seq_v8_len){
    // 记得delete[]
    int* seq_row = new int[seq_v8_len];
    for(auto i = 0; i < seq_v8_len; i++)
        seq_row[i] = seq_bitmap[seq_v8[i]];
    return seq_row;
}

void gen_wseq(int* seq_row, vector<int> seq_dict,vector<int> &seq_input){
    for(auto order:seq_dict)
        seq_input.push_back(seq_row[order]);
}

vector<int> gen_wseq(int* seq_row, unordered_map<int, int> seq_dict,vector<int> seq_order){
    vector<int> seq_wb;
    for(auto order:seq_order)
        seq_wb.push_back(seq_row[order]);
    return seq_wb;
}

void SerialSort_block_vec(int* seq_bitmap, int* seq_v8,int end_seq_v8, vector<vector<int>> seq_list,vector<int> &seq_input){
    int* seq_row = gen_rseq(seq_bitmap, 0,seq_v8,end_seq_v8);
    int* rseq = SeqReverse(seq_row, end_seq_v8);
    delete[] seq_row; seq_row = NULL;
    int cnt = 0;
    // gen_wseq可以改成int* 返回wseq_list:vector<int*>

    for(auto seq:seq_list){
        gen_wseq(rseq, seq,seq_input);
    }
    delete[] rseq; rseq = NULL;
}

#endif