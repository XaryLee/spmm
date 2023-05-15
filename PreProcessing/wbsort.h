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

int* gen_rseq(int* seq_bitmap, int** seq_v8, int *seq_v8_len,int *bsize,int region_size,int seq_len){
    // 记得delete[]
    int* seq_row = new int[seq_len];
    int tail = 0;
    for(int k=0;k<region_size;k++){
        for(auto i = 0; i < seq_v8_len[k]; i++)
            seq_row[tail++] = seq_bitmap[seq_v8[k][i] + bsize[k]];
        }
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

void SerialSort_block_vec(int* seq_bitmap, int** seq_v8,int *end_seq_v8, int* seq_list,int end_seq_list,int* &seq_input,int &end_seq_input,int *bsize,int region_size,int seq_len){
    int* seq_row = gen_rseq(seq_bitmap,seq_v8,end_seq_v8,bsize,region_size,seq_len);
    int* rseq = SeqReverse(seq_row, seq_len);
    delete[] seq_row; seq_row = NULL;
    int cnt = 0;
    // gen_wseq可以改成int* 返回wseq_list:vector<int*>
    for(int i=0;i<end_seq_list;i++)
        seq_input[i]=rseq[seq_list[i]];
    end_seq_input=end_seq_list;
    // for(auto seq:seq_list){
    //     for(auto order:seq)
    //     seq_input.push_back(rseq[order]);
    // }
    delete[] rseq; rseq = NULL;
}

#endif