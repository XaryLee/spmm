#ifndef _WBSORT_H_
#define _WBSORT_H_

#include <iostream>
#include <vector>
#include <map>
#include "csr.h"

using namespace std;

vector<int> SeqReverse(vector<int> seq){
    auto length = seq.size();
    vector<int> reseq(length);
    for(int i = 0; i < (int)length; i++)
        reseq[seq[i]] = i;
    return reseq;
}

vector<int> gen_rseq(vector<int> seq_bitmap, vector<int> seq_v8){
    vector<int> seq_row;
    for(auto each:seq_v8)
        seq_row.push_back(seq_bitmap[each]);
    return seq_row;
}

vector<int> gen_wseq(vector<int> seq_row, unordered_map<int, int> seq_dict,vector<int> seq_order){
    vector<int> seq_wb;
    for(auto order:seq_order)
        seq_wb.push_back(seq_row[order]);
    return seq_wb;
}

// vector<int> SerialSort(vector<int> seq_bitmap, vector<int> seq_v8, unordered_map<int, int> seq_dict){
//     vector<int> rseq = SeqReverse(gen_rseq(seq_bitmap, seq_v8));
//     vector<int> seq_wb = gen_wseq(rseq, seq_dict);
//     return seq_wb;
// }

vector<vector<int>> SerialSort_block(vector<int> seq_bitmap, vector<int> seq_v8, vector<unordered_map<int, int>> seq_list,vector<vector<int>> seq_order){
    vector<vector<int>> wseq_list;
    vector<int> rseq = SeqReverse(gen_rseq(seq_bitmap, seq_v8));
    int cnt = 0;
    for(auto seq:seq_list){
        vector<int> wseq = gen_wseq(rseq, seq,seq_order[cnt++]);
        wseq_list.push_back(wseq);
    }
    return wseq_list;
}

#endif