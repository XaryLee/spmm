#ifndef _BITMAP_H_
#define _BITMAP_H_

#include <iostream>
#include <math.h>
#include <algorithm>
#include <numeric>
#include "csr.h"

using namespace std;

int** get_scoreboard(SpM &mr, int sectsize){
    // cout << "get_scoreboard" << endl;
    int shape = mr.shape[0];
    int** board = new int*[shape];
    int sectnum = ceil(mr.shape[1]/sectsize);
    int* nn;
    // cout<<mr.shape[0]<<mr.shape[1]<<mr.shape[2];
    cout << "Bitmap Sections = " << sectnum << endl;
    // cout<<"000"<<mr.shape[0]<<endl;
    for(int i = 0; i < mr.shape[0]; i++){
        // cout<<1<<endl;
        // size(score) = sectnum
        int* score = new int[sectnum]();
        if(mr.indptr[i] == mr.indptr[i+1]){
            board[i] = score;
            continue;
        }
        // cout<<2<<endl;
        int beg, end;
        // for n in mr.indices[mr.indptr[i]:mr.indptr[i+1]]
        beg = mr.indptr[i];
        end = mr.indptr[i+1];
        nn = &(mr.indices[beg]);
        // cout<<3<<endl;
        for(int j = 0; j < end-beg; j++){
            int n = nn[j];
            score[(int)(n/sectsize)] = score[(int)(n/sectsize)] + 1;
        }
        // cout<<4<<endl;
        board[i] = score;
        // cout<<5<<endl;
        // if(i%10000==0) cout<<i<<endl;
    }
    // cout << "end for" << endl;
    return board;
}

int blur(int* score, int sectsize, SpM &mr){
    // cout << "blur" << endl;
    int sectnum = ceil(mr.shape[1]/sectsize);
    // size(score) = sectnum
    int index, max = 0;
    for( int i = 0; i < sectnum; i++ )
        if(score[i] > max){
            max = score[i];
            index = i;
        }
    if(max)
        return index + 1;
    else
        return 0;
}

int* argsort(int* v, int len){
    int* idx = new int[len];
    iota(idx, idx+len, 0);
    sort(idx, idx+len, [&v](int i1, int i2){return v[i1] < v[i2];});
    return idx;
}

int* gen_bitorder(int** scores, int sectsize, SpM &mr){
    // cout << "gen_bitorder" << endl;
    int nrows = mr.shape[0]; // 行数
    int* bitseq = new int[nrows];
    // cout<<0<<endl;
    for(int i = 0; i < nrows; i++){
        // cout<<1<<endl;
        auto score = scores[i];
        int bint = blur(score, sectsize, mr);
        // cout<<2<<endl;
        // cout << bint << endl;
        // cout<<3<<endl;
        //cout<< bint<<" ";
        bitseq[i] = bint;
        // cout<<4<<endl;
        }
    // cout << "end for" << endl;
    return argsort(bitseq, nrows);
}

int* bitmap_reorder(SpM &mr, int sectsize){
    // cout << "bitmap_reorder" << endl;
    auto scores = get_scoreboard(mr, sectsize);
    auto row_seq = gen_bitorder(scores, sectsize, mr);
    auto nrows = mr.shape[0];
    for(int i = 0; i < nrows; i++){
        delete[] scores[i];
        scores[i] = NULL;
    }
    delete[] scores;
    scores = NULL;
    // int sectnum = ceil(mr.shape[1]/sectsize);
    // for(auto &score:scores){
    //     delete[] score;
    //     score = NULL;
    // }
    // cout<<"length:"<<row_seq.size()<<endl;
    return row_seq;
}

#endif