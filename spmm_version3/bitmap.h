#ifndef _BITMAP_H_
#define _BITMAP_H_

#include <iostream>
#include <vector>
#include <math.h>
#include <algorithm>
#include <numeric>
#include <ctime>
#include "csr.h"

using namespace std;

// int blur(int* score, int sectnum, SpM &mr){
//     int index, max = 0;
//     for( int i = 0; i < sectnum; i++ )
//         if(score[i] > max){
//             max = score[i];
//             index = i;
//         }
//     if(max)
//         return index + 1;
//     else
//         return 0;
// }

int* argsort(int* v, int len){
    // clock_t beg = clock();
    int* idx = new int[len];
    iota(idx, idx+len, 0);
    sort(idx, idx+len, [&v](int i1, int i2){return v[i1] < v[i2];});
    // clock_t end = clock();
    // cout << "argsort time: " << (double)(end - beg) / CLOCKS_PER_SEC << endl;
    return idx;
}

vector<int> argsort(const vector<int> v){
    vector<int> idx(v.size());
    iota(idx.begin(), idx.end(), 0);
    sort(idx.begin(), idx.end(), [&v](int i1, int i2){return v[i1] < v[i2];});
    return idx;
}

int* bitmap_reorder(SpM &mr, int sectsize){
    int nrows = mr.shape[0];
    int sectnum = ceil(mr.shape[1]/sectsize);
    int* nn;
    int score;
    int beg, end;
    int* bitseq = new int[nrows];
    int index = 0, max_index = -1, max = 0;
    int* idx = new int[nrows];
    int* board = new int[sectnum+1]();
    for(int i = 0; i < nrows; i++){
        max = 0;
        max_index = -1;
        index = 0;
        score = 0;
        beg = mr.indptr[i];
        end = mr.indptr[i+1];
        nn = mr.indices + beg;
        for(int j = 0; j < end-beg; j++){
            int n = nn[j];
            if((int)(n/sectsize) > index){
                if(score > max){
                    max = score;
                    max_index = index;
                    score = 0;
                }
            }
            index = (int)(n/sectsize);
            score++;
        }
        bitseq[i] = max_index + 1;
        board[bitseq[i]]++;
    }
    int* total = new int[sectnum+1]();
    for(int i = 1; i < sectnum + 1; i++){
        total[i] = total[i-1] + board[i-1];
    }
    delete[] board;
    board = new int[sectnum+1]();
    for(int i = 0; i < nrows; i++){
        idx[total[bitseq[i]] + board[bitseq[i]]] = i;
        board[bitseq[i]]++;
    }
    delete[] board, total;
    return idx;
}

// SpM bitmap_reorder(SpM &mr, int sectsize){
//     int nrows = mr.shape[0];
//     int sectnum = ceil(mr.shape[1]/sectsize);
//     int* nn;
//     int score;
//     int beg, end;
//     int* bitseq = new int[nrows];
//     int index = 0, max_index = -1, max = 0;
//     int* idx = new int[nrows];
//     int* board = new int[sectnum+1]();
//     SpM new_mr(mr.shape);
//     int tail = 0;
//     for(int i = 0; i < nrows; i++){
//         max = 0;
//         max_index = -1;
//         index = 0;
//         score = 0;
//         beg = mr.indptr[i];
//         end = mr.indptr[i+1];
//         nn = mr.indices + beg;
//         for(int j = 0; j < end-beg; j++){
//             int n = nn[j];
//             if((int)(n/sectsize) > index){
//                 if(score > max){
//                     max = score;
//                     max_index = index;
//                     score = 0;
//                 }
//             }
//             index = (int)(n/sectsize);
//             score++;
//         }
//         bitseq[i] = max_index + 1;
//         board[bitseq[i]]++;
//     }
//     int* total = new int[sectnum+1]();
//     for(int i = 1; i < sectnum + 1; i++){
//         total[i] = total[i-1] + board[i-1];
//     }
//     delete[] board;
//     board = new int[sectnum+1]();
//     for(int i = 0; i < nrows; i++){
//         // idx[total[bitseq[i]] + board[bitseq[i]]] = i;
//         int index = total[bitseq[i]] + board[bitseq[i]];
//         board[bitseq[i]]++;
//         new_mr.indptr[index+1] = (new_mr.indptr[index] + (mr.indptr[i+1] - mr.indptr[i]));
//         for( int j = mr.indptr[i]; j < mr.indptr[i+1]; j++ ){
//             // nnz[ptr+j-beg] = mtx.data[j];
//             // colidx[ptr+j-beg] = mtx.indices[j];
//             new_mr.data[tail] = mr.data[j];
//             new_mr.indices[tail] = mr.indices[j];
//             tail++;
//         }
//     }
//     delete[] board, total;
//     return idx;
// }

#endif