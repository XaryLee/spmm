#ifndef _BITMAP_H_
#define _BITMAP_H_

#include <iostream>
#include <vector>
#include <math.h>
#include <algorithm>
#include <numeric>
#include "csr.h"

using namespace std;

vector<vector<int>> get_scoreboard(SpM mr, int sectsize, int thres=0){
    cout << "get_scoreboard" << endl;
    vector<vector<int>> board;
    int sectnum = ceil(mr.shape[1]/sectsize);
    // cout<<mr.shape[0]<<mr.shape[1]<<mr.shape[2];
    cout << "Bitmap Sections = " << sectnum << endl;
    for( int i = 0; i < mr.shape[0]; i++){
        vector<int> score(sectnum);
        // size(score) = sectnum
        if(mr.indptr[i] == mr.indptr[i+1]){
            board.push_back(score);
            continue;
        }
        vector<int>::const_iterator beg, end;
        // for n in mr.indices[mr.indptr[i]:mr.indptr[i+1]]
        beg = mr.indices.begin() + mr.indptr[i];
        end = mr.indices.begin() + mr.indptr[i+1];
        vector<int> nn(beg, end);
        for(int n:nn)
            score[(int)(n/sectsize)] = score[(int)(n/sectsize)] + 1;
        board.push_back(score);
    }
    cout << "end for" << endl;
    vector<vector<int>>(board).swap(board);
    return board;
}

int blur(vector<int> score, int sectsize, SpM mr){
    int sectnum = ceil(mr.shape[1]/sectsize);
    // size(score) = sectnum
    int index, max = 0;
    for( int i = 0; i < sectnum; i++ )
        if(score[i] > max){
            max = score[i];
            index = i;
        }
    if(max)
        return index;
    else
        return 0;
}

vector<int> argsort(const vector<int> v){
    vector<int> idx(v.size());
    iota(idx.begin(), idx.end(), 0);
    sort(idx.begin(), idx.end(), [&v](int i1, int i2){return v[i1] < v[i2];});
    return idx;
}

vector<int> gen_bitorder(vector<vector<int>> scores, int sectsize, SpM mr){
    cout << "gen_bitorder" << endl;
    vector<int> bitseq;
    for(auto score:scores)
        bitseq.push_back(blur(score, sectsize, mr));
    return argsort(bitseq);
}

vector<int> bitmap_reorder(SpM mr, int sectsize){
    cout << "bitmap_reorder" << endl;
    auto scores = get_scoreboard(mr, sectsize);
    auto row_seq = gen_bitorder(scores, sectsize, mr);
    // cout<<"length:"<<row_seq.size()<<endl;
    return row_seq;
}

// int main(){
//     return 0;
// }

#endif