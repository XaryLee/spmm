#include <iostream>
#include <vector>
#include <math.h>
#include <algorithm>
#include <numeric>

#ifndef _CSR_H_
#define _CSR_H_
#include "csr.cpp"
#endif

using namespace std;

vector<int*> get_scoreboard(SpM mr, int sectsize, int thres=0){
    vector<int*> board;
    int sectnum = ceil(mr.shape[1]/sectsize);
    cout << "Bitmap Sections = " << sectnum << endl;
    for( int i = 0; i < mr.shape[0]; i++){
        int *score = new int[sectnum]();
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
    return board;
}

int blur(int* score, int sectsize, SpM mr){
    int sectnum = ceil(mr.shape[1]/sectsize);
    // size(score) = sectnum
    int index, max = 0;
    bool flag = 0;
    for( int i = 0; i < sectnum; i++ )
        if(*(score + i) > max){
            flag = 1;
            index = i;
        }
    if(flag)
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

vector<int> gen_bitorder(vector<int*> scores, int sectsize, SpM mr){
    vector<int> bitseq;
    for(auto score:scores)
        bitseq.push_back(blur(score, sectsize, mr));
    return argsort(bitseq);
}

vector<int> bitmap_reorder(SpM mr, int sectsize){
    auto scores = get_scoreboard(mr, sectsize);
    auto row_seq = gen_bitorder(scores, sectsize, mr);
    return row_seq;
}

// int main(){
//     return 0;
// }