#include<iostream>
#include<vector>
#include<string>
#include<fstream>
#include<map>
#include "csr.h"
#include "bitmap.h"
#include "transmat.h"
#include "v8sort.h"
#include "wbsort.h"

bool origin = 1;

using namespace std;
string get_mname(string origin_name){
    string mname;
    for(int i = 0;i<origin_name.length();i++){
        if(origin_name[i]!=46){
            mname.push_back(origin_name[i]);
        }
        else break;
    }
    return mname;
}

SpM csr_matrix(string path){
    ifstream file(path);
    int num_row, num_col, num_lines;
    while (file.peek() == '%') file.ignore(2048, '\n');

    file >> num_row>> num_col >> num_lines;
    vector<int> shape(3);
    shape[0]=num_row;
    shape[1]=num_col;
    shape[2]=num_lines;

    vector<int> indices;
    vector<int> indptr = {0};
    vector<double> data(num_lines,1);
    int cnt = 1;
    int times = 0;
    for (int l = 0; l < num_lines; l++)
    {
        int row, col;
        file >> row >> col;
        indices.push_back(row);
        if(cnt != int(col)||l == num_lines - 1){
            if(cnt != int(col)&&cnt != int(col) - 1){
                int div = col - cnt - 1;
                for(int i = 0;i < div;i++){
                    indptr.push_back(indptr.back());
                    cnt++;
                }
            }
            indptr.push_back(indptr.back() + times);
            times = 0;
            cnt++;
        }
        times++;
    }
    indptr.push_back(shape[2]);
    file.close();
    // cout<<mr.shape[0]<<mr.shape[1];
    return SpM(data,indices,indptr,shape);
}

void gen_serail_origin(SpM mr,SpM &new_mr,map<int,int> &seq_dict){
    int seq_cnt = 0;
    auto ind_iter = mr.indices.begin();
    while(ind_iter != mr.indices.end()){
        auto key_iter = seq_dict.find(*ind_iter);
        if(key_iter == seq_dict.end()){
            seq_dict[*ind_iter] = seq_cnt;
            seq_cnt++;
            if(seq_cnt == mr.shape[1]) break;
        }
        ind_iter++;
    }
    vector<int> new_seq;
    ind_iter = mr.indices.begin();
    while(ind_iter != mr.indices.end()){
        new_seq.push_back(seq_dict[*ind_iter]);
        ind_iter++;
    }
    new_mr = SpM(mr.data,new_seq,mr.indptr,mr.shape);
}
//是不是有点问题
void gen_serail(SpM mr,SpM &new_mr,map<int,int> &seq_dict){
    int seq_cnt = 0;
    vector<int> exist_colidx(mr.shape[1],0);
    for(auto iter = mr.indices.begin();iter != mr.indices.end();iter++) exist_colidx[*iter] = 1;
    for(int i=0;i<mr.shape[1];i++){
        if(exist_colidx[i] == 1){
            seq_dict[i] = seq_cnt;
            seq_cnt = seq_cnt + 1;
        }
    }
    vector<int> new_seq;
    for(auto colidx=mr.indices.begin();colidx != mr.indices.end();colidx++) new_seq.push_back(seq_dict[*colidx]);
    new_mr = SpM(mr.data,new_seq,mr.indptr,mr.shape);
}

void gen_trace_formats(SpM mr,vector<int> &seq_input,vector<int> &rseq,vector<int> &sseq,int &bnum,bool bitmap= 0,bool v8= 0,bool v8_sort= 0,bool serial= 0){
    vector<int> seq_bitmap;
    SpM mr_bitmap;
    if(bitmap){
        int sect = 4;
        seq_bitmap = bitmap_reorder(mr,sect);
        // for(auto i:seq_bitmap) cout<<i<<" ";
        mr_bitmap = reorder_row(mr,seq_bitmap);
        for(auto i:mr_bitmap.indptr) cout<<i<<" ";
    }
    else for(int i = 0;i < seq_bitmap.size();i++) seq_bitmap[i]=i;


    vector<int> seq_v8_block;
    int offset =0;
    vector<map<int,int>> bseq_list;
    vector<int> bserial_colidx;
    vector<double> bserial_data;
    vector<int> bserial_indptr = {0};
    int indptr_offset = 0;

    int panel_size = 2048;

    vector<SpM> regions;
    vector<int> bsize_list;
    gen_new_panels(mr_bitmap,regions,bsize_list,bnum);

    int cnt = 0;
    vector<vector<int>> spv8_lists;
    vector<vector<int>> panelsize_list;
    for(auto r = regions.begin();r!=regions.end();r++){
        vector<int> seq_v8;
        vector<int> spv8_list;
        if(v8){
            vector<int> add_panelsize_list = gen_panel_list(*r);
            panel_sort_nnz(seq_v8,spv8_list,*r,add_panelsize_list);
            SpM tmp_r = reorder_row(*r,seq_v8);
            *r = transpose_spv8_nnz(tmp_r,spv8_list,add_panelsize_list);
            vector<int> ex_in;
            ex_in.assign(add_panelsize_list.begin(),add_panelsize_list.end()-1);
            panelsize_list.push_back(ex_in);
        }
        else{
            if(v8_sort){
                panel_sort(seq_v8,spv8_list,*r,panel_size);
                *r = reorder_row(*r,seq_v8);
            }
            else{
                spv8_list = {};
                seq_v8.resize(r->indptr.size()-1);
                for(int i = 0;i < seq_v8.size();i++) seq_v8[i] = i;
            }
        }
        for(auto iter = seq_v8.begin();iter!=seq_v8.end();iter++) seq_v8_block.push_back(*iter + bsize_list[cnt]);
        for(auto iter = r->data.begin();iter != r->data.end();iter++) bserial_data.push_back(*iter);
        for(int i = 1;i < r->indptr.size();i++) bserial_indptr.push_back(r->indptr[i]+indptr_offset);
        indptr_offset = indptr_offset + r->indptr[r->indptr.size()-1];
        spv8_lists.push_back(spv8_list);
        cnt = cnt + 1;

        SpM smr;
        map<int,int> seq;
        if(serial){
            if(origin) gen_serail_origin(*r,smr,seq);
            else gen_serail(*r,smr,seq);
        }
        else{
            smr = *r;
            seq = {};
            for(int i = 0;i < r->shape[1];i++) seq[i] = i;
        }

        bseq_list.push_back(seq);

        for(auto cc = smr.indices.begin();cc != smr.indices.end();cc++) bserial_colidx.push_back(*cc + offset);
        offset = offset + seq.size();
    }
    //这里是要一维的数组
    vector<vector<int>> seq_temp = SerialSort_block(seq_bitmap,seq_v8_block,bseq_list);
    seq_input.resize(0);
    for(auto i = seq_temp.begin();i != seq_temp.end();i++){
        for(auto j = i->begin();j != i->end();j++){
            seq_input.push_back(*j);
        }
    }
    rseq = SeqReverse(gen_rseq(seq_bitmap,seq_v8_block));
    sseq.resize(0);
    for(auto c = bseq_list.begin();c != bseq_list.end();c++){
        for(auto k = c->begin();k != c->end();k++){
            sseq.push_back(k->first);
        }
    }
}

int main(){
    vector<string> mlist;
    ifstream file("matrix.txt");
    string s;
    while(getline(file,s)){
        mlist.push_back(s);
    }
    file.close();

    // vector<int> bnum_list={1,4,8,16,24,32};
    // auto bnum_iter = bnum_list.begin();
    // auto bnum_end = bnum_list.end();
    // while(bnum_iter != bnum_end){
    //     int bnum = *bnum_iter;
        
    //     ++bnum_iter;
    // }

    for(auto mlist_iter = mlist.begin();mlist_iter != mlist.end();mlist_iter++){
        string mname = *mlist_iter;
        mname = get_mname(mname);
        string path1 = "mat/mtx/";
        string path2 = "/";
        string path3 = ".mtx";
        string mpath = path1+mname+path2+mname+path3;
        SpM mr = csr_matrix(mpath);
        // gen_serail_origin(mr,new_mr,seq_dict);
        vector<int> seq;
        vector<int> sseq;
        vector<int> rseq;
        int bnum;
        gen_trace_formats(mr,seq,rseq,sseq,bnum,1,1,0,1);
        cout<<"done"<<endl;
        if(sseq.empty()) cout<<"empty";
        else{
        for(int i = 0;i<sseq.size();i++){
            cout<<sseq[i]<<" ";
        }}
    }

    return 0;
}