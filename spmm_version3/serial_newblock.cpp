#include<iostream>
#include<vector>
#include<string>
#include<fstream>
#include<map>
#include<unordered_map>
#include "csr.h"
#include "bitmap.h"
#include "transmat.h"
#include "v8_sort.h"
// #include "wbsort.h"

#define MAXLENGTH 10000000
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
    int shape[3];
    shape[0]=num_row;
    shape[1]=num_col;
    shape[2]=num_lines;

    int* indices = new int[num_lines];
    int* indptr = new int[num_row];
    indptr[0] = 0;
    double* data = new double[num_lines];
    int cnt = 1;
    int times = 0;
    int index = 1;
    for (int l = 0; l < num_lines; l++)
    {
        data[l] = 1;
        int row, col;
        file >> row >> col;
        indices[l] = row - 1;
        if(cnt != int(col)||l == num_lines - 1){
            if(cnt != int(col)&&cnt != int(col) - 1){
                int div = col - cnt - 1;
                for(int i = 0;i < div;i++){
                    indptr[index++] = indptr[index - 1];
                    cnt++;
                }
            }
            indptr[index++] = indptr[index - 1] + times;
            times = 0;
            cnt++;
        }
        times++;
    }
    indptr[index++] = shape[2];
    file.close();
    auto mr = SpM(data,indices,indptr,shape);
    delete[] data; data = NULL;
    delete[] indices; indices = NULL;
    delete[] indptr; indptr = NULL;
    // cout<<mr.shape[0]<<mr.shape[1];
    return mr;
}

void gen_serail_origin(SpM &mr,SpM &new_mr,unordered_map<int,int> &seq_dict,vector<int> &seq_order){
    int seq_cnt = 0;
    for(int i = 0;i<mr.shape[2];i++){
        int ind_iter = mr.indices[i];
        auto key_iter = seq_dict.find(ind_iter);
        if(key_iter == seq_dict.end()){
            // seq_dict.insert(make_pair(ind_iter,seq_cnt));
            seq_order.push_back(ind_iter);
            seq_dict[ind_iter] = seq_cnt;
            seq_cnt++;
            if(seq_cnt == mr.shape[1]) break;
        }
    }
    int new_seq[mr.shape[2]];
    for(int i = 0;i<mr.shape[2];i++){
        new_seq[i] = seq_dict[mr.indices[i]];
    }
    new_mr = SpM(mr.data,new_seq,mr.indptr,mr.shape);
}

void gen_serail(SpM &mr,SpM &new_mr,unordered_map<int,int> &seq_dict){
    int seq_cnt = 0;
    int exist_colidx[mr.shape[1]];
    for(int i = 0;i<mr.shape[2];i++) exist_colidx[mr.indices[i]] = 1;
    for(int i=0;i<mr.shape[1];i++){
        if(exist_colidx[i] == 1){
            seq_dict[i] = seq_cnt;
            seq_cnt = seq_cnt + 1;
        }
    }
    int new_seq[mr.shape[2]];
    for(int i = 0;i<mr.shape[2];i++) new_seq[i] = seq_dict[mr.indices[i]];
    new_mr = SpM(mr.data,new_seq,mr.indptr,mr.shape);
}

void gen_trace_formats(SpM &mr,vector<int> &seq_input,vector<int> &rseq,vector<int> &sseq,SpM &smr_out,int &bnum,bool bitmap= 0,bool v8= 0,bool v8_sort= 0,bool serial= 0){
    int* seq_bitmap = new int;
    SpM mr_bitmap;
    if(bitmap){
        int sect = 8;
        seq_bitmap = bitmap_reorder(mr,sect);
        mr_bitmap = reorder_row(mr,seq_bitmap);
    }
    else for(int i = 0;i < sizeof(seq_bitmap)/sizeof(*seq_bitmap);i++) seq_bitmap[i]=i;
    // mr_bitmap.check();
    // for(int i = 0;i < mr.shape[0];i++) cout << seq_bitmap[i]<<" ";

    vector<int> seq_v8_block;
    int offset =0;
    vector<unordered_map<int,int>> bseq_list;
    vector<int> bserial_colidx;
    vector<double> bserial_data;
    vector<int> bserial_indptr = {0};
    int indptr_offset = 0;

    int panel_size = 2048;

    vector<SpM> regions;
    vector<int> bsize_list{0};

    cout<<"enter gen_new_panels"<<endl;
    gen_new_panels(mr_bitmap,regions,bsize_list,bnum);
    cout<<"out gen_new_panels"<<endl;
    regions[0].check();

    int cnt = 0;
    vector<vector<int>> spv8_lists;
    vector<vector<int>> panelsize_list;
    vector<vector<int>> seq_order;

    for(int index = 0;index < regions.size();index++){
        vector<int> seq_v8;
        vector<int> spv8_list;
        if(v8){
            vector<int> add_panelsize_list;
            //返回一个add_panelsize_list的长度
            add_panelsize_list = gen_panel_list(regions[index]);

            panel_sort_nnz(seq_v8,spv8_list,regions[index],add_panelsize_list);

            int *seq_v8_arr = new int[seq_v8.size()];
            for(int i = 0;i < seq_v8.size();i++) seq_v8_arr[i] = seq_v8[i];

            SpM tmp_r = reorder_row(regions[index],seq_v8_arr);

            delete[] seq_v8_arr;seq_v8_arr = NULL;

            regions[index] = transpose_spv8_nnz(tmp_r,spv8_list,add_panelsize_list);
            tmp_r.check();

            vector<int> ex_in;            
            ex_in.assign(add_panelsize_list.begin(),add_panelsize_list.end()-1);
            panelsize_list.push_back(ex_in);
        }
        else{
            if(v8_sort){
                panel_sort(seq_v8,spv8_list,regions[index],panel_size);
                int *seq_v8_arr = new int[seq_v8.size()];
                for(int i = 0;i < seq_v8.size();i++) seq_v8_arr[i] = seq_v8[i];
                regions[index] = reorder_row(regions[index],seq_v8_arr);
                delete[] seq_v8_arr;seq_v8_arr = NULL;
            }
            else{
                spv8_list = {};
                seq_v8.resize(regions[index].shape[0]);
                for(int i = 0;i < seq_v8.size();i++) seq_v8[i] = i;
            }
        }
        for(auto iter = seq_v8.begin();iter!=seq_v8.end();iter++) seq_v8_block.push_back(*iter + bsize_list[cnt]);
        for(int i = 0;i < regions[index].shape[2];i++) bserial_data.push_back(regions[index].data[i]);
        for(int i = 1;i < regions[index].shape[0] - 1;i++) bserial_indptr.push_back(regions[index].indptr[i]+indptr_offset);
        indptr_offset = indptr_offset + regions[index].indptr[regions[index].shape[0]];
        spv8_lists.push_back(spv8_list);
        cnt = cnt + 1;

        SpM smr;
        unordered_map<int,int> seq;
        vector<int> order;
        if(serial){
            if(origin) gen_serail_origin(regions[index],smr,seq,order);
            else gen_serail(regions[index],smr,seq);
        }
        else{
            smr = regions[index];
            seq = {};
            for(int i = 0;i < regions[index].shape[1];i++) seq[i] = i;
        }

        seq_order.push_back(order);
        bseq_list.push_back(seq);

        for(int cc = 0;cc < smr.shape[2]; cc++) bserial_colidx.push_back(smr.indices[cc] + offset);
        offset = offset + seq.size();
    }
    double *bserial_data_arr = new double[bserial_data.size()];
    for(int i = 0;i < bserial_data.size();i++) bserial_data_arr[i] = bserial_data[i];
    int *bserial_colidx_arr = new int[bserial_colidx.size()];
    for(int i = 0;i < bserial_colidx.size();i++) bserial_colidx_arr[i] = bserial_colidx[i];
    int *bserial_indptr_arr = new int[bserial_indptr.size()];
    for(int i = 0;i < bserial_indptr.size();i++) bserial_indptr_arr[i] = bserial_indptr[i];

    smr_out = SpM(bserial_data_arr,bserial_colidx_arr,bserial_indptr_arr,mr.shape);

    delete[] bserial_data_arr; bserial_data_arr = NULL;
    delete[] bserial_colidx_arr; bserial_colidx_arr = NULL;
    delete[] bserial_indptr_arr; bserial_indptr_arr = NULL;
    // //这里是要一维的数组
    // vector<vector<int>> seq_temp = SerialSort_block(seq_bitmap,seq_v8_block,bseq_list,seq_order);
    // seq_input.resize(0);
    // for(auto i = seq_temp.begin();i != seq_temp.end();i++){
    //     for(auto j = i->begin();j != i->end();j++){
    //         seq_input.push_back(*j);
    //     }
    // }
    // rseq = SeqReverse(gen_rseq(seq_bitmap,seq_v8_block));
    // sseq.resize(0);
    // cnt = 0;
    // for(auto k : seq_order[cnt]){
    //     sseq.push_back(k);
    // }

    delete[] seq_bitmap; seq_bitmap = NULL;

}

int main(){
    vector<string> mlist;
    ifstream file("D:\\Users\\Desktop\\learning in XJTU\\VSCodeC++\\Spmm\\matrix.txt");
    string s;
    while(getline(file,s)){
        mlist.push_back(s);
    }
    file.close();
    for(auto mlist_iter = mlist.begin();mlist_iter != mlist.end();mlist_iter++){
        string mname = *mlist_iter;
        mname = get_mname(mname);
        string path1 = "D:\\Users\\Desktop\\learning in XJTU\\VSCodeC++\\Spmm\\mat\\mtx\\";
        string path2 = "/";
        string path3 = ".mtx";
        string mpath = path1+mname+path2+mname+path3;
        SpM mr = csr_matrix(mpath);

        // cout << "check csr_matrix" << endl;
        // mr.check();

        // gen_serail_origin(mr,new_mr,seq_dict);
        vector<int> seq;
        vector<int> sseq;
        vector<int> rseq;
        int bnum;
        SpM smr;
        gen_trace_formats(mr,seq,rseq,sseq,smr,bnum,1,1,0,1);
        cout<<"done"<<endl;
        if(sseq.empty()) cout<<"empty";
        else{
        for(int i = 0;i<sseq.size();i++){
            cout<<sseq[i]<<" ";
        }}
    }
    system("pause");

    return 0;
    }