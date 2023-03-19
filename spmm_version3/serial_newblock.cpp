#include<iostream>
#include<vector>
#include<string>
#include<fstream>
#include<map>
#include<unordered_map>
#include "csr.h"
#include "bitmap.h"
// #include "transmat.h"
// #include "v8sort.h"
// #include "wbsort.h"

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

    int indices[num_lines];
    int indptr[num_row];
    indptr[0] = 0;
    double data[num_lines];
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
    // cout<<mr.shape[0]<<mr.shape[1];
    return SpM(data,indices,indptr,shape);
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
    int* seq_bitmap;
    SpM mr_bitmap;
    if(bitmap){
        int sect = 8;
        seq_bitmap = bitmap_reorder(mr,sect);
        // mr_bitmap = reorder_row(mr,seq_bitmap);
    }
    else for(int i = 0;i < sizeof(seq_bitmap)/sizeof(*seq_bitmap);i++) seq_bitmap[i]=i;
    for(int i = 0;i < sizeof(seq_bitmap)/sizeof(*seq_bitmap);i++) cout << seq_bitmap[i];

    // vector<int> seq_v8_block;
    // int offset =0;
    // vector<unordered_map<int,int>> bseq_list;
    // vector<int> bserial_colidx;
    // vector<double> bserial_data;
    // vector<int> bserial_indptr = {0};
    // int indptr_offset = 0;

    // int panel_size = 2048;

    // vector<SpM> regions;
    // vector<int> bsize_list{0};

    // // cout << "check mr_bitmap" << endl;
    // // mr_bitmap.check();

    // gen_new_panels(mr_bitmap,regions,bsize_list,bnum);
    // // for(auto i : regions[0].indices) cout<<i<<" ";
    // // for(auto i: regions[0].indptr) cout<<i<<" ";

    // int cnt = 0;
    // vector<vector<int>> spv8_lists;
    // vector<vector<int>> panelsize_list;
    // vector<vector<int>> seq_order;
    // for(auto r = regions.begin();r!=regions.end();r++){
    //     // for(auto i : (*r).indices) cout<<i<<" ";
    //     // for(auto i: (*r).indptr) cout<<i<<" ";
    //     vector<int> seq_v8;
    //     vector<int> spv8_list;
    //     if(v8){
    //         vector<int> add_panelsize_list = gen_panel_list(*r);

    //         // cout << "check *r" << endl;
    //         // (*r).check();

    //         panel_sort_nnz(seq_v8,spv8_list,*r,add_panelsize_list);
            
    //         // cout << "check spv8_list\n";
    //         // for(auto each:spv8_list) cout << each << ' ';
    //         // cout << endl;

    //         // cout << "check panelsize_list\n";
    //         // for(auto each:add_panelsize_list) cout << each << ' ';
    //         // cout << endl;
            

    //         SpM tmp_r = reorder_row(*r,seq_v8);

    //         // cout << "check tmp_r\n";
    //         // tmp_r.check();

    //         // cout << "tmp_r_indices_size:\n";
    //         // cout << tmp_r.indices.size() << '\n';
    //         // cout << "tmp_r_rowptr:\n";
    //         // for(auto each:tmp_r.indptr) cout << each << ' '; cout << endl;

    //         //transpose_spv8_nnz:
    //         //bug here?输出的indice_size为75，应该和输入一样是85
    //         //
    //         *r = transpose_spv8_nnz(tmp_r,spv8_list,add_panelsize_list);
    //         // cout<<(*r).indptr.size()<<" "<<(*r).indices.size()<<endl;
            
    //         // cout << "check transpose\n";
    //         // (*r).check();

    //         vector<int> ex_in;            
    //         ex_in.assign(add_panelsize_list.begin(),add_panelsize_list.end()-1);
    //         panelsize_list.push_back(ex_in);
    //     }
    //     else{
    //         if(v8_sort){
    //             panel_sort(seq_v8,spv8_list,*r,panel_size);
    //             *r = reorder_row(*r,seq_v8);
    //         }
    //         else{
    //             spv8_list = {};
    //             seq_v8.resize(r->indptr.size()-1);
    //             for(int i = 0;i < seq_v8.size();i++) seq_v8[i] = i;
    //         }
    //     }
    //     for(auto iter = seq_v8.begin();iter!=seq_v8.end();iter++) seq_v8_block.push_back(*iter + bsize_list[cnt]);
    //     for(auto iter = r->data.begin();iter != r->data.end();iter++) bserial_data.push_back(*iter);
    //     for(int i = 1;i < r->indptr.size();i++) bserial_indptr.push_back(r->indptr[i]+indptr_offset);
    //     indptr_offset = indptr_offset + r->indptr[r->indptr.size()-1];
    //     spv8_lists.push_back(spv8_list);
    //     cnt = cnt + 1;

    //     SpM smr;
    //     unordered_map<int,int> seq;
    //     vector<int> order;
    //     if(serial){
    //         if(origin) gen_serail_origin(*r,smr,seq,order);
    //         else gen_serail(*r,smr,seq);
    //     }
    //     else{
    //         smr = *r;
    //         seq = {};
    //         for(int i = 0;i < r->shape[1];i++) seq[i] = i;
    //     }

    //     seq_order.push_back(order);
    //     bseq_list.push_back(seq);

    //     for(auto cc = smr.indices.begin();cc != smr.indices.end();cc++) bserial_colidx.push_back(*cc + offset);
    //     offset = offset + seq.size();
    // }
    // smr_out = SpM(bserial_data,bserial_colidx,bserial_indptr,mr.shape);
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