#include<iostream>
#include<vector>
#include<string>
#include<fstream>
#include<map>
#include <chrono>
#include<unordered_map>
#include "csr.h"
#include "bitmap.h"
#include "transmat.h"
#include "v8sort.h"
#include "wbsort.h"

#define MAXLENGTH 10000000
#define SECT 2048
#define MAXREGION 20

using namespace std;
using namespace chrono;
double time_bitmap=0;
double time_v8=0;
double time_serial = 0;
double time_wbsort = 0;

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
    ifstream fin;
    fin.open(path);
    while (fin.peek() == '%') fin.ignore(2048, '\n');
    fin.ignore(2048,'\n');
    char ch = 0;
    int count=1;
    while(ch!='\n')
    {
      ch=fin.get();
      if(ch==' ') {count++;}
    }
    //cout<<"count="<<count<<endl;
    fin.close();

    fin.open(path);
    int num_row, num_col, num_lines;
    while (fin.peek() == '%') fin.ignore(2048, '\n');
    fin >> num_row>> num_col >> num_lines;
    int shape[3];
    shape[0]=num_row;
    shape[1]=num_col;
    shape[2]=num_lines;
    int* indices = new int[num_lines];
    int* row_ = new int [num_lines];
    int* col_ = new int [num_lines];
    int* indptr = new int[num_row+1];
    int* rowlen = new int [num_row];
    double* data = new double[num_lines];
    
    for (int i = 0; i < num_row; i++)rowlen[i] = 0;
    indptr[0] = 0;
    if(count==2)
    {
        for (int l = 0; l < num_lines; l++)
        {
            data[l] = 1;
            int row, col;
            fin >> row >> col;
            rowlen[row-1] +=1;
            row_[l] = row-1;
            col_[l] = col-1;
        }
    }
    else
    {
        for (int l = 0; l < num_lines; l++)
        {
            data[l] = 1;
            int row, col;
            fin >> row >> col>>count;
            rowlen[row-1] +=1;
            row_[l] = row-1;
            col_[l] = col-1;
        }
    }
    fin.close();
    for(int i = 1; i < num_row; i++){indptr[i] = rowlen[i-1]+indptr[i-1];rowlen[i-1] = 0;}
    rowlen[num_row-1]=0;
    indptr[num_row] = shape[2];
    for (int i =0; i<num_lines ; i++)
    {
        indices[indptr[row_[i]]+rowlen[row_[i]]]=col_[i];
        rowlen[row_[i]] +=1;
    }
    //for(int i=0;i<(num_lines);i++) cout<<indices[i]<<" ";//this line is checked
    auto mr = SpM(data,indices,indptr,shape);
    delete[] data; data = NULL;
    delete[] indices; indices = NULL;
    delete[] rowlen; rowlen = NULL;
    delete[] col_; col_=NULL;
    delete[] indptr; indptr = NULL;
    delete[] row_; row_=NULL;
    // cout <<endl<< row_[40] << ' ' << shape[2] << endl;
    // 
    return mr;
}

// void gen_serial_origin(SpM &mr,SpM &new_mr,unordered_map<int,int> &seq_dict,vector<int> &seq_order){
//     int seq_cnt = 0;
//     vector<bool> bitmap(mr.shape[1],0);
    
//     int new_seq[mr.shape[2]];
//     seq_dict.reserve(mr.shape[1]);
//     seq_dict.max_load_factor(0.5);
//     for(int i = 0;i<mr.shape[2];i++){
//         int ind_iter = mr.indices[i];
//         if(bitmap[ind_iter] == 0){
//             bitmap[ind_iter] = 1;
//             new_seq[i] = seq_cnt;
//             seq_dict.emplace(ind_iter,seq_cnt);
//             seq_cnt += 1;
//         }
//         else{
//             new_seq[i] = seq_dict[ind_iter];
//         }
//     }
//     new_mr = SpM(mr.data,new_seq,mr.indptr,mr.shape);
// }

void gen_serial_origin(SpM &mr, SpM &new_mr, unordered_map<int,int> &seq_dict, vector<int> &seq_order) {
    int seq_cnt = 0;
    int n_cols = mr.shape[1];
    int n_vals = mr.shape[2];

    bool bitmap[n_cols] = {0};
    int* new_seq = new int[n_vals];

    seq_dict.reserve(n_cols);
    seq_dict.max_load_factor(0.5);

    for (int i = 0; i < n_vals; i++) {
        int ind_iter = mr.indices[i];

        if (bitmap[ind_iter] == 0) {
            // auto start1 = chrono::high_resolution_clock::now();
            bitmap[ind_iter] = 1;
            new_seq[i] = seq_cnt;
            seq_dict.emplace(ind_iter, seq_cnt);
            seq_cnt++;
            // auto end1 = chrono::high_resolution_clock::now();
            // auto duration = chrono::duration_cast<chrono::microseconds>(end1 - start1);
            // time1 += duration.count();
        } else {
            // auto start2 = chrono::high_resolution_clock::now();
            new_seq[i] = seq_dict[ind_iter];
            // auto end2 = chrono::high_resolution_clock::now();
            // auto duration = chrono::duration_cast<chrono::microseconds>(end2 - start2);
            // time2 += duration.count();
        }
        // if(bitmap[ind_iter] == 0) seq_dict.emplace(ind_iter, seq_cnt++);
        // bitmap[ind_iter] = 1;
        // new_seq[i] =  seq_dict[ind_iter];
    }

    new_mr = SpM(mr.data, new_seq, mr.indptr, mr.shape);
    delete[] new_seq;
}

void gen_serial_origin_vec(SpM &mr,SpM &new_mr,vector<int> &seq_key){
    int seq_cnt = 0;
    vector<int> bitmap(mr.shape[1],-1);
    int* new_seq = new int[mr.shape[2]];
    for(int i = 0;i<mr.shape[2];i++){
        int ind_iter = mr.indices[i];
        if(bitmap[ind_iter] == -1){
            bitmap[ind_iter] = seq_cnt;
            seq_key.push_back(ind_iter);
            seq_cnt += 1;
        }
        new_seq[i] = bitmap[ind_iter];
    }
    new_mr = SpM(mr.data,new_seq,mr.indptr,mr.shape);
    delete[] new_seq;new_seq=nullptr;
}

// void gen_serial(SpM &mr,SpM &new_mr,unordered_map<int,int> &seq_dict){
//     int seq_cnt = 0;
//     int exist_colidx[mr.shape[1]];
//     for(int i = 0;i<mr.shape[2];i++) exist_colidx[mr.indices[i]] = 1;
//     for(int i=0;i<mr.shape[1];i++){
//         if(exist_colidx[i] == 1){
//             seq_dict[i] = seq_cnt;
//             seq_cnt = seq_cnt + 1;
//         }
//     }
//     int new_seq[mr.shape[2]];
//     for(int i = 0;i<mr.shape[2];i++) new_seq[i] = seq_dict[mr.indices[i]];
//     new_mr = SpM(mr.data,new_seq,mr.indptr,mr.shape);
// }

void gen_trace_formats(SpM &mr,vector<int> &seq_input,vector<int> &rseq,vector<int> &sseq,SpM &smr_out,int &bnum){

    int* seq_bitmap = new int[mr.shape[0]];

    
    
    // SpM mr_bitmap;

    auto begin_bitmap= high_resolution_clock::now();

    int sect = SECT;
    //cout << "mr_indptr: " << mr.indptr[mr.shape[0]] << endl;//test
    //mr.check(false);//test

    // seq_bitmap = bitmap_reorder(mr,sect);

    // SpM mr_bitmap = reorder_row(mr,seq_bitmap);
    // mr_bitmap.check();
    // cout << mr.shape[0] << ' ' << mr.shape[1] << ' ' << mr.shape[2] << endl;
    // cout << reorder_row(mr,seq_bitmap).shape[0] << ' ' << reorder_row(mr,seq_bitmap).shape[1] << ' ' << reorder_row(mr,seq_bitmap).shape[2] << endl;
    // cout << mr_bitmap.shape[0] << ' ' << mr_bitmap.shape[1] << ' ' << mr_bitmap.shape[2] << endl;
    // auto end_bitmap= high_resolution_clock::now();
    // auto duration_bitmap = duration_cast<microseconds>(end_bitmap - begin_bitmap);
    // time_bitmap +=  double(duration_bitmap.count());

    SpM mr_bitmap = bitmap_reorder(mr, sect, seq_bitmap);
    // cout << "done" << endl;
    
    auto end_bitmap= high_resolution_clock::now();
    auto duration_bitmap = duration_cast<microseconds>(end_bitmap - begin_bitmap);
    time_bitmap +=  double(duration_bitmap.count());
    
    auto begin_v8= high_resolution_clock::now();
    // mr_bitmap.check();
    // for(int i = 0;i < mr.shape[0];i++) cout << seq_bitmap[i]<<" ";

    vector<int> seq_v8_block;
    int offset =0;
    // vector<unordered_map<int,int>> bseq_list;
    vector<vector<int>> bseq_list_key;
    vector<vector<int>> bseq_list_value;
    vector<int> bserial_colidx;
    vector<double> bserial_data;
    vector<int> bserial_indptr = {0};
    int indptr_offset = 0;
    int panel_size = 2048;

    SpM* regions = new SpM[MAXREGION];
    int regions_length = 0;
    vector<int> bsize_list{0};

    // cout<<"enter gen_new_panels"<<endl;
    // cout << mr_bitmap.shape[0] << ' ' << mr_bitmap.shape[1] << ' ' << mr_bitmap.shape[2] << endl;
    regions_length = gen_new_panels(mr_bitmap,regions,bsize_list,bnum);
    // cout << "done" << endl;
    // cout<<"out gen_new_panels"<<endl;
    // regions[0].check();

    int cnt = 0;
    vector<vector<int>> spv8_lists;
    vector<vector<int>> panelsize_list;
    vector<int> seq_v8;
    vector<int> spv8_list;
    vector<int> add_panelsize_list;
    vector<vector<int>> seq_key_list;
    // cout << regions_length << endl;
    for(int index = 0;index < regions_length;index++){
        vector<int>().swap(seq_v8);
        vector<int>().swap(spv8_list);
        add_panelsize_list = gen_panel_list(regions[index]);

        panel_sort_nnz(seq_v8,spv8_list,regions[index],add_panelsize_list);

        int *seq_v8_arr = new int[seq_v8.size()];
        for(int i = 0;i < seq_v8.size();i++) seq_v8_arr[i] = seq_v8[i];

        SpM tmp_r = reorder_row(regions[index],seq_v8_arr);

        regions[index] = transpose_spv8_nnz(tmp_r,spv8_list,add_panelsize_list);

        vector<int> ex_in;            
        ex_in.assign(add_panelsize_list.begin(),add_panelsize_list.end()-1);

        panelsize_list.push_back(ex_in);

        for(auto iter = seq_v8.begin();iter!=seq_v8.end();iter++) seq_v8_block.push_back(*iter + bsize_list[cnt]);
        for(int i = 0;i < regions[index].shape[2];i++) bserial_data.push_back(regions[index].data[i]);
        for(int i = 1;i < regions[index].shape[0] - 1;i++) bserial_indptr.push_back(regions[index].indptr[i]+indptr_offset);
        indptr_offset = indptr_offset + regions[index].indptr[regions[index].shape[0]];
        spv8_lists.push_back(spv8_list);
        cnt = cnt + 1;
        
        auto begin_serial= high_resolution_clock::now();

        SpM smr;
        vector<int> seq_key;
        // unordered_map<int,int> seq;

        // gen_serial_origin(regions[index],smr,seq,order);
        gen_serial_origin_vec(regions[index],smr,seq_key);

        auto end_serial= high_resolution_clock::now();
        auto duration_serial = duration_cast<microseconds>(end_serial - begin_serial);
        time_serial += double(duration_serial.count());

        // bseq_list_key.push_back(seq_key);

        for(int cc = 0;cc < smr.shape[2]; cc++) bserial_colidx.push_back(smr.indices[cc] + offset);
        seq_key_list.push_back(seq_key);
        offset = offset + seq_key.size();
        
    }
    auto end_v8= high_resolution_clock::now();
    auto duration_v8 = duration_cast<microseconds>(end_v8 - begin_v8);
    time_v8 = double(duration_v8.count());
    // cout << "end for" << endl;
    auto begin_wbsort = high_resolution_clock::now();
    double *bserial_data_arr = new double[bserial_data.size()];
    for(int i = 0;i < bserial_data.size();i++) bserial_data_arr[i] = bserial_data[i];
    int *bserial_colidx_arr = new int[bserial_colidx.size()];
    for(int i = 0;i < bserial_colidx.size();i++) bserial_colidx_arr[i] = bserial_colidx[i];
    int *bserial_indptr_arr = new int[bserial_indptr.size()];
    for(int i = 0;i < bserial_indptr.size();i++) bserial_indptr_arr[i] = bserial_indptr[i];

    smr_out = SpM(bserial_data_arr,bserial_colidx_arr,bserial_indptr_arr,mr.shape);
    // smr_out.check();
    delete[] bserial_data_arr; bserial_data_arr = NULL;
    delete[] bserial_colidx_arr; bserial_colidx_arr = NULL;
    delete[] bserial_indptr_arr; bserial_indptr_arr = NULL;
    //这里是要一维的数组
    vector<vector<int>> seq_temp = SerialSort_block_vec(seq_bitmap,seq_v8_block,bseq_list_key);
    seq_input.resize(0);
    for(auto i = seq_temp.begin();i != seq_temp.end();i++){
        for(auto j = i->begin();j != i->end();j++){
            seq_input.push_back(*j);
        }
    }
    int* seq_row = gen_rseq(seq_bitmap,seq_v8_block);
    int* rseq_arr;
    rseq_arr = SeqReverse(seq_row, seq_v8_block.size());
    delete[] seq_row; seq_row = NULL;
    rseq = vector<int>(rseq_arr, rseq_arr + seq_v8_block.size());
    delete[] rseq_arr; rseq_arr = NULL;
    sseq.resize(0);
    for(int i=0;i<seq_key_list.size();i++){
        for(auto k : seq_key_list[i]){
            sseq.push_back(k);
        }
    }

    delete[] seq_bitmap; seq_bitmap = NULL;
    auto end_wbsort = high_resolution_clock::now();
    auto duration_wbsort = duration_cast<microseconds>(end_wbsort - begin_wbsort);
    time_wbsort += double(duration_wbsort.count());

}

int main(){
    vector<string> mlist;
    ifstream file("version4/matrix.txt");
    string s;
    while(getline(file,s)){
        mlist.push_back(s);
    }
    file.close();
    
    for(auto mlist_iter = mlist.begin();mlist_iter != mlist.end();mlist_iter++){
        auto begin = high_resolution_clock::now();
        auto begin_read = high_resolution_clock::now();
        string mname = *mlist_iter;
        mname = get_mname(mname);
        string path1 = "version4/mat/mtx/";
        string path2 = "/";
        string path3 = ".mtx";
        string mpath = path1+mname+path2+mname+path3;
        SpM mr = csr_matrix(mpath);
        auto end_read = high_resolution_clock::now();
        auto duration_read = duration_cast<microseconds>(end_read - begin_read);
        
        //mr.check(false);//test
        
        // cout << "check csr_matrix" << endl;
        // mr.check();

        // gen_serial_origin(mr,new_mr,seq_dict);
        vector<int> seq;
        vector<int> sseq;
        vector<int> rseq;
        int bnum;
        SpM smr;
        
        auto begin_transform = high_resolution_clock::now();
        gen_trace_formats(mr,seq,rseq,sseq,smr,bnum);
        auto end_transform = high_resolution_clock::now();
        auto end = high_resolution_clock::now();
        auto duration_transform = duration_cast<microseconds>(end_transform - begin_transform);
        auto duration = duration_cast<microseconds>(end- begin);
        double read = double(duration_read.count());
        double transform = double(duration_transform.count());
        double total = double(duration.count());

        cout<<"----name:"<<mname<<"----"<<endl;
        cout<<"----time to read:" <<read<<"----"<<endl;
        cout<<"--------time for bitmap:"<<time_bitmap<<"----"<<endl;
        cout<<"--------time for v8:"<<time_v8<<"----"<<endl; // v8+serial
        cout<<"--------time for serial:"<<time_serial<<"----"<<endl;
        cout<<"--------time for wbsort:"<<time_wbsort<<"----"<<endl;
        cout<<"----sum up:"<<time_bitmap+time_v8+time_wbsort<<"----"<<endl;
        cout<<"----time for entire transformation:"<<transform<<"----"<<endl;
        cout<<"----total time used:" <<total<<"----"<<endl;
        
        // cout<<"done"<<endl;//test
        // if(sseq.empty()) cout<<"empty";//test

        // else{
        // for(int i = 0;i<sseq.size();i++){
        //     cout<<sseq[i]<<" ";
        // }}
    }
    //system("pause");
    
    return 0;
    }