#include<iostream>
#include<vector>
#include<string>
#include<fstream>
#include<map>
#include <chrono>
#include<unordered_map>
#include <omp.h>

//#define CORENUM 8

#include "csr.h"
#include "bitmap.h"
#include "transmat.h"
#include "v8sort.h"
#include "wbsort.h"

#define MAXLENGTH 1000000
#define SECT 2048
#define MAXREGION 20

using namespace std;
using namespace chrono;
double time_bitmap=0;
double time_v8=0;
double time_serial = 0;
double time_wbsort = 0;
double time_b = 0;
double time_a = 0;
double time_c = 0;
double time_d = 0;
double time_e = 0;
double time_f = 0;
double time_g = 0;
double time_gpnl=0;
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

// void gen_serial_origin(SpM &mr, SpM &new_mr, unordered_map<int,int> &seq_dict, vector<int> &seq_order) {
//     int seq_cnt = 0;
//     int n_cols = mr.shape[1];
//     int n_vals = mr.shape[2];

//     bool bitmap[n_cols] = {0};
//     int* new_seq = new int[n_vals];

//     seq_dict.reserve(n_cols);
//     seq_dict.max_load_factor(0.5);

//     for (int i = 0; i < n_vals; i++) {
//         int ind_iter = mr.indices[i];

//         if (bitmap[ind_iter] == 0) {
//             // auto start1 = chrono::high_resolution_clock::now();
//             bitmap[ind_iter] = 1;
//             new_seq[i] = seq_cnt;
//             seq_dict.emplace(ind_iter, seq_cnt);
//             seq_cnt++;
//             // auto end1 = chrono::high_resolution_clock::now();
//             // auto duration = chrono::duration_cast<chrono::microseconds>(end1 - start1);
//             // time1 += duration.count();
//         } else {
//             // auto start2 = chrono::high_resolution_clock::now();
//             new_seq[i] = seq_dict[ind_iter];
//             // auto end2 = chrono::high_resolution_clock::now();
//             // auto duration = chrono::duration_cast<chrono::microseconds>(end2 - start2);
//             // time2 += duration.count();
//         }
//         // if(bitmap[ind_iter] == 0) seq_dict.emplace(ind_iter, seq_cnt++);
//         // bitmap[ind_iter] = 1;
//         // new_seq[i] =  seq_dict[ind_iter];
//     }

//     new_mr = SpM(mr.data, new_seq, mr.indptr, mr.shape);
//     delete[] new_seq;
// }

void gen_serial_origin_vec(int *shape,int *indices,int *new_seq,int* &seq_key,int &end_seq_key){
    int seq_cnt = 0;
    vector<int> bitmap(shape[1],-1);
    // new_seq = new int[mr.shape[2]];
    for(int i = 0;i<shape[2];i++){
        int ind_iter = indices[i];
        int bit_value = bitmap[ind_iter];
        if(bit_value == -1){
            bitmap[ind_iter] = seq_cnt;
            seq_key[end_seq_key++]=ind_iter;
            new_seq[i] = seq_cnt;
            seq_cnt += 1;
        }
        else{
            new_seq[i] = bit_value;
        }
    }
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

void gen_trace_formats(SpM &mr,int* &seq_input,int &end_seq_input,vector<int> &rseq,int &bnum,int **bserial_colidx,double **bserial_data,int **bserial_indptr,int *shape,int *bsize_list,int *seq_offset){

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

    // vector<unordered_map<int,int>> bseq_list;
    int* bseq_list_key=new int[mr.shape[2]];
    int end_bseq_list_key=0;
    int indptr_offset = 0;
    int panel_size = 2048;

    int regions_length = 0;
    bsize_list[0]=0;
    seq_offset[0]=0;
    int end_bsize_list=1;
    // cout<<"enter gen_new_panels"<<endl;
    // cout << mr_bitmap.shape[0] << ' ' << mr_bitmap.shape[1] << ' ' << mr_bitmap.shape[2] << endl;
    auto begin_gpnl= high_resolution_clock::now();
    // regions_length = gen_new_panels(mr_bitmap,regions,bsize_list,end_bsize_list,bnum);//gen_panel
    regions_length = gen_new_panels(mr_bitmap,bsize_list,end_bsize_list,bnum);
    auto end_gpnl= high_resolution_clock::now();
    auto duration_gpnl = duration_cast<microseconds>(end_gpnl - begin_gpnl);
    time_gpnl +=  double(duration_gpnl.count());
    // cout << "done" << endl;
    // cout<<"out gen_new_panels"<<endl;
    // regions[0].check();
    //cout<<regions_length<<endl;
    int cnt = 0;
    vector<vector<int>> spv8_lists;
    // int* seq_v8 =new int[mr.shape[0]];
    // int end_seq_v8=0;
    vector<int> spv8_list;
    vector<int> add_panelsize_list;
    // cout << regions_length << endl;
    int end_add_panelsize_list = 0;
    end_seq_input=0;

    int** seq_v8_block=new int*[regions_length];
    for(int i=0;i<regions_length;i++) seq_v8_block[i] = new int[mr.shape[0]];
    int *end_seq_v8_block=new int[regions_length]();

    bserial_data=new double*[regions_length];
    for(int i=0;i<regions_length;i++) bserial_data[i] = new double[mr.indptr[mr.shape[0]]];
    int *end_bserial_data=new int[regions_length]();

    bserial_indptr = new int*[regions_length];
    int *end_bserial_indptr= new int[regions_length];
    for(int i=0;i<regions_length;i++) {
        bserial_indptr[i] = new int[mr.shape[0]+1];
        bserial_indptr[i][0]=0;
        end_bserial_indptr[i] = 1;
    }

    bserial_colidx = new int*[regions_length];
    for(int i=0;i<regions_length;i++) bserial_colidx[i] = new int[mr.indptr[mr.shape[0]]];
    int *end_bserial_colidx=new int[regions_length];

    // cout<<regions_length<<endl;
    #pragma omp parallel for num_threads(CORENUM)
    for(int index = 0;index < regions_length;index++){
        //vector<int>().swap(seq_v8);
        // end_seq_v8=0;
        vector<int>().swap(spv8_list);
        int* add_panelsize_list=new int[mr.indptr[mr.shape[0]]];
        end_add_panelsize_list=0;

        auto begin_a= high_resolution_clock::now();
        // gen_panel_list(regions[index], add_panelsize_list, end_add_panelsize_list);
        gen_panel_list(mr_bitmap, bsize_list[index], bsize_list[index+1], add_panelsize_list,end_add_panelsize_list);//a
        auto end_a= high_resolution_clock::now();
        auto duration_a = duration_cast<microseconds>(end_a - begin_a);
        time_a +=  double(duration_a.count());

        auto begin_b= high_resolution_clock::now();

        // panel_sort_nnz(seq_v8,end_seq_v8,spv8_list,regions[index],add_panelsize_list,end_add_panelsize_list);//b
        panel_sort_nnz(seq_v8_block[index], end_seq_v8_block[index], spv8_list, mr_bitmap, bsize_list[index], bsize_list[index+1], add_panelsize_list, end_add_panelsize_list);
        auto end_b= high_resolution_clock::now();
        auto duration_b = duration_cast<microseconds>(end_b - begin_b);
        time_b +=  double(duration_b.count());
        // int *seq_v8_arr = new int[seq_v8.size()];
        // for(int i = 0;i < seq_v8.size();i++) seq_v8_arr[i] = seq_v8[i];

        auto begin_c= high_resolution_clock::now();
        end_bserial_indptr[index] = bsize_list[index+1] - bsize_list[index] + 1;
        //SpM tmp_r = reorder_row(mr_bitmap, bsize_list[index], bsize_list[index+1], seq_v8_block[index],bserial_indptr[index]);//c
        
        auto indptr = mr_bitmap.indptr + bsize_list[index];
        int offset = indptr[0];
        auto data = mr_bitmap.data + indptr[0];
        auto indices = mr_bitmap.indices + indptr[0];
        int shape[3];
        shape[0] = bsize_list[index+1] - bsize_list[index];
        shape[1] = mr_bitmap.shape[1];
        shape[2] = indptr[shape[0]] - indptr[0];

        SpM tmp_r(shape);
        int tail = 0;
        for(int i = 0; i<shape[0]; i++){
            int s = seq_v8_block[index][i];
            int ind = (tmp_r.indptr[i]+(indptr[s+1]-indptr[s]));
            tmp_r.indptr[i+1] = ind;
            bserial_indptr[index][i+1] = ind;
            for( int j = indptr[s]-offset; j < indptr[s+1]-offset; j++ ){
                tmp_r.data[tail] = data[j];
                tmp_r.indices[tail] = indices[j];
                tail++;
            }
        }
        int *rcolidx=tmp_r.indices;
        end_bserial_data[index] = shape[2];
        bserial_data[index]=tmp_r.data;
        int bserial_data_tail = 0;
        int base = 0;
        for(int i = 0; i < (int)(spv8_list.size()); i++){
            int count = spv8_list[i];
            int panelsize = add_panelsize_list[i+1] - add_panelsize_list[i];
            for(int row_index = base; row_index < base + count; row_index += 8){
                int rowlen = tmp_r.indptr[row_index + 1] - tmp_r.indptr[row_index];
                int nnz_index = tmp_r.indptr[row_index];
                double *add_data = new double[rowlen*8];
                int *add_colidx = new int[rowlen*8];
                for(int row = 0; row < rowlen; row++)
                    for(int col = 0; col < 8; col++){
                        add_data[8 * row + col] = data[nnz_index + col * rowlen + row];
                        add_colidx[8 * row + col] = tmp_r.indices[nnz_index + col * rowlen + row];}
                for(int j = 0; j < rowlen*8; j++){
                    bserial_data[index][bserial_data_tail+j] = add_data[j];
                    rcolidx[bserial_data_tail+j] = add_colidx[j];
                }
                bserial_data_tail+=rowlen*8;
                delete[] add_data;
                delete[] add_colidx;
            }
            if(i == (int)(spv8_list.size()-1)){
                int row_begin = base + count;
                int nnz_begin = tmp_r.indptr[row_begin];
                bserial_data_tail+=shape[2]-nnz_begin;
            }
            else{
                int row_begin = base + count;
                int row_end = base + panelsize;
                int nnz_begin = tmp_r.indptr[row_begin];
                int nnz_end = tmp_r.indptr[row_end];
                bserial_data_tail+=nnz_end-nnz_begin;
            }
            base += panelsize;
        }
        
        auto end_c= high_resolution_clock::now();
        auto duration_c = duration_cast<microseconds>(end_c - begin_c);
        time_c +=  double(duration_c.count());

        //int *rcolidx = new int[tmp_r.shape[2]];
        // double *rdata = new double[tmp_r.shape[2]];
        //end_bserial_data[index] = tmp_r.shape[2];

        // auto begin_d= high_resolution_clock::now();
        // transpose_spv8_nnz(tmp_r,spv8_list,add_panelsize_list,rcolidx,bserial_data[index]);//d
        // auto end_d= high_resolution_clock::now();
        // auto duration_d = duration_cast<microseconds>(end_d - begin_d);
        // time_d +=  double(duration_d.count());
      
        
        delete[] add_panelsize_list;add_panelsize_list=NULL;

        auto begin_f= high_resolution_clock::now();
        // cout<<end_seq_v8<<" "<<tmp_r.shape[2]<<" "<<tmp_r.shape[0] - 1<<endl;
        // for(int i=0;i<end_seq_v8;i++) seq_v8_block[end_seq_v8_block++]=seq_v8[i] + bsize_list[cnt];//f
        // for(int i = 0;i < tmp_r.shape[2];i++) bserial_data[end_bserial_data++]=rdata[i];
        // for(int i = 1;i < tmp_r.shape[0] - 1;i++) bserial_indptr[end_bserial_indptr++]=tmp_r.indptr[i]+indptr_offset;
        auto end_f= high_resolution_clock::now();
        auto duration_f = duration_cast<microseconds>(end_f - begin_f);
        time_f +=  double(duration_f.count());
        
        indptr_offset = indptr_offset + tmp_r.indptr[tmp_r.shape[0]];
        spv8_lists.push_back(spv8_list);
        cnt = cnt + 1;
        
        auto begin_serial= high_resolution_clock::now();

        // int *seq = new int[tmp_r.shape[2]];
        //vector<int> seq_key;
        // unordered_map<int,int> seq;

        // gen_serial_origin(regions[index],smr,seq,order);
        end_bserial_colidx[index] = tmp_r.shape[2];
        gen_serial_origin_vec(tmp_r.shape,rcolidx,bserial_colidx[index],bseq_list_key,end_bseq_list_key);
        //bseq_list_key.push_back(seq_key);

        auto end_serial= high_resolution_clock::now();
        auto duration_serial = duration_cast<microseconds>(end_serial - begin_serial);
        time_serial += double(duration_serial.count());

        

        // for(int cc = 0;cc < tmp_r.shape[2]; cc++) bserial_colidx[end_bserial_colidx++]=seq[cc] + offset;
        //依赖项
        seq_offset[index + 1] = end_bseq_list_key;
        // delete[] seq;
        
    }
    // delete[] seq_v8;seq_v8=NULL;
    //cout<<mr.shape[0]<<" "<<mr.shape[1]<<" "<<mr.shape[2]<<" "<<end_bseq_list_key<<endl;
    auto end_v8= high_resolution_clock::now();
    auto duration_v8 = duration_cast<microseconds>(end_v8 - begin_v8);
    time_v8 = double(duration_v8.count());
    // cout << "end for" << endl;
    auto begin_wbsort = high_resolution_clock::now();

    auto begin_g= high_resolution_clock::now();

    shape[0] = mr.shape[0];
    shape[1] = mr.shape[1];
    shape[2] = mr.shape[2];
    // smr_out.check();
    //这里是要一维的数组
    int seq_len = 0;
    for(int i=0;i<regions_length;i++) seq_len += end_seq_v8_block[i];
    SerialSort_block_vec(seq_bitmap,seq_v8_block,end_seq_v8_block,bseq_list_key,end_bseq_list_key,seq_input,end_seq_input,bsize_list,regions_length,seq_len);
    auto end_g= high_resolution_clock::now();
    auto duration_g = duration_cast<microseconds>(end_g - begin_g);
    time_g =  double(duration_g.count());
    // seq_input.resize(0);
    // for(auto i = seq_temp.begin();i != seq_temp.end();i++){
    //     for(auto j = i->begin();j != i->end();j++){
    //         seq_input.push_back(*j);
    //     }
    // }
    int* seq_row = gen_rseq(seq_bitmap,seq_v8_block,end_seq_v8_block,bsize_list,regions_length,seq_len);
    int* rseq_arr= SeqReverse(seq_row, seq_len);
    delete[] seq_row; seq_row = NULL;
    rseq = vector<int>(rseq_arr, rseq_arr + seq_len);
    
    delete[] rseq_arr; rseq_arr = NULL;
    // sseq.resize(0);
    // for(int i=0;i<seq_key_list.size();i++){
    //     for(auto k : seq_key_list[i]){
    //         sseq.push_back(k);
    //     }
    // }
    //delete[] seq_v8_block; seq_v8_block=NULL;
    delete[] seq_bitmap; seq_bitmap = NULL;
    auto end_wbsort = high_resolution_clock::now();
    auto duration_wbsort = duration_cast<microseconds>(end_wbsort - begin_wbsort);
    time_wbsort += double(duration_wbsort.count());

}

int main(){
    vector<string> mlist;
    ifstream file("matrix.txt");
    string s;
    while(getline(file,s)){
        mlist.push_back(s);
    }
    file.close();
    
    ofstream outfile;
    string filepath = "result.txt";
    outfile.open(filepath);
    if(!outfile){
        cout<<"result.txt not exist";
        exit(1);
    }
    for(auto mlist_iter = mlist.begin();mlist_iter != mlist.end();mlist_iter++){
        auto begin = high_resolution_clock::now();
        auto begin_read = high_resolution_clock::now();
        string mname = *mlist_iter;
        mname = get_mname(mname);
        string path1 = "mat/mtx/";
        string path2 = "/";
        string path3 = ".mtx";
        string mpath = path1+mname+path2+mname+path3;
        SpM mr = csr_matrix(mpath);
        auto end_read = high_resolution_clock::now();
        auto duration_read = duration_cast<microseconds>(end_read - begin_read);
        
        int* seq=new int[mr.shape[2]];
        int end_seq=0;
        vector<int> rseq;
        int bnum;
        int **bserial_colidx;
        double **bserial_data;
        int **bserial_indptr;
        int *shape = new int[3];
        //to be test
        int *bsize_list = new int[1000];
        int *seq_offset = new int[1000];
        
        auto begin_transform = high_resolution_clock::now();
        gen_trace_formats(mr,seq,end_seq,rseq,bnum,bserial_colidx,bserial_data,bserial_indptr,shape,bsize_list,seq_offset);
        auto end_transform = high_resolution_clock::now();
        auto end = high_resolution_clock::now();
        auto duration_transform = duration_cast<microseconds>(end_transform - begin_transform);
        auto duration = duration_cast<microseconds>(end- begin);
        double read = double(duration_read.count());
        double transform = double(duration_transform.count());
        double total = double(duration.count());

        outfile<<mname<<" "<<(time_bitmap+time_v8+time_wbsort)/1000<<"ms"<<endl;
        cout<<"******************"<<endl;
        cout<<"----name:"<<mname<<"----"<<endl;
        cout<<"----time to read:" <<read<<"----"<<endl;
        cout<<"--------time for bitmap:"<<time_bitmap<<"----"<<endl;
        cout<<"--------time for v8:"<<time_v8<<"----"<<endl; // v8+serial
        cout<<"--------time for serial:"<<time_serial<<"----"<<endl;
        cout<<"--------time for wbsort:"<<time_wbsort<<"----"<<endl;
        cout<<"----sum up:"<<time_bitmap+time_v8+time_wbsort<<"----"<<endl;
        cout<<"----time for entire transformation:"<<transform<<"----"<<endl;
        cout<<"----total time used:" <<total<<"----"<<endl;
        cout<<"--------time for gen_new_panels:"<<time_gpnl<<"----"<<endl;
        cout<<"--------time for a:"<<time_a<<"----"<<endl;
        cout<<"--------time for b:"<<time_b<<"----"<<endl;
        cout<<"--------time for c:"<<time_c<<"----"<<endl;
        cout<<"--------time for d:"<<time_d<<"----"<<endl;
        cout<<"--------time for e:"<<time_e<<"----"<<endl;
        cout<<"--------time for f:"<<time_f<<"----"<<endl;
        cout<<"--------time for g:"<<time_g<<"----"<<endl;
        cout<<"----total for v8:"<<time_a+time_b+time_c+time_d+time_e+time_f+time_gpnl<<"----"<<endl;


        // cout<<"done"<<endl;//test
        // if(sseq.empty()) cout<<"empty";//test

        // else{
        // for(int i = 0;i<sseq.size();i++){
        //     cout<<sseq[i]<<" ";
        // }}
    }
    outfile.close();
    //system("pause");
    
    return 0;
    }
