import numpy as np

def gen_seq_txt(seq, size,path=''):
    FileName ='seq'
    if path:
        path += '/'
    sseq = str(seq).replace(', ', ' ').replace('[', '').replace(']', '')
    f = open('{}.txt'.format(path + FileName), 'w')
    f.write(sseq)
    f.close()
    f = open('{}/info.txt'.format(path), 'a')
    f.writelines(str(size)+'\n')
    f.close()

def SeqReverse(seq):
    length = len(seq)
    reseq = np.zeros(length)
    for i in range(length):
        reseq[seq[i]] = i
    return reseq

def gen_rseq(seq_bitmap, seq_v8):
    seq_row = seq_bitmap[seq_v8]#得到2次行重排最终的行重排顺序
    # print("seq_row =", seq_row)
    return seq_row

def gen_wseq(seq_row, seq_dict):
    seq_wb = []
    for key in seq_dict.keys():
        seq_wb.append(seq_row[key])
    # print(seq_wb)
    wseq = np.array(seq_wb)
    return wseq

def SerialSort(seq_bitmap, seq_v8, seq_dict, path=''):
    seq_row = gen_rseq(seq_bitmap, seq_v8)
    # print(seq_row)
    rseq=SeqReverse(seq_row)
    # print(rseq)
    seq_wb = list(gen_wseq(rseq, seq_dict).astype(int))
    gen_seq_txt(seq_wb,len(seq_wb),path)
    return seq_wb

def SerialSort_block(seq_bitmap, seq_v8, seq_list, path=''):
    wseq_list = []
    seq_row = gen_rseq(seq_bitmap, seq_v8)
    rseq=SeqReverse(seq_row)
    for seq in seq_list:
        wseq = list(gen_wseq(rseq, seq).astype(int))
        wseq_list.append(wseq)
    len_seq=0
    for i in range(len(wseq_list)):
        len_seq=len_seq+len(wseq_list[i])
        # print(len(wseq_list[i]))
    gen_seq_txt(wseq_list, len_seq, path)
    return wseq_list

# def main():
#     seq_bitmap = np.array([0,1,2,3])
#     seq_v8 = np.array([2,0,3,1])
#     seq_dict = {1:0,2:1,0:2,3:3}
#     seq_dict1 = {0:1, 1:2, 2:0, 3:4, 4:3}
#     seq_list = [seq_dict, seq_dict1]
#     seq_wb = SerialSort(seq_bitmap, seq_v8, seq_dict)
#     # wseq_list = SerialSort_block(seq_bitmap, seq_v8, seq_list)
#     print(seq_wb)
#     # print(wseq_list)

# main()
