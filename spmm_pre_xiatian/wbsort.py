import numpy as np

def gen_seq_txt(seq, type='', path=''):
    FileName = 'Sequence' + type
    if path:
        path += '/'
    sseq = str(seq).replace(', ', ' ').replace('[', '').replace(']', '')
    f = open('{}.txt'.format(path + FileName), 'w')
    f.write(sseq)
    f.close()

def SeqReverse(seq):
    length = len(seq)
    reseq = np.zeros(length)
    for i in range(length):
        reseq[seq[i]] = i
    return reseq

def gen_rseq(seq_bitmap, seq_v8):
    seq_row = seq_bitmap[seq_v8]#得到2次行重排最终的行重排顺序
    #print("seq_row =", seq_row)
    return SeqReverse(seq_row)

def gen_wseq(seq_row, seq_dict):
    seq_wb = []
    for key in seq_row:
        if key in seq_dict.keys():
            seq_wb.append(seq_dict[key])#依据seq_dict做列重排写回顺序
        else:
            seq_wb.append(seq_row.shape[0] - 1)#直接写回到末尾
    wseq = SeqReverse(seq_wb)
    return wseq

def SerialSort(seq_bitmap, seq_v8, seq_dict, path=''):
    seq_row = gen_rseq(seq_bitmap, seq_v8)
    seq_wb = list(gen_wseq(seq_row, seq_dict).astype(int))
    gen_seq_txt(seq_wb,path=path)
    return seq_wb

def SerialSort_block(seq_bitmap, seq_v8, seq_list, path=''):
    wseq_list = []
    seq_row = gen_rseq(seq_bitmap, seq_v8)
    for seq in seq_list:
        wseq = list(gen_wseq(seq_row, seq).astype(int))
        wseq_list.append(wseq)
    gen_seq_txt(wseq_list, 'Block', path)
    return wseq_list

def gen_ioseq(seq_bitmap, seq_v8, seq_dict, path=''):
    size = seq_bitmap.shape[0]
    seq = np.array(list(range(size)))
    iseq = list(gen_wseq(seq, seq_dict).astype(int))
    oseq = list(SeqReverse(gen_rseq(seq_bitmap, seq_v8)).astype(int))
    gen_seq_txt(iseq, 'Input', path)
    gen_seq_txt(oseq, 'Output', path)
    return iseq, oseq

def gen_ioseq_block(seq_bitmap, seq_v8, seq_list, path=''):
    size = seq_bitmap.shape[0]
    seq = np.array(list(range(size)))
    iseq_list = []
    for s in seq_list:
        iseq = list(gen_wseq(seq, s).astype(int))
        iseq_list.append(iseq)
    gen_seq_txt(iseq_list, 'BlockInput', path)
    oseq = list(SeqReverse(gen_rseq(seq_bitmap, seq_v8)).astype(int))
    gen_seq_txt(oseq, 'BlockOutput', path)
    return iseq_list, oseq

def main():
    seq_bitmap = np.array([1, 2, 0, 4, 3])
    seq_v8 = np.array([2, 0, 3, 4, 1])
    seq_dict = {0:1, 1:2, 2:0, 3:4, 4:3}
    seq_dict1 = {0:2, 1:0, 2:3, 4:1}
    seq_list = [seq_dict, seq_dict1]
    seq_wb = SerialSort(seq_bitmap, seq_v8, seq_dict)
    wseq_list = SerialSort_block(seq_bitmap, seq_v8, seq_list)
    #print(seq_wb)
    gen_ioseq(seq_bitmap, seq_v8, seq_dict)
    gen_ioseq_block(seq_bitmap, seq_v8, seq_list)
    #print(wseq_list)

if __name__ == '__main__':
    main()