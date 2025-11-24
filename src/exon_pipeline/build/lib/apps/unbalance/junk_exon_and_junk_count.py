#-*-coding:utf-8-*-

'''
功能1：给定一批样本(训练集)，得到这批样本的junk exon数据集
功能2：给定另一批样本(测试集)，得到这批样本阳性结果在junk exon数据集的比例

#step1 输入文库号，批次号, 得到阳性结果
#step2 统计阳性结果的个数
#step3 得到junk exon结果
#step4 统计某个样本得到的阳性结果出现在junk exon中的比例

输入文件至少有三列，前三列分别是批次号，文库号，文件路径(flag_info_filter.txt)

author: housy
'''

import glob
import re
import fire
import random
import subprocess
from operator import *
from intervaltools.core import *
#from mod_tools import *

class IntervalTag(IntervalBase):
    def __init__(self, chr, start, end, tag):
        super().__init__(chr, start, end)
        self.tag = tag

class IO:
    def load_file(self, f):
        res = []
        with open(f, "r") as fo:
            for i in fo:
                if not i.startswith("#"):
                    res.append(i.split('\n')[0].split('\t'))
        return res
    def write_file(self, content, f):
        with open(f, "w") as fw:
            for i in content:
                i = list(map(str, i))
                print("\t".join(i), file = fw)

def load_data(data, key=itemgetter(0, 1, 2)):
    #data = load_file(path, '\t', '#')
    idx = BucketIndexIntervalList(data, default_base=IntervalBase, key=key)
    idx.make_index(10000)
    return idx

#key1 = itemgetter(0, 1, 2, 11)
#key2 = itemgetter(0, 1, 2)

#def load_data2(data, key=key1):
    #data = load_file(path, '\t', '#')
    #idx = BucketIndexIntervalList(data, default_base=IntervalTag, key=key)
    #idx.make_index(10000)
    #return idx

def get_intersect(idx1, idx2 ):
    """文件1和文件2区间取交集，输出文件1存在交集的部分"""
    res = []
    for ex in idx1:
        if ex.chr != "chrX" and ex.chr != "chrY":
            intersect_intervals = idx2.find_all(ex, action='intersect')
            if intersect_intervals:
                res.append([*list(ex.interval())])
    res = list(map(list, list(set([tuple(t) for t in res]))))
    return res

def get_report_exon(batch, wkcode, path, exon_bed_file):
    '''得到一个文件的阳性外显子结果'''
    final = ''
    if glob.glob(path):
        res_path = glob.glob(path)[0]
        idx2 = []
        obj = IO()
        idx2 = obj.load_file(res_path)
        idx1 = obj.load_file(exon_bed_file)
        idx1 = load_data(idx1)
        idx2 = load_data(idx2)
        res = get_intersect(idx1, idx2)
        final = [[*i, batch, wkcode] for i in res]
    else:
        print("FileNotFound:" + batch + "_" + wkcode + " 没有阳性结果文件")
    return final

def get_all_report_exon(batch_wkcode_file, exon_bed_file, output_file = ''):
    '''得到输入所有文件的阳性结果
    '''
    obj = IO()
    batch_wkcode_data = obj.load_file(batch_wkcode_file)
    all_report_exon_lis = []
    #p = ProgressBar(len(batch_wkcode_data))
    for i in batch_wkcode_data:
        batch = i[0]
        wkcode = i[1]
        path = i[2]
        res = get_report_exon(batch, wkcode, path, exon_bed_file)
        for ex in res:
            all_report_exon_lis.append([*ex])
        #p.read()
    if output_file == "":
        return all_report_exon_lis
    else:
        obj.write_file(all_report_exon_lis, output_file)
        return output_file

def get_all_report_exon_count(all_report_exon_lis, exon_bed_file):
    '''输入训练数据集中所有文库得到的阳性外显子，二维列表
    输出每个外显子出现的次数
    '''
    all_report_exon_lis = [[i[0], i[1], i[2]] for i in all_report_exon_lis]
    dic = {}
    obj = IO()
    exon_bed = obj.load_file(exon_bed_file)
    #p = ProgressBar(len(all_report_exon_lis))
    for i in all_report_exon_lis:
        i = list(map(str, i))
        if "\t".join(i) not in dic.keys():
            dic["\t".join(i)] = 1
        else:
            dic["\t".join(i)] += 1
        #p.read()
    for i in exon_bed:
        key = '\t'.join(i[:3])
        if key not in dic:
            dic[key] = 0
    return dic

def get_all_report_exon_count2(all_report_exon_lis_file, exon_bed_file):
    '''输入训练数据集中所有文库得到的阳性外显子，二维列表
    输出每个外显子出现的次数
    '''
    info = subprocess.Popen(['less ' + all_report_exon_lis_file + '|cut -f1-3 | sort | uniq -c | sed "s/^[ \t]*//g" | sed "s/ /\t/g" '], shell = True, stdout = subprocess.PIPE)
    array_info = info.stdout.read().decode('utf-8').split('\n')[:-1]
    obj = IO()
    exon_bed = obj.load_file(exon_bed_file)
    dic = {}
    #p = ProgressBar(len(array_info))
    for i in array_info:
        key = "\t".join(i.split('\t')[1:])
        dic[key] = int(i.split('\t')[0])
        #p.read()
    for i in exon_bed:
        key = '\t'.join(i[:3])
        if key not in dic:
            dic[key] = 0
    return dic

def get_junk_exon(all_report_exon_lis, exon_bed_file, count = 50000):
    '''从1/2所在的值开始向上取值，直到取到50000个，例如如果在20时超过50000，则在20时随机选取部分'''
    count = int(count)
    if isinstance(all_report_exon_lis, str):
        dic = get_all_report_exon_count2(all_report_exon_lis, exon_bed_file)
    else:
        dic = get_all_report_exon_count(all_report_exon_lis, exon_bed_file)
    value = []
    for i in dic:
        value.append(int(dic[i]))
    value = sorted(value)
    start = value[int(len(value)/2)] + 1
    train_dic = {}
    #p = ProgressBar(len(dic.keys()))
    for i in dic.keys():
        if dic[i] in train_dic.keys():
            train_dic[dic[i]].append([*i.split('\t')])
        else:
            train_dic[dic[i]] = [[*i.split('\t')]]
        #p.read()
    res = []
    while len(res) < count and start in train_dic.keys():
        if len(res) >= count or start not in train_dic.keys():
            break
        if len(res) < count and len(train_dic[start]) <= (count - len(res)):
            res.extend([str(start), *ex] for ex in train_dic[start])
        elif len(res) < count and len(train_dic[start]) > (count - len(res)):
            residual = count - len(res)
            random.shuffle(train_dic[start])
            res.extend([str(start), *ex] for ex in train_dic[start][:residual])
        start += 1
    if len(res) != count:
        print("Count value is Error")
    return res

def get_junk_count(junk_exon, report_exon, key=itemgetter(0, 1, 2)):
    '''得到外显子在junk exon中的个数'''
    idx1 = BucketIndexIntervalList(junk_exon, default_base=IntervalBase, key=key)
    idx1.make_index(10000)
    idx2 = BucketIndexIntervalList(report_exon, default_base=IntervalBase, key=key)
    idx2.make_index(10000)
    count = 0
    for ex in idx1:
        intervals = idx2.find_all(ex, action='__eq__')
        if intervals:
            count += 1
    return count

def main1(batch_wkcode_file, output_file, count=50000):
    '''输入文库批次文件，得到junk exon结果
    '''
    exon_bed_file = "/mnt/rhd/EXdev_chen/WES_pipe_v4_chen/db/K_cap_bed/NT01T_target.bed"
    all_report_exon_lis = get_all_report_exon(batch_wkcode_file, exon_bed_file, output_file + "_report_exon")
    junk_exon = get_junk_exon(all_report_exon_lis, exon_bed_file)
    IO().write_file(junk_exon, output_file)

def main2(batch_wkcode_file, junk_exon_file, output_file):
    '''输入文件：文库批次文件，junk exon文件，得到每个样本在junk_exon数据集下的个数
    '''
    exon_bed_file = "/mnt/rhd/EXdev_chen/WES_pipe_v4_chen/db/K_cap_bed/NT01T_target.bed"
    final = []
    obj = IO()
    sample_data = obj.load_file(batch_wkcode_file) 
    junk_exon_data = obj.load_file(junk_exon_file)
    junk_exon = [[i[1], i[2], i[3]] for i in junk_exon_data]
    all_junk = len(junk_exon)
    #p = ProgressBar(len(sample_data))
    for i in sample_data:
        batch = i[0]
        wkcode = i[1]
        path = i[2]
        report_exon = get_report_exon(batch, wkcode, path, exon_bed_file)
        report_exon_count = len(report_exon)
        count1 = get_junk_count(junk_exon, report_exon, key=itemgetter(0, 1, 2))
        frac = str(round(int(count1)/all_junk, 4))
        final.append([*i, report_exon_count, count1, frac])
        #p.read()
    IO().write_file(final, output_file)

if __name__ == '__main__':
    fire.Fire({ "junk_exon": main1,
            "junk_count": main2})

