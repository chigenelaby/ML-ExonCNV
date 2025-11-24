#-*-coding:utf-8-*-

'''
输出指定样本的核型信息和性别
用法示例：
python xxx.py aneuploid.txt output_karyotype.txt 
'''

__author__ = 'housy'

import os
import re
import fire
import sys

X_info_FM = {"loss1": "X", "loss2": "", "gain1": "XXX", "gain2": "XXXX", "NA": "XX"}
X_info_M = {'loss': "", "gain": 'XX', "NA": "X"}
Y_info = {"loss": "", "gain": "YY", "NA": "Y"}

class IO:
    def load_file(self, f):
        res = []
        with open(f, "r") as fo:
            for i in fo:
                if not i.startswith("#"):
                    res.append(i.split('\n')[0])
        return res
    def write_file(self, content, f):
        with open(f, "w") as fw:
            for i in content:
                print("\t".join(i), file = fw)

def get_karyotype_info(sex_qc_data):
    '''得到染色体数目,核型信息
    res type: ["chr1_1_249250621_wild_wild_C" "chr1" 1 249250621	0.9872818997454369 "NA"]
    输出列表，两个元素，例如[45,XY,-21  XY]
    '''
    dic = {}
    for line in sex_qc_data:
        if line.split("\t")[1] not in dic.keys():
            dic[line.split("\t")[1]] = line.split("\t")[5]
        else:
            print(line)
    if "chrY" in dic.keys() and "chrX" in dic.keys():
        karyotype = X_info_M[dic["chrX"]] + Y_info[dic["chrY"]]
    elif "chrY" not in dic.keys() and "chrX" in dic.keys():
        karyotype = X_info_FM[dic["chrX"]]
    else:
        karyotype = "No Sex chromosomes"
        
    normal_num = 44
    ind = []
    for i in dic.keys():
        if i != "chrX" and i != "chrY":
            if re.findall(r'\d+', i) and re.match(r'loss', dic[i]):
                ind.append("-" + re.findall(r'\d+', i)[0])
            elif re.findall(r'\d+', i) and re.match(r'gain', dic[i]):
                ind.append("+" + re.findall(r'\d+', i)[0])
            if dic[i] == "loss2":
                normal_num -= 2
            elif dic[i] == "loss1":
                normal_num -= 1
            elif dic[i] == "gain1":
                normal_num += 1
            elif dic[i] == "gain2":
                normal_num += 2
            else:
                continue
    normal_num = normal_num + len(karyotype)
    index = ",".join(ind)
    if index == '':
        res = [str(normal_num) + "," + karyotype, karyotype]
    else:
        res = [str(normal_num) + "," + karyotype + ',' + index, karyotype]
    return res

def get_karyotype(sex_qc_file, output_file):
    '''得到该样本的核型信息
    参数: sex_qc_file, output_file
    '''   
    obj = IO()
    sex_qc_data = obj.load_file(sex_qc_file)
    res = get_karyotype_info(sex_qc_data)
    obj.write_file([res], output_file)

if __name__ == "__main__":
    fire.Fire(get_karyotype)
