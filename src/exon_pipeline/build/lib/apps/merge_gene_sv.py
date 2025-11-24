#-*-coding:utf-8-*-
"""
sv数据 整合wes exon
1. 以sv结果为基础 如果有相交区段则标注相交证据
2. 补充wes exons项目单独的数据
3. 整合家系信息
    - 如果wgs和wes一致，则对应列数，整合为wgs文库号
    - 如果wgs和wes数据不一致，不对应的列数单独成列，0/0=0 wild 相对可靠
"""
import re
from itertools import *
from mod_tools.mod_data_IO import *
from intervaltools.core import IntervalBase, BucketIndexIntervalList
from exon_pipeline.utils.conv_tools import conv_filter
from exon_pipeline.apps.merge_gene_file import *

reg = re.compile('[0-9A-Z][A-Z]{2}\d{3}')

def change_line(indexs_items):
    """返回line的改变函数
    indexs_items raw_indexs, res_indexs 将原始数据的index转换成结果数据索引
    """
    def func(line):
        line = list(line)
        tmp = list(line)
        for x, y in zip(*indexs_items):
            line[y] = tmp[x]
        return line
    return func

def line_add_info(infos_map, wkcode_index_dict, n):
    count = 0
    for i, j in infos_map:
        if i in wkcode_index_dict:
            pass
        elif j in wkcode_index_dict:
            if not i:
                pass
            else:
                v = wkcode_index_dict.pop(j)
                wkcode_index_dict[i] = v
        else:
            count += 1
            i = i if i else j
            wkcode_index_dict[i] = [n, n+1, n+2]
            n += 3
    return count 

def get_add_data_func(count, raw=None):
    """返回给line增加count组数据"""
    if not raw:
        raw = ['0/0=1', 'wild', '相对可靠']
    def func(line):
        line = list(line)
        for i in range(count):
            line.extend(raw)
        return line
    return func

def get_header(header, h_dict):
    """生成新的header"""
    h = []
    junk = list(chain(*h_dict.values()))
    for index, i in enumerate(header):
        if index in junk:
            break
        h.append(i)
    
    for wk, v in sorted(h_dict.items(), key=lambda x:x[1]):
        hs = ['%s_exon'%wk, '%s_tag'%wk, '%s_result'%wk]
        h.extend(hs)
    return h

def get_is_deep_tag_func(h_dict1):
    def is_deep_tag(line):
        for wk, indexs in h_dict1.items():
            if line[indexs[1]][:4] in {'loss', 'gain'}:
                return True
        return False
    return is_deep_tag


def merge_gene_sv(sv_gene_file, wes_gene_file=None, wgs_map_file=None):
    """
    """
    # 读取文件
    headers1, res1 = load_tsv(sv_gene_file)
    headers2, res2 = load_tsv(wes_gene_file)
    infos_map = load_file(wgs_map_file, '\t', '#')
    
    # 预处理
    h_dict1 = get_wkcode_index_from_headers(headers1)
    h_dict2 = get_wkcode_index_from_headers(headers2)

    infos_map = [(i, j) for i, j in infos_map if i in h_dict1 or j in h_dict2]
    conv1 = BucketIndexIntervalList(res1, key=IntervalLine, default_base=IntervalLine)
    conv2 = BucketIndexIntervalList(res2, key=IntervalLine, default_base=IntervalLine)
    conv1.make_index()
    conv2.make_index()
    add_index2 = []
    c1 = line_add_info(infos_map, h_dict1, len(headers1))
    c2 = line_add_info(infos_map, h_dict2, len(headers2))
    # 填充成员结果: sv报出的结果,如缺少家系成员的检测结果,但在wes结果中已存在(例如单人wgs+家系wes),则补充wes家系成员的检测结果;否则填充未覆盖;反之同理-wangxf-20230602
    for i, j in infos_map:
        if not i and j:
            add_index2.extend(h_dict2[j])
    for i in conv1:
        raw = []
        conv_com = conv2.find_one(i)
        if add_index2 and conv_com:
            for index in add_index2:
                raw.append(conv_com.line[index])
        else:
            raw = ['0/0=0', 'Uncov', 'NA']*c1
        tmp_i = list(i.line)
        tmp_i.extend(raw)
        i.line = tuple(tmp_i)
    add_index1 = []
    for i, j in infos_map:
        if not j and i:
            add_index1.extend(h_dict1[j])
    for i in conv2:
        raw = []
        conv_com = conv1.find_one(i)
        if add_index1 and conv_com:
            for index in add_index1:
                raw.append(conv_com.line[index])
        else:
            raw = ['0/0=0', 'Uncov', 'NA']*c2
        tmp_i = list(i.line)
        tmp_i.extend(raw)
        i.line = tuple(tmp_i)

    # 过滤
    is_deep_tag = get_is_deep_tag_func(h_dict1)
    conv_filter(conv2, key=lambda x: not conv1.find_one(x, action='intersect') or not is_deep_tag(x.line))
    
    #c1 = line_add_info(infos_map, h_dict1, len(headers1))
    #c2 = line_add_info(infos_map, h_dict2, len(headers2))
    	
    
    #func1 = get_add_data_func(c1)
    #func2 = get_add_data_func(c2)
    
    #conv1.add_attr('line', key=lambda x:func1(x.line))
    #conv2.add_attr('line', key=lambda x:func2(x.line))
    
    # 转换位置
    infos = [i if i else j for i, j in infos_map]
    change1 = []
    change2 = []
    for i in infos:
        change1.extend(h_dict1[i])
        change2.extend(h_dict2[i])
    
    #func_change = change_line((change2, change1))
    func_change = change_line((change1, change2))
    conv2.add_attr('line', key=lambda x:func_change(x.line))

    merge_conv(conv1, conv2)
    res = [i.line for i in conv1]
    
    headers = get_header(headers1, h_dict1)
    
    return headers, res
