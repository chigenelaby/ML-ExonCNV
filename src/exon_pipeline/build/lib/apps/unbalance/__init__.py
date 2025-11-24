#-*-coding:utf-8-*-
"""不均衡系数"""
from operator import itemgetter
from mod_tools.mod_data_IO import IO_pickle
from exon_pipeline.apps.reliability import *
from exon_pipeline.utils.database_store import get_index
import exon_pipeline.parameters as opt
import collections

UNBALANCE_LIMIT1 = opt.UNBALANCE_LIMIT1
UNBALANCE_LIMIT2 = opt.UNBALANCE_LIMIT2


# def whether_roughovlap(a, b):
#     bool1 = a.chr == b.chr and a.start < b.end and b.start < a.end
#     bool2 = False
#     if bool1:
#         start_o = max(a.start, b.start)
#         end_o = min(a.end, b.end)
#         if end_o - start_o + 1 > 50:
#             bool2 = True

#     return bool2


# def whether_ovlap(a, b):
#     tmpbool = a.chr == b.chr and a.start <= b.end and b.start <= a.end

#     return tmpbool


def group_chr(complex):
    d = collections.defaultdict(list)
    for c in complex:
        d[c.chr].append(c)
    return d


def calculate_unbalance(filter_convs, junk_exons):
    """输入样本的初步过滤的convs和junk_exons, 返回不均衡系数"""

    # total_unbalreg = 100179
    # count = 0

    # dict_chr_junkexons = group_chr(junk_exons)

    # for conv in filter_convs:
    #     for ex in conv.exons:
    #         for junk_e in dict_chr_junkexons[ex.chr]:
    #             if whether_ovlap(junk_e, ex):
    #                 count += 1
    #                 break

    # print("Runnning calculate unbalance {} {} {}".format(count, total_unbalreg, count/total_unbalreg))
    
    # return count / total_unbalreg

    total_lenght = len(junk_exons)
    count = 0

    # dict_chr_convs = group_chr(filter_convs)
    
    # for junk_e in junk_exons:
    #     for conv in filter_convs:
    #         for ex in conv.exons:
    #             if junk_exons.find_one(ex):
    #                 count += 1

    # for junk_e in junk_exons:
    #     print("junk_e", junk_e.chr, junk_e.start, junk_e.end)
    #     b = False
    #     for conv in dict_chr_convs[junk_e.chr]:
    #         for ex in conv.exons:
    #             # print("ex", ex.chr, ex.start, ex.end)
    #             if junk_e.chr == ex.chr and junk_e.start == ex.start and junk_e.end == ex.end:
    #                 count += 1
    #                 b = True
    #                 break

    #         if b:
    #             break

    for conv in filter_convs:
        for ex in conv.exons:
            if junk_exons.find_one(ex):
                # print("ex {} {} {}".format(ex.chr, ex.start, ex.end))
                count += 1

    # print(total_lenght, count, count/total_lenght)

        # b = False
        
        # for conv in dict_chr_convs[junk_e.chr]:
        #     # print(">>>", conv.chr)
        #     for ex in conv.exons:
        #         # print("ex", ex.chr, ex.start, ex.end)
        #         if whether_roughovlap(junk_e, ex):
        #             b = True
        #             break            
        #     if b:
        #         break

        # if b:
        #     count += 1
    
    # i = 0
    # for exon in junk_exons:
    #     print("junk", exon.chr, exon.start, exon.end)
    #     if i>5:
    #         break
    #     i += 1
    return count / total_lenght


def get_unbalance_string(unbalance_value):
    """返回不均衡文档"""
    if unbalance_value < UNBALANCE_LIMIT1:
        return """样本正常"""
    elif unbalance_value < UNBALANCE_LIMIT2:
        return """轻微不均衡"""
    else:
        return """严重不均衡"""

def get_severity_unbalance_tag(unbalance_value):
    """返回严重不均衡判断字段"""
    if unbalance_value >= UNBALANCE_LIMIT2:
        return 'NO'
    else:
        return 'YES'

def fmt_unbalance_item(unbalance_value):
    """
    不均衡质控格式化
    degradation coefficient(~0.03~0.06~), 样本正常|YES(0.1071)
    degradation coefficient(~0.03~0.06~), 样本正常|YES(0.10712)
    """
    unbalance_string = get_unbalance_string(unbalance_value)
    severity_unbalance_tag = get_severity_unbalance_tag(unbalance_value)
    header = 'degradation coefficient(~%s~%s~)'%(UNBALANCE_LIMIT1, UNBALANCE_LIMIT2)
    # info = '%s|%s(%.4f)'%(unbalance_string, severity_unbalance_tag, unbalance_value)
    info = '%s|%s(%.5f)'%(unbalance_string, severity_unbalance_tag, unbalance_value)
    return (header, info)

def get_unbalance_item(filter_convs_file, junk_exons_file):
    """
    filter_convs_file 过滤后conv文件
    junk_exons_file 不均衡标注exons文件
    """
    convs = IO_pickle.load_pickle(filter_convs_file)
    junk_exons = get_index(junk_exons_file, key=itemgetter(1,2,3))
    unbalance_value = calculate_unbalance(convs, junk_exons)
    unbalance_item = fmt_unbalance_item(unbalance_value)
    return unbalance_item
