#-*-coding:utf-8-*-
"""非整倍体验证"""

import os
import logging
import fire
from operator import itemgetter
from itertools import repeat
import numpy as np
from exon_pipeline.core.conv import *


def get_tolomere_tag_flag(tolomere_conv, convs):
    """相同tag比例大于60%"""
    if tolomere_conv.tag == 'NA':
        return False
    match_convs = convs.find_all(tolomere_conv, action='intersect')
    counts = [exon for i in match_convs if i.tag[:4]==tolomere_conv.tag[:4] for exon in i.exons if tolomere_conv.include(exon)]
    ratio = len(counts) / len(tolomere_conv.exons)
    # print(tolomere_conv.chr, tolomere_conv.start, tolomere_conv.end, tolomere_conv.tag, ratio, len(counts), len(tolomere_conv.exons))
    #if tolomere_conv.chr == "chr13":
    #    print(ratio, len(counts), len(tolomere_conv.exons), tolomere_conv.tag)
    # if tolomere_conv.tag == 'NA':
    #    return False
    # print(ratio, tolomere_conv, convs)

    thld = 0.85 if tolomere_conv.chr != "chrY" else 0.68
    # return ratio > 0.85
    return ratio > thld

def get_mean_delte_diff_flag(conv):
    """获取非正倍体平均delte_diff flag标签
    非性染色体。mean_delte_diff 要小于0.15    
    """
    if conv.chr in {'chrX', 'chrY'}:
        return True
    delte_diff = [abs(i.diff - j.diff) for i, j in zip(conv.exons, conv.exons[1:])]
    mean_delte_diff = np.mean(delte_diff)
    # print(conv.chr, conv.start, conv.end, mean_delte_diff)
    # if conv.chr == "chr13":
    #     print(mean_delte_diff)
    # print(mean_delte_diff)

    # return mean_delte_diff < 0.2
    return mean_delte_diff < 0.25

