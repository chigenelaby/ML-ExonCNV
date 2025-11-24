#-*-coding:utf-8-*-
"""可靠性评级"""

import os
import logging
from operator import attrgetter


logger = logging.getLogger(__name__)


def get_is_high_freq(freq_conv, junk_intervals=None):
    """高频"""
    if junk_intervals and junk_intervals.find_one(freq_conv, action='intersect'):
        return False
    if freq_conv.freq > 0.05 or freq_conv.all_freq > 0.15:
        return True
    else:
        return False

def figure_score(conv, flags):
    """计算分"""
    score = 0
    for flag_type in flags:
        flag = getattr(conv, flag_type)
        score += flag.score
    return score

def is_contradiction(conv, flags):
    """存在矛盾证据"""
    res_score = [getattr(conv, flag_type).score for flag_type in flags]
    tag1 = [s for s in res_score if s>0]
    tag2 = [s for s in res_score if s<0]
    return len(tag1) and len(tag2)

def fmt_reliability_value(all_value, min_limit, max_limit):
    """格式化分数"""
    if all_value < min_limit:
        all_value = min_limit
    if all_value > max_limit:
        all_value = max_limit
    return all_value

def get_reliability_show(freq_conv, seq=';', flags=None,):
    if flags is None:
        # flags = ['flag_cnv_family', 'flag_vaf', 'flag_cnvseq', 
        #           'flag_diff', 'flag_is_cnv', 'flag_length', 'flag_reliable_exon']
        flags = ['flag_diff', 'flag_is_cnv', 'flag_length', 'flag_reliable_exon']
    return list(filter(lambda x:x, map(lambda x:getattr(x, 'show_string'), attrgetter(*flags)(freq_conv))))

def get_reliability_value(freq_conv, junk_intervals=None):
    # print(0)
    """输入是包含所有标签的conv"""
    max_limit = 3
    min_limit = 0
    all_flag = {i for i in freq_conv.__dict__ if i.startswith('flag')}
    other_flags = {'flag_cnv_family', 'flag_vaf', 'flag_cnvseq'}
    base_flags = {'flag_diff', 'flag_is_cnv', 'flag_length', 'flag_reliable_exon'}
    is_high_freq = get_is_high_freq(freq_conv, junk_intervals=junk_intervals)
    if is_high_freq:
        max_limit = 1
    base_value = figure_score(freq_conv, base_flags)
    other_value = figure_score(freq_conv, other_flags)
    all_value = base_value + other_value
    
    # 返回未知 
    # 1.矛盾数据
    # 2.非loss2 和 cnv 低频 无辅助证据
    # 3. cnv value>=3 高频
    if is_contradiction(freq_conv, other_flags) and all_value >= 2:
        return None
    if not freq_conv.flag_is_cnv.classify and freq_conv.tag != 'loss2':
        if not is_high_freq and not other_value:
            #return None
            if all_value >= 2:
                return None
            else:
                pass
    if freq_conv.flag_is_cnv.classify and is_high_freq and all_value >= 3:
        return None
    
    if len(freq_conv.exons) == 1:
        if freq_conv.tag[:4] == 'loss':
            all_value += [0, -1][is_high_freq]  ## 20211009 by housy
        else:
            all_value += [0, -1][is_high_freq]
    elif len(freq_conv.exons) <= 3:
        all_value += [0, -1][is_high_freq]
    else:
        all_value += [0, -2][is_high_freq]
    
    # loss2 低频 无矛盾证据 
    if freq_conv.tag == 'loss2':
        if not is_high_freq and not is_contradiction(freq_conv, other_flags) and all_value >= 2:
            return 3
    
    level = fmt_reliability_value(all_value, min_limit, max_limit)
    
    # 区段长度大于10M
    if freq_conv.__len__() >= 10_000_000 and level is not None and level < 3:
        level = 3
    
    return level

def transform_value(value):
    """转化可靠性"""
    d = {None: "未知",
         0: "不可靠",
         1: "可能不可靠",
         2: "相对可靠",
         3: "特别可靠"}
    return d[value]

def get_reliability(freq_conv, junk_intervals=None):
    value = get_reliability_value(freq_conv, junk_intervals=junk_intervals)
    return transform_value(value)
