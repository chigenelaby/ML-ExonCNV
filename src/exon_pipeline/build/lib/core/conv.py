#-*-coding:utf-8-*-
"""连接"""

import os
import re
import logging
from operator import attrgetter

import numpy as np

from intervaltools.core import *
from exon_pipeline.apps.ExonConv import *
from exon_pipeline.utils.io import *
from exon_pipeline.utils.conv_tools import get_conv_diff

class FlagLoss2Condition(StaticCondition1):
    """遇到明确的loss2, 非loss2将中断"""
    def __init__(self, exs):
        self.samples = exs
    
    def condition(self, exon, tag=None):
        """条件判断
        type exon: Exon实例 
        rtype: bool
        """
        flag = True
        if tag == 'loss2':
            return flag
        try:
            chr = exon.chr
            start = exon.start
            end = exon.end
            ex = self.samples.find_one(IntervalBase(chr, start, end))
            if ex.flag_loss2:
                flag = False
        except:
            pass
        return flag

def is_male(sexuality='XX'):
    """是否是男性"""
    return sexuality == 'XY'

def set_conv(exons, cnvs, sexuality):
    """将外显子配置到cnvs区段中，并配置相关属性"""
    flag_is_male = is_male(sexuality)
    
    cnvs.sexuality = sexuality
    cnvs.flag_is_male = flag_is_male
    
    cnvs.add_attr('exons', key=lambda x: exons.find_all(x, action='be_included'))
    # cnvs.check_attr('exons', errors='strict')
    cnvs.add_attr('diff', key=get_conv_diff)
    # cnvs.check_attr('diff', errors='strict')
    cnvs.add_attr('is_male_sex_chr', key=lambda x:(flag_is_male and is_sex_chr(x.chr)))
    cnvs.add_attr('tag', key=lambda x:get_tag_from_diff(x.diff, x.is_male_sex_chr, get_is_autosomal_region(x)))
    cnvs.check_attr('tag', errors='strict')
    
    # 信息继承
    cnvs.wkcode = exons.wkcode


def get_conv_data(exons, base=IntervalBase, sexuality='XX'):
    """外显子连接
    """
    r = show(exons, 'chr start end diff contrast.flag_bad_capture'.split())
    r = [(*i[:-1], ['NA', 'bad_capture'][i[-1]]) for i in r]
    
    conv_obj = ExonConvCNV(r, key=None, sexuality=sexuality)
    
    #染色体号条件
    chr_conditon = ChrCondition()
    conv_obj.add_condition(chr_conditon, condition_type='static_end')
    
    #外显子个数条件
    ex_num_conditon = ExonNumCondition(0.5)
    conv_obj.add_condition(ex_num_conditon, condition_type='dynamic')
    
    #diff值条件
    dff_conditon = DiffConvCondition()
    conv_obj.add_condition(dff_conditon, condition_type='static_conv')
    
    # 内部cnv条件，内部不应该出现准确的cnv，静止终止条件
    neg_num_conditon = NegativeNumCondition()
    conv_obj.add_condition(neg_num_conditon, condition_type='static_end')
    
    # 否定个数条件。连续10个非tag外显子即中断.
    inner_cnv_conditon = InnerCnvCondition()
    conv_obj.add_condition(inner_cnv_conditon, condition_type='static_end')
    
    # 外显子间距条件：外显子间距大于500k时，单侧长度必须大于间距0.33或单侧外显子个数不小于10连接
    exon_spacing_conditon = ExonSpacingCnvCondition()
    conv_obj.add_condition(exon_spacing_conditon, condition_type='static_conv')
    
    # loss2强制判断
    flag_loss2_conditon = FlagLoss2Condition(exons)
    conv_obj.add_condition(flag_loss2_conditon, condition_type='static_end')
    
    cnvs_res = conv_obj.run()
    
    cnvs_inv = [(i[0].chr, i[0].start, i[-1].end) for i in cnvs_res]
    
    cnvs = BucketIndexIntervalList(cnvs_inv, default_base=IntervalBase)
    cnvs.make_index()
    
    #将外显子配置到cnvs区段
    set_conv(exons, cnvs, sexuality)
    
    return cnvs
