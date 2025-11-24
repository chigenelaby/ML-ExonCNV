#-*-coding:utf-8-*-
"""可靠性模块"""
import os
import logging
import re
from operator import *

import numpy as np

from mod_tools.mod_data_IO import *
from intervaltools.core import *
from exon_pipeline.core.contrast import *
from exon_pipeline.core.controlled import *
from exon_pipeline.core.conv import *
from exon_pipeline.utils.io import fmt_write_file, show
from exon_pipeline.utils.flags import FlagInfoBase, FlagScoreBase

from exon_pipeline.apps.reliability.delta_diff import diff_condition
from exon_pipeline.apps.reliability.exon_reliability import *
from exon_pipeline.apps.reliability.flag_family import *
from exon_pipeline.apps.vaf import *
from exon_pipeline.apps.reliability.flag_cnv_seq import *

logger = logging.getLogger(__name__)
##flag_score

flag_diff_True = FlagScoreBase('flag_diff', classify=True, score=1, show_string="满足diff值条件")
flag_diff_None = FlagScoreBase('flag_diff', classify=None, score=0, show_string="")
flag_diff_False = FlagScoreBase('flag_diff', classify=False, score=-2, show_string="")

flag_is_cnv_True = FlagScoreBase('flag_is_cnv', classify=True, score=0, show_string="")
flag_is_cnv_False = FlagScoreBase('flag_is_cnv', classify=False, score=0, show_string="")

flag_length_True = FlagScoreBase('flag_length', classify=True, score=1, show_string="满足外显子个数条件")
flag_length_False = FlagScoreBase('flag_length', classify=False, score=0, show_string="")

flag_reliable_exon_True = FlagScoreBase('flag_reliable_exon', classify=True, score=1, show_string="满足可靠外显子比例条件")
flag_reliable_exon_False = FlagScoreBase('flag_reliable_exon', classify=False, score=0, show_string="")

flag_vaf_True = FlagScoreBase('flag_vaf', classify=True, score=1, show_string="vaf支持")
flag_vaf_None = FlagScoreBase('flag_vaf', classify=None, score=0, show_string="")
flag_vaf_False = FlagScoreBase('flag_vaf', classify=False, score=-2, show_string="vaf不支持")

flag_cnv_family_True = FlagScoreBase('flag_cnv_family', classify=True, score=1, show_string="存在至少一个家系成员支持")
flag_cnv_family_None = FlagScoreBase('flag_cnv_family', classify=None, score=0, show_string="")

flag_cnvseq_True = FlagScoreBase('flag_cnvseq', classify=2, score=2, show_string="wgs支持")
flag_cnvseq_oppo2 = FlagScoreBase('flag_cnvseq', classify=-2, score=-2, show_string="wgs结论相反")
flag_cnvseq_oppo1 = FlagScoreBase('flag_cnvseq', classify=-1, score=-1, show_string="wgs弱不支持")
flag_cnvseq_None = FlagScoreBase('flag_cnvseq', classify=0, score=0, show_string="")
####
reliable_exon_ratio_dict1 = {2: lambda x: True, 1: lambda x: x >= 0.6}
reliable_exon_ratio_dict2 = {2: lambda x: x > 0.3, 1: lambda x: x >= 0.7}

reliable_exon_condition1 = ExonReliabilityCondition(ratio_dict=reliable_exon_ratio_dict1)
reliable_exon_condition2 = ExonReliabilityCondition(ratio_dict=reliable_exon_ratio_dict2)


class GetFlagDiffCondition:
    flag_type = 'flag_diff'
    """diff值flag"""
    def get_flag(self, conv):
        tag = conv.tag
        diff = conv.diff
        if tag == 'NA':
            return flag_diff_False
        res = diff_condition.forward(tag, diff)
        if res >=1 :
            return flag_diff_True
        elif res < 0:
            return flag_diff_False
        else:
            return flag_diff_None

class GetFlagIsCnv:
    flag_type = 'flag_is_cnv'
    """cnv标签"""
    def get_flag(self, conv):
        """cnv标签"""
        tag = conv.tag
        diff = conv.diff
        if tag == 'NA':
            return flag_is_cnv_False
        res = diff_condition.forward(tag, diff)
        if res <0:
            return flag_is_cnv_False
        length = len(conv.exons)
        flag = flag_is_cnv_False
        if tag == 'loss2' and length >= 5:
            flag = flag_is_cnv_False
        elif length >= 7:
            flag = flag_is_cnv_True
        return flag


class GetFlagLengthCondition:
    flag_type = 'flag_length'
    """获取长度标签"""
    def __init__(self, reliable_exon_condition=reliable_exon_condition1):
        self.reliable_exon_condition = reliable_exon_condition
    
    def get_flag(self, conv):
        """长度标签需要满足reliable_exon_condition1"""
        if not self.reliable_exon_condition.forward(conv):
            return flag_length_False
        tag = conv.tag
        length = len(conv.exons)
        flag = flag_length_False
        if tag == 'loss2' and length >= 3:
            flag = flag_length_True
        elif tag == 'loss1' and length >= 5:
            flag = flag_length_True
        elif length >= 10:
            flag = flag_length_True
        return flag


class GetReliableExonCondition:
    """获取长度标签"""
    flag_type = 'flag_reliable_exon'
    def __init__(self, reliable_exon_condition=reliable_exon_condition2):
        self.reliable_exon_condition = reliable_exon_condition
    
    def get_flag(self, conv):
        """长度标签需要满足reliable_exon_condition1"""
        if self.reliable_exon_condition.forward(conv):
            return flag_reliable_exon_True
        else:
            return flag_reliable_exon_False

class GetCNVFamilyCondition(FlagCNVFamily):
    """输入family_path1, family_path2"""
    flag_type = 'flag_cnv_family'
    def get_flag(self, conv):
        """
        flag True 家系相交且tag同
        flag None 无家系数据
        flag False 家系数据无支持数据
        """
        length = len(conv.exons)
        if length < 2:  ## 20211009 by housy
            flag = None
        else:
            flag = super().get_flag(conv)
        if flag:
            return flag_cnv_family_True
        else:
            return flag_cnv_family_None


class GetVafCondition:
    flag_type = 'flag_vaf'
    def __init__(self, vcf_path=None, vaf_db_save_path=None, repeat_path=None):
        self.vcf_path = vcf_path
        self.vaf_db_save_path = vaf_db_save_path
        self.repeat_path = repeat_path
        self.vaf_obj = None
        self.initialize()
    
    def initialize(self):
        """初始化"""
        try:
            vcf = fmt_load_vcf_obj(self.vcf_path)
            ds = DatabaseStore()
            ds.load(self.vaf_db_save_path)
            repeat_obj = IO_pickle.load_pickle(self.repeat_path)
            vaf_obj = FlagExonVafBase(repeat_obj=repeat_obj, vcf_obj=vcf, vcf_db=ds)
        except:
            raise
            logger.info('获取vaf数据error')
            vaf_obj = None
        self.vaf_obj = vaf_obj
    
    def get_flag(self, conv):
        """获取vaf标签"""
        if self.vaf_obj is None:
            return flag_vaf_None
        flag = self.vaf_obj.get_flag(conv)
        if flag.classify[:1] == '1':
            return flag_vaf_None
        elif flag.classify[:1] == '2':
            return flag_vaf_False
        elif flag.classify[:1] == '3':
            return flag_vaf_True
        return flag_vaf_None


class GetCnvSeqCondition:
    flag_type = 'flag_cnvseq'
    def __init__(self, wgs_res_path=None, wgs_medi_path=None):
        self.wgs_res_path = wgs_res_path
        self.wgs_medi_path = wgs_medi_path
        self.cnv_seq_condition_obj = None
        self.initialize()
        self.flag_conv_score_dict = {'1.1.0': flag_cnvseq_oppo2, 
                                     '1.1.1': flag_cnvseq_True, 
                                     '1.1.2': flag_cnvseq_oppo1, 
                                     '1.2.1.0': flag_cnvseq_oppo2, 
                                     '1.2.1.1': flag_cnvseq_True, 
                                     '1.2.2.0': flag_cnvseq_oppo1, 
                                     '1.2.2.1': flag_cnvseq_None, 
                                     '1.3.0': flag_cnvseq_None, 
                                     '1.3.1': flag_cnvseq_None, 
                                     '1.9': flag_cnvseq_None, 
                                     }
    
    def initialize(self):
        try:
            fs = GetFlagCnvSeq(self.wgs_res_path, self.wgs_medi_path)
        except:
            logger.info('获取cnvseq数据error')
            fs = None
        self.cnv_seq_condition_obj = fs
    
    def get_flag(self, conv):
        """获取校正标签"""
        obj = self.cnv_seq_condition_obj
        if obj is None:
            return flag_cnvseq_None
        flag = obj.get_flag(conv)
        return self.flag_conv_score_dict[flag.classify]


class StatisticsFlagBase:
    """设置和统计flag条件"""
    def __init__(self, exons_conv):
        self.exons_conv = exons_conv
        self.conditions = []
        self.flags = []
    
    def add_flag_condition(self, condition):
        """添加flag条件
        flag_condition 应当有get_flag接口
        """
        self.conditions.append(condition)
        self.flags.append(condition.flag_type)
        logger.info('添加条件 %s'%condition.flag_type)
    
    def get_one_score(self, conv):
        """计算分"""
        score = 0
        for flag_type in self.flags:
            flag = getattr(conv, flag_type)
            score += flag.score
        return score
    
    def forward(self):
        """运行，计算条件"""
        exons_conv = self.exons_conv
        for condition in self.conditions:
            flag_type = condition.flag_type
            func = getattr(condition, 'get_flag')
            exons_conv.add_attr(flag_type, key=func)
        exons_conv.add_attr('score', key=self.get_one_score)
        

