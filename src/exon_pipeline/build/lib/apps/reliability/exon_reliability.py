#-*-coding:utf-8-*-
"""单个外显子可靠性"""
from itertools import *
from exon_pipeline.apps.reliability.delta_diff import diff_condition

class ExonReliabilityCondition:
    """单个外显子可靠性"""
    def __init__(self, diff_condition=diff_condition, ratio_dict=None):
        self.diff_condition = diff_condition
        self.ratio_dict = ratio_dict
    
    def forward_one(self, tag, exon):
        """计算外显子可靠性值"""
        res = self.diff_condition.forward(tag, exon.diff)
        cv = exon.contrast.cv
        flag_outlier = exon.flag_outlier
        flag_bad_contrast = exon.contrast.flag_bad_contrast
        if res == -2:
            return 0
        if flag_bad_contrast:
            return 1
        if res > 0 and flag_outlier:
            return 2
        return 1
    
    def _reliable_exon_ratio_format(self, ratio_dict=None, levels=[0,1,2]):
        if ratio_dict is None:
            ratio_dict = dict(zip(levels, repeat(None)))
        for key, value in ratio_dict.items():
            if value is None:
                ratio_dict[key] = lambda x: True
        return ratio_dict
    
    def forward(self, conv, ratio_dict=None):
        """单个外显子可靠性比例条件
        ratio_dict是一个等级: 函数的字典
            等级是单个外显子可靠性等级。
            函数是参数为符合该等级(大于等于)的比例.
        """
        tag = conv.tag
        if tag == 'NA':
            return False
        if ratio_dict is None:
            ratio_dict = self.ratio_dict
        exon_num = len(conv.exons)
        ratio_dict = self._reliable_exon_ratio_format(ratio_dict)
        flags = []
        for key, v_func in ratio_dict.items():
            levels = []
            for ex in conv.exons:
                level = self.forward_one(tag, ex)
                levels.append(level)
            
            ratio = len([level for level in levels if level >= key]) / exon_num
            flags.append(v_func(ratio))
        return all(flags)