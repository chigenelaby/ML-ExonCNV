#-*-coding:utf-8-*-
"""对照库使用，定义了对照库系数"""

import os
import re
import logging

from intervaltools.core import *

CONTRAST_DEPTH_FILTER = 30
CONTRAST_CV_FILTER = 0.3
CONTRAST_SIZE_FILTER = 10

class ControlledBase:
    """对照类
    输入的应该是对照库contrast BucketIndexIntervalList(IntervalContrast)的实例
    """
    def __init__(self, contrast):
        self.contrast = contrast
        contrast.initialized = False
    
    def get_flag_bad_capture(self, exon):
        """输入单个对照库exons信息
        返回捕获情况bad_capture标签
        """
        if exon.dps_mean < CONTRAST_DEPTH_FILTER:
            return True
        if exon.size < CONTRAST_SIZE_FILTER:
            return True
        return False
    
    def get_flag_bad_contrast(self, exon):
        """输入单个对照库exons信息
        返回捕获情况bad_contrast
        """
        if self.get_flag_bad_capture(exon):
            return True
        if exon.cv > CONTRAST_CV_FILTER:
            return True
        return False
    
    def get_diff_weight(self, exon):
        """获取diff系数 即该对照库depth/diff 应该乘的权重"""
        try:
            weight = 1 / exon.mean
        except ZeroDivisionError:
            weight = 1
        return weight
    
    def forwark(self, reload=False):
        """将对照库加上flag标签和diff_coefficient"""
        contrast = self.contrast
        if not reload and getattr(contrast, 'initialized', None):
            return None
        contrast.add_attr('weight', key=self.get_diff_weight)
        contrast.add_attr('flag_bad_capture', key=self.get_flag_bad_capture)
        contrast.add_attr('flag_bad_contrast', key=self.get_flag_bad_contrast)
        contrast.initialized = True
    
    def get_diff(self, exons):
        """输入一个样本exons情况
        exons BucketIndexIntervalList(IntervalDepth) 
        添加 属性
        diff
        flag_bad_capture
        flag_bad_contrast
        sample_value
        z_value
        """
        self.forwark()
        contrast = self.contrast
        exons.add_attr('contrast', key=contrast.find_one)
        exons.add_attr('sample_value', key=lambda x: x.depth / exons.sample_depth)
        exons.add_attr('diff', key=lambda x: x.sample_value * x.contrast.weight)
        exons.add_attr('flag_outlier', key=lambda x: (x.sample_value < x.contrast.min_critical or
                                                      x.sample_value > x.contrast.max_critical))
        return exons

