#-*-coding:utf-8-*-

"""
频率库
"""

import os
import copy
import logging
import itertools
import functools
from operator import itemgetter
from collections import OrderedDict
from exon_pipeline.core.conv import *
from exon_pipeline.utils.io import fmt_write_file, show
from exon_pipeline.apps.reliability import *


class FreqBase(DatabaseStore):
    """
    频率库
        e.g.
        >>> save_path = '/share/chg2master/prod/Other/tanggb/task/task20200701_exon_rebuild/exon_src/test/freq_data/XX_freq.pic'
        >>> fq = FreqBase()
        >>> fg.load(save_path)
        >>> fq.get_freq('chr1', 69090, 622034)
        OrderedDict([('loss1', 0.05555555555555555), ('loss2', 0.0), ('gain1', 0.05555555555555555), ('gain2', 0.0)])
        >>> fq.get_freq_count('chr1', 69090, 622034)
        OrderedDict([('loss1', 1), ('loss2', 0), ('gain1', 1), ('gain2', 0)])
        mode=2：表示外来的数据，计算频率时判断存在交集的区段占比>=0.9且非交集区段<1000
        mode=1: 表示正常WES数据，计算频率时判断包含关系
    """
    def get_freq_infos(self, chr, start, end, mode=1):
        """获取区间tag详情信息"""
        tag_exon = IntervalBase(chr, start, end)
        freq_intervals = []
        for intervals in self._store:
            if str(mode) == "2":
                interval = intervals.find_one(tag_exon, action='check_intersect', ratio_key=(None, 0.9))
                if interval:
                    no_intersect_len = tag_exon.end - tag_exon.start - (min(interval.end, tag_exon.end) - max(interval.start, tag_exon.start))
                    if no_intersect_len > 5000:
                        interval = None
            else:
                interval = intervals.find_one(tag_exon, action='include')
            if interval:
                freq_intervals.append(interval)
        
        freq_st_dict = {}
        for interval in freq_intervals:
            freq_st_dict.setdefault(interval.tag, []).append(interval)
        return freq_st_dict
    
    def get_freq_count(self, chr, start, end, mode):
        """获取频次信息"""
        freq_st_dict = self.get_freq_infos(chr, start, end, mode)
        
        freq_tags = ['loss1', 'loss2', 'gain1', 'gain2']
        freq_dict = OrderedDict()
        for tag in freq_tags:
            freq = len(freq_st_dict.get(tag, []))
            freq_dict[tag] = freq
        return freq_dict
    
    def get_freq(self, chr, start, end, mode=1):
        """获取频率信息"""
        total = len(self._store)
        freq_dict = self.get_freq_count(chr, start, end, mode)
        freq_dict1 = OrderedDict()
        for tag, fq in freq_dict.items():
            freq_dict1[tag] = fq / total
        return freq_dict1
