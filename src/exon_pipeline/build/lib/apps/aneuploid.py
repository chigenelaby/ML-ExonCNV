#-*-coding:utf-8-*-
"""非整倍体"""

import os
import logging
import copy
from operator import itemgetter
from itertools import repeat
from intervaltools.core import IntervalBase, BucketIndexIntervalList
from exon_pipeline.utils.database_store import get_index
from exon_pipeline.utils.conv_tools import conv_filter
from exon_pipeline.utils.sample_tools import fmt_id_string
from exon_pipeline.utils.fmt_tools import get_fmt_tag
from exon_pipeline.core.conv import set_conv
from exon_pipeline.apps.aneuploid_check import get_tolomere_tag_flag, get_mean_delte_diff_flag

logger = logging.getLogger(__name__)


class TolomereConv:
    """端粒合并"""
    def __init__(self, tolomere_file, exons, convs):
        self.tolomere_file = tolomere_file
        self.exons = exons
        self.convs = convs
        self.sexuality = self.convs.sexuality
        self.show_dict = {True: '非整倍体验证通过', False: '非整倍体验证不通过'}
        self.tolomere_convs = self.get_tolomere_conv()
        self.aneuploid = self.get_aneuploid_qc()
    
    def get_tolomere_conv(self,):
        """构建端粒conv"""
        tolomere_file = self.tolomere_file
        exons = self.exons
        sexuality = self.sexuality
        tolomere_convs = get_index(tolomere_file, key=itemgetter(0,1,2))
        set_conv(exons, tolomere_convs, sexuality)
        return tolomere_convs

    def is_reliable_tolomere(self, tolomere_conv):
        """是可靠的非整倍体
        1. 满足tolomere_tag_flag 标签
        2. 满足mean_delte_diff_flag 标签    
        """
        # print("tolomere_conv", tolomere_conv, tolomere_conv.tag)

        if tolomere_conv.tag == "NA":
            return False
        convs = self.convs
        tolomere_tag_flag = get_tolomere_tag_flag(tolomere_conv, convs)
        mean_delte_diff_flag = get_mean_delte_diff_flag(tolomere_conv)
        # if tolomere_conv.chr == "chr21":
        #     print(tolomere_conv.start, tolomere_conv.end, tolomere_tag_flag, mean_delte_diff_flag)
        #if tolomere_conv.tag == "NA":
        #    return False
        # print("tolomere_tag_flag", tolomere_tag_flag, "mean_delte_diff_flag", mean_delte_diff_flag)
        
        return tolomere_tag_flag and mean_delte_diff_flag
    
    def get_aneuploid_qc(self, ):
        """获取非整倍体质控"""
        exons = self.exons
        sexuality = self.sexuality
        tolomere_convs = self.tolomere_convs
        res = {}
        for conv in tolomere_convs:
            res.setdefault(conv.chr, []).append(conv)
        
        c = []
        for chr, conv_list in res.items():
            if len(conv_list) == 2:
                conv_list = sorted(conv_list)
                c1, c2 = conv_list
                merge_conv = IntervalBase(chr, c1.start, c2.end)
                c.append(merge_conv)
                continue
            c.extend(conv_list)
        c = BucketIndexIntervalList(sorted(c))
        c.make_index()
        set_conv(exons, c, sexuality)
        if sexuality == 'XX':
            conv_filter(c, lambda x:x.chr!='chrY')
        c.add_attr('is_reliable_tolomere', key=self.is_reliable_tolomere)
        c.add_attr('fmt_tag', key=lambda x:get_fmt_tag(x) if x.is_reliable_tolomere else 'NA')
        c.add_attr('id_string', key=lambda x: fmt_id_string(x.chr, x.start, x.end, x.fmt_tag, x.fmt_tag))
        show_dict = self.show_dict
        c.add_attr('show_string', key=lambda x: show_dict[x.is_reliable_tolomere])
        return c
    
    
    def merge_tolomere(self, tolomere_convs):
        """端粒conv合并"""
        res = {}
        for conv in tolomere_convs:
            res.setdefault(conv.chr, []).append(conv)
        
        c = []
        for chr, conv_list in res.items():
            if len(conv_list) == 2:
                conv_list = sorted(conv_list)
                c1, c2 = conv_list
                if c1.tag == c2.tag:
                    merge_conv = IntervalBase(chr, c1.start, c2.end)
                    c.append(merge_conv)
                    continue
            c.extend(conv_list)
        c = BucketIndexIntervalList(sorted(c))
        c.make_index()
        return c
    
    def conv_merge_(self, tolomere_convs, convs):
        """合并conv和tolomere_convs"""
        if not tolomere_convs:
            return
        # conv_filter(convs, key=lambda x: not tolomere_convs.find_one(x, action='include'))  # 暂时不过滤，后续get_freq_reliability程序中过滤 by housy 202030706
        tolomere_convs.add_attr('is_aneuploid', values=repeat(True))
        convs._data.extend(tolomere_convs)
        convs._data.sort()
        convs.make_index()
        convs.check_attr('is_aneuploid', attr_default=False)
        
    def forward(self, ):
        """构建非整倍体矫正后的conv"""
        exons = self.exons
        sexuality = self.sexuality
        tolomere_convs = self.tolomere_convs
        
        # for a in tolomere_convs:
        #     print("tolomere_convs1", a, a.tag, a.diff)
        
        conv_filter(tolomere_convs, key=self.is_reliable_tolomere)
        
        # for a in tolomere_convs:
        #     print("tolomere_convs2", a, a.tag, a.diff)
        
        merge_convs = self.merge_tolomere(tolomere_convs)
        set_conv(exons, merge_convs, sexuality)
        self.conv_merge_(merge_convs, self.convs)
    



