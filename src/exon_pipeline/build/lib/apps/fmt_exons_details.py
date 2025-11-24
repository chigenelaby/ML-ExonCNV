#-*-coding:utf-8-*-
"""生成单个外显子信息"""

from exon_pipeline.core.format import *
from exon_pipeline.core.conv import *
from exon_pipeline.apps.reliability import *
from exon_pipeline.utils.database_store import get_filter_index  ## by housy


class ExonDetails:
    def __init__(self, exons, convs, gene_intervals_path, low_coverage_path):
        self.exons = exons
        self.convs = convs
        # self.exon_cds_list = get_index(gene_intervals_path, default_base=ExonCDS, key=itemgetter(0,1,2,3,4))
        self.exon_cds_list = get_filter_index(gene_intervals_path, low_coverage_path, default_base=ExonCDS, key=itemgetter(0,1,2,3,4))  ## by housy
        self.headers = self.get_headers()
    
    def forward(self,):
        """
        - chr	染色体号
        - start	起始位置
        - end	终止位置
        - cds	cds区
        - gene	基因
        - diff	diff值
        - depth	样本深度
        - sample_value	样本值
        - flag_outlier	是否离群标签
        - z_value	z值
        - contrast.mean	对照库均值
        - contrast.cv	对照库cv
        - contrast.std	对照库mean
        - contrast.dps_mean	对照库深度
        - contrast.flag_bad_capture 对照库bad_capture标签
        - tag	tag值
        """
        exons = self.exons
        convs = self.convs
        exon_cds_list = self.exon_cds_list
        res = []
        exon_info_attrs = 'chr start end cds gene'.split()
        exon_values_attrs = 'diff depth sample_value flag_outlier'.split()
        contrast_values_attrs = 'contrast.mean contrast.cv contrast.std contrast.dps_mean contrast.flag_bad_capture'.split()
        for ex in exon_cds_list:
            exon_info = attrgetter(*exon_info_attrs)(ex)
            exon = exons.find_one(ex)
            if not exon:        #忽略前期过滤掉的区间 yzh
                continue
            exon_values = attrgetter(*exon_values_attrs)(exon)
            z_score = self.get_z_score(exon)
            contrast_values = attrgetter(*contrast_values_attrs)(exon)
            conv = convs.find_one(exon, action='include')
            tag = conv.tag if conv else 'NA'
            res.append((*exon_info, *exon_values, z_score, *contrast_values, tag))
        return res
    
    def get_headers(self):
        return ['chr', 'start', 'end', 'cds', 'gene', 
                'diff', 'depth', 'sample_value', 'flag_outlier', 'z_score', 
                'contrast.mean', 'contrast.cv', 'contrast.std', 
                'contrast.dps_mean', 'contrast.flag_bad_capture', 'tag']
        
    def get_z_score(self, exon):
        """z值"""
        if not exon.contrast.std:
            return 0.0
        return (exon.sample_value - exon.contrast.mean) / exon.contrast.std
        
        