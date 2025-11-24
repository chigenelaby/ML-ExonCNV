#-*-coding:utf-8-*-
"""质控文件输出"""

import logging
from operator import itemgetter
import fire

from exon_pipeline import __version__ as vs
from exon_pipeline.utils.sample_tools import get_info, read_config
from exon_pipeline.utils.sample_info import get_wkcode, WkCodeBase
from exon_pipeline.core.conv import *
from exon_pipeline.utils.io import fmt_write_file, show, round_table
from exon_pipeline.utils.database_store import get_filter_index 
from exon_pipeline.core.format import *
from exon_pipeline.core.conv import *
from exon_pipeline.apps.reliability import *

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
        

logger = logging.getLogger(__name__)

# 输入
NT_TARGET_FILE = 'nt_target'
EXONS_INFOS_FILE = 'exons_infos_file'
FREQ_FILE = 'test_freq_file'
# 输出

EXONS_DETAILS_FILE = 'exons_details_file'

LOW_COVERAGE_FILE = 'low_coverage_file'

def make_exons_details(wkcode, analyse_path, 
                       config_file, outdir=None):
    """
    生成单人质控
    """
    logger.info('配置参数...')
    config = read_config(config_file)
    if outdir is None:
        outdir = analyse_path
    wk = get_wkcode(wkcode, analyse_path, config=config)
    out_wk = get_wkcode(wkcode, outdir, config=config)
    
    logger.info('读取样本数据...')
    path = wk.get_filename(FREQ_FILE)
    convs = IO_pickle.load_pickle(path)
    
    exons_path = wk.get_filename(EXONS_INFOS_FILE)
    exons = IO_pickle.load_pickle(exons_path)
    gene_intervals_path = wk.get_filename(NT_TARGET_FILE)
    
    low_coverage_path = wk.get_filename(LOW_COVERAGE_FILE)
    ed = ExonDetails(exons, convs, gene_intervals_path, low_coverage_path)
    res = ed.forward()
    headers = ed.headers

    output = out_wk.get_filename(EXONS_DETAILS_FILE)
    fmt_write_file(round_table(res, 4), output, headers=headers, infos={'version': vs})
    
if __name__ == '__main__':
    # logging.basicConfig(level=logging.DEBUG)
    fire.Fire(make_exons_details)
