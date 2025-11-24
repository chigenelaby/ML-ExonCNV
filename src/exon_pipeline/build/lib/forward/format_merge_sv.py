#-*-coding:utf-8-*-
"""
sv数据 整合wes exon
1. 以sv结果为基础 如果有相交区段则标注相交证据
2. 补充wes exons项目单独的数据
3. 整合家系信息
    - 如果wgs和wes一致，则对应列数，整合为wgs文库号
    - 如果wgs和wes数据不一致，不对应的列数单独成列，0/0=0 wild 相对可靠
"""
import logging
from operator import itemgetter
import fire
from exon_pipeline import __version__ as vs
from exon_pipeline.utils.sample_tools import get_info, read_config
from exon_pipeline.utils.sample_info import get_wkcode, WkCodeBase
from exon_pipeline.apps.merge_gene_sv import merge_gene_sv
from exon_pipeline.utils.io import fmt_write_file, show

logger = logging.getLogger(__name__)

FORMAT_GENE_INFO = 'format_gene_info'

def format_merge_sv(sv_gene_file, wes_gene_file=None, wgs_map_file=None, output=None):
    """
    sv数据 整合wes exon
    1. 以sv结果为基础 如果有相交区段则标注相交证据
    2. 补充wes exons项目单独的数据
    3. 整合家系信息
        - 如果wgs和wes一致，则对应列数，整合为wgs文库号
        - 如果wgs和wes数据不一致，不对应的列数单独成列，0/0=0 wild 相对可靠
    """
    h, r = merge_gene_sv(sv_gene_file, wes_gene_file, wgs_map_file)
    
    fmt_write_file(r, output, headers=h, 
                   infos={'version': vs, })