#-*-coding:utf-8-*-
"""
NQB数据
只有父母没有先证者
合并父母相互矫正的gene文件
"""
import logging
from operator import itemgetter
import fire
from exon_pipeline import __version__ as vs
from exon_pipeline.utils.sample_tools import get_info, read_config
from exon_pipeline.utils.sample_info import get_wkcode, WkCodeBase
from exon_pipeline.apps.merge_gene_file import merge_gene_file
from exon_pipeline.utils.io import fmt_write_file, show

logger = logging.getLogger(__name__)

FORMAT_GENE_INFO = 'format_gene_info_new'  ## 兼容外显子缺失重复可靠性校正 by housy 20230306

def format_merge_gene_file(wkcode1, analyse_path1, 
                           wkcode2, analyse_path2, 
                           output, config_file, no_filter=False):
    """
    携带者筛查format gene文件合并
    将夫妻携带者双方相互矫正的gene文件进行merge
    
    --wkcode1 文库号1
    --analyse-path1 wkcode-分析目录1
    --wkcode2 文库号2
    --analyse-path wkcode-分析目录2
    --output 合并后的gene文件路径
    --config-file 配置文件
    """
    config = read_config(config_file)
    wk1 = get_wkcode(wkcode1, analyse_path1, config=config)
    wk2 = get_wkcode(wkcode2, analyse_path2, config=config)
    
    gene_file1 = wk1.get_filename(FORMAT_GENE_INFO)
    gene_file2 = wk2.get_filename(FORMAT_GENE_INFO)
    
    h, r = merge_gene_file(gene_file1, gene_file2, no_filter=no_filter)
    
    fmt_write_file(r, output, headers=h, 
                   infos={'version': vs, 'wkcode1': wkcode1, 'wkcode2': wkcode2})
