#-*-coding:utf-8-*-
"""质控文件输出"""

import logging
from operator import itemgetter
import fire

from exon_pipeline import __version__ as vs
from exon_pipeline.utils.sample_tools import get_info, read_config
from exon_pipeline.utils.sample_info import get_wkcode, WkCodeBase
from exon_pipeline.core.conv import *
from exon_pipeline.utils.io import fmt_write_file, show
from exon_pipeline.apps.unbalance import get_unbalance_item

logger = logging.getLogger(__name__)
# 输入
FREQ_FILE = 'test_freq_file'
UNBALANCE_EXONS_FILE = 'unbalance_exons_file'
# 输出

QC_OUTPUT = 'qc_output'


def make_qc(wkcode, analyse_path, 
            config_file, outdir=None):
    """
    生成单人质控
    """
    # 配置
    logger.info('配置参数...')
    config = read_config(config_file)
    if outdir is None:
        outdir = analyse_path
    wk = get_wkcode(wkcode, analyse_path, config=config)
    out_wk = get_wkcode(wkcode, outdir, config=config)
    
    # 样本
    logger.info('读取样本数据...')
    filter_conv_file = wk.get_filename(FREQ_FILE)

    unbalance_exons_file = wk.get_filename(UNBALANCE_EXONS_FILE)
    qc_item = []
    unbalance_item = get_unbalance_item(filter_conv_file, unbalance_exons_file)
    qc_item.append(unbalance_item)
    output = out_wk.get_filename(QC_OUTPUT)
    fmt_write_file(qc_item, output, infos={'version': vs})
    
