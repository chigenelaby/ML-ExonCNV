#-*-coding:utf-8-*-
"""test"""
import logging
import fire
from exon_pipeline.core.contrast import *
from exon_pipeline.core.controlled import *
from exon_pipeline.core.conv import *

from exon_pipeline.utils.io import fmt_write_file, show
from exon_pipeline.utils.load_qc import get_sexuality_from_qc_report

logger = logging.getLogger(__name__)

DEPTH_FILE_STR = 'depth_file'
STAT_FILE_STR = 'stat_file'
NT_TARGET_STR = 'nt_target'

EXONS_INFOS_FILE = 'exons_infos_file'
EXONS_CONV_FILE = 'exons_conv_file'
QC_REPORT_FILE = 'qc_report_file'

# 文件中字段
STAT_DEPTH_STR = 'Average_rmdupdepth_point' 



def run(wkcode, analyse_path, contrast_file, config_file, outdir=None):
    """"""
    # 配置
    logger.info('配置参数...')
    config = read_config(config_file)
    if outdir is None:
        outdir = analyse_path
    
    wk = get_wkcode(wkcode, analyse_path, config=config)
    out_wk = get_wkcode(wkcode, outdir, config=config)
    
    # 对照库
    logger.info('读取对照库...')
    contrast = IO_pickle.load_pickle(contrast_file)
    cb = ControlledBase(contrast)
    
    # 样本
    logger.info('读取样本数据...')
    path = wk.get_filename(DEPTH_FILE_STR)
    idx = get_index(path, default_base=IntervalDepth, key=itemgetter(0,1,2,6))
    sample_depth = get_samlpe_depth(wk)
    idx.wkcode = wk
    wk_str = str(wk)
    idx.add_attr('wkcode', values=repeat(wk_str))
    idx.sample_depth = sample_depth
    
    logger.info('计算输出...')
    exons = cb.get_diff(idx)
    
    qc_report_file = wk.get_filename(QC_REPORT_FILE)
    sexuality = get_sexuality_from_qc_report(qc_report_file)
    exons_conv = get_conv_data(exons, sexuality=sexuality)
    
    # 输出
    exons_infos_file =  out_wk.get_filename(EXONS_INFOS_FILE)
    exons_conv_file = out_wk.get_filename(EXONS_CONV_FILE)
    IO_pickle.write_pickle(exons, exons_infos_file)
    IO_pickle.write_pickle(exons_conv, exons_conv_file)



if __name__ == '__main__':
    # logging.basicConfig(level=logging.NOTSET)
    logging.basicConfig(level=logging.DEBUG)
    # logging.basicConfig(level=logging.INFO)
    # logging.basicConfig(level=logging.WARNING)
    # logging.basicConfig(level=logging.ERROR)
    # logging.basicConfig(level=logging.CRITICAL)
    fire.Fire(run)