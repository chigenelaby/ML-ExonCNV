#-*-coding:utf-8-*-
"""构建exons信息"""
import logging
import fire
from exon_pipeline.core.contrast import *
from exon_pipeline.core.controlled import *
from exon_pipeline.core.conv import *

from exon_pipeline import __version__ as vs
import exon_pipeline.parameters as opt
from exon_pipeline.utils.io import fmt_write_file, show
from exon_pipeline.utils.load_qc import get_sexuality_from_qc_report
from exon_pipeline.utils.fmt_tools import get_fmt_tag

from mod_tools.mod_data_IO import *
from intervaltools.core import IntervalBase, BucketIndexIntervalList 
from exon_pipeline.utils.database_store import get_filter_index 

logger = logging.getLogger(__name__)

DEPTH_FILE_STR = 'depth_file'

EXONS_INFOS_FILE = 'exons_infos_file'
EXONS_CONV_FILE = 'exons_conv_file'
QC_REPORT_FILE = 'qc_report_file'
LOW_COVERAGE_FILE = 'low_coverage_file'


def make_exons(wkcode, analyse_path, contrast_file, config_file,sex, outdir=None):
    """
    --wkcode 文库号
    --analyse-path wkcode-分析目录
    --contrast-file 对照库构建文件，wkcode-分析目录
    --config-file 配置文件
    --outdir 输出分析目录 默认为analyse_path
    """
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
        
    low_coverage_path = wk.get_filename(LOW_COVERAGE_FILE)
    idx = get_filter_index(path, low_coverage_path ,default_base=IntervalDepth, key=itemgetter(0,1,2,6))

    sample_depth = get_samlpe_depth(wk)
    idx.wkcode = wk
    wk_str = str(wk)
    idx.add_attr('wkcode', values=repeat(wk_str))
    idx.sample_depth = sample_depth
    
    logger.info('计算输出...')
    exons = cb.get_diff(idx)

    qc_report_file = wk.get_filename(QC_REPORT_FILE)
    exons_conv = get_conv_data(exons, sexuality=sex)
    exons_infos_file =  out_wk.get_filename(EXONS_INFOS_FILE)
    exons_conv_file = out_wk.get_filename(EXONS_CONV_FILE)
    IO_pickle.write_pickle(exons, exons_infos_file)
    IO_pickle.write_pickle(exons_conv, exons_conv_file)
    
    exons_conv.add_attr('fmt_tag', key=get_fmt_tag)
    attrs = 'chr start end fmt_tag diff'.split()
    r = show(exons_conv, attrs=attrs)
    output2 = os.path.splitext(exons_conv_file)[0] + '.txt'
    fmt_write_file(r, output2, headers=attrs, infos={'version': vs, 'wkcode': exons_conv.wkcode})


if __name__ == '__main__':
    fire.Fire(make_exons)
