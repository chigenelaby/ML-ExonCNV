#-*-coding:utf-8-*-
"""test2 cnvseq标签测试"""
import logging
import fire
from exon_pipeline.core.contrast import *
from exon_pipeline.core.controlled import *
from exon_pipeline.core.conv import *

from exon_pipeline.utils.io import fmt_write_file, show
from exon_pipeline.apps.reliability.flag_cnv_seq import *

logger = logging.getLogger(__name__)

EXONS_CONV_FILE = 'exons_conv_file'
CNV_SEQ_FLAG_FILE = 'cnv_seq_flag_file'



def run(wkcode, analyse_path, config_file, wgs_res_path, wgs_medi_path, outdir=None):
    """"""
    # 配置
    logger.info('配置参数...')
    config = read_config(config_file)
    if outdir is None:
        outdir = analyse_path
    wk = get_wkcode(wkcode, analyse_path, config=config)
    out_wk = get_wkcode(wkcode, outdir, config=config)
    
    # 样本
    logger.info('读取样本数据...')
    path = wk.get_filename(EXONS_CONV_FILE)
    convs = IO_pickle.load_pickle(path)
    
    fs = GetFlagCnvSeq(wgs_res_path, wgs_medi_path)
    convs.add_attr('cnv_seq_flag', key=fs.get_flag)
    
    # 过滤
    cc = [i for i in convs if i.cnv_seq_flag.classify!='1.9' and i.cnv_seq_flag.classify!='1.3.0' and
                              i.cnv_seq_flag.classify!='1.3.1']
    
    h = 'chr start end diff tag cnv_seq_flag.classify cnv_seq_flag.show_string'.split()
    res = show(cc, h)
    
    output = out_wk.get_filename(CNV_SEQ_FLAG_FILE)
    fmt_write_file(res, output, headers=h)



if __name__ == '__main__':
    # logging.basicConfig(level=logging.NOTSET)
    logging.basicConfig(level=logging.DEBUG)
    # logging.basicConfig(level=logging.INFO)
    # logging.basicConfig(level=logging.WARNING)
    # logging.basicConfig(level=logging.ERROR)
    # logging.basicConfig(level=logging.CRITICAL)
    fire.Fire(run)