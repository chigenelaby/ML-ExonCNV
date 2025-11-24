#-*-coding:utf-8-*-
"""test标签测试"""

import logging
import fire
from exon_pipeline.core.contrast import *
from exon_pipeline.core.controlled import *
from exon_pipeline.core.conv import *

from exon_pipeline.utils.load_qc import get_sexuality_from_qc_report
from exon_pipeline.utils.io import fmt_write_file, show
from exon_pipeline.apps.reliability import *
from exon_pipeline.apps.freq import *
from exon_pipeline.apps.freq.build_freq_data import *

logger = logging.getLogger(__name__)

FLAG_INFO_FILE = 'flag_info_file'
QC_REPORT_FILE = 'qc_report_file'
FREQ_DB = 'freq_db'
FREQ_FILE = 'test_freq_file'
FREQ_FILE2 = 'test_freq_file2'

def run(wkcode, analyse_path, 
        config_file, outdir=None):
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
    path = wk.get_filename(FLAG_INFO_FILE)
    convs = IO_pickle.load_pickle(path)
    
    filter_key = lambda x: (x.tag!='NA' and x.flag_diff.classify != False)
    filter_key2 = lambda i:(i.tag!='NA' and i.flag_diff.classify != False 
                    and (i.score>0 or (len(i.exons)>2 and 
                    len([ex for ex in i.exons if ex.contrast.cv<0.2])/len(i.exons)>0.5) or i.flag_is_cnv.classify))
    func_filter(convs, filter_key2)
    
    qc_report_file = wk.get_filename(QC_REPORT_FILE)
    sexuality = get_sexuality_from_qc_report(qc_report_file)
    
    freq_db_file = wk.get_filename(FREQ_DB).format(sexuality)
    freq_db = FreqBase()
    freq_db.load(freq_db_file)
    
    convs.add_attr('freq', key=lambda x:freq_db.get_freq(*x.interval()).get(x.tag))
    convs.add_attr('all_freq', key=lambda x:sum(freq_db.get_freq(*x.interval()).values()))
    
    output = out_wk.get_filename(FREQ_FILE)
    h = 'chr start end diff tag flag_diff.classify flag_is_cnv.classify flag_length.classify flag_reliable_exon.classify flag_cnv_family.classify flag_vaf.classify flag_cnvseq.classify score all_freq freq'.split()
    res = show(convs, h)
    fmt_write_file(res, output, headers=h)
    
    output2 = out_wk.get_filename(FREQ_FILE2)
    IO_pickle.write_pickle(convs, output2)



if __name__ == '__main__':
    # logging.basicConfig(level=logging.NOTSET)
    logging.basicConfig(level=logging.DEBUG)
    # logging.basicConfig(level=logging.INFO)
    # logging.basicConfig(level=logging.WARNING)
    # logging.basicConfig(level=logging.ERROR)
    # logging.basicConfig(level=logging.CRITICAL)
    fire.Fire(run)