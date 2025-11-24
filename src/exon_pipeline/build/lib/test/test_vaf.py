#-*-coding:utf-8-*-
"""test2 cnvseq标签测试"""
import logging
import fire
from exon_pipeline.core.contrast import *
from exon_pipeline.core.controlled import *
from exon_pipeline.core.conv import *

from exon_pipeline.utils.io import fmt_write_file, show
from exon_pipeline.apps.vaf import *

logger = logging.getLogger(__name__)

EXONS_CONV_FILE = 'exons_conv_file'
VCF_PATH = 'vcf_path'
VAF_FLAG_FILE = 'vaf_flag_file'

VAF_DB = 'vaf_db'
REPEAT_OBJ_FILE = 'repeat_obj_file'

def run(wkcode, analyse_path, config_file, outdir=None):
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
    
    vcf_path = wk.get_filename(VCF_PATH)
    vcf = fmt_load_vcf_obj(vcf_path, bed_path=None, output_fmt_vcf=None)
    
    logger.info('读取对照库数据...')
    db_save_path = wk.get_filename(VAF_DB)
    
    ds = DatabaseStore()
    ds.load(db_save_path)
    
    repeat_obj_file = wk.get_filename(REPEAT_OBJ_FILE)
    repeat_obj = IO_pickle.load_pickle(repeat_obj_file)
    
    vaf_obj = FlagExonVafBase(repeat_obj=repeat_obj, vcf_obj=vcf, vcf_db=ds)
    
    convs.add_attr('vaf_flag', key=vaf_obj.get_flag)
    
    # 过滤
    # cc = [i for i in convs if i.cnv_seq_flag.classify!='1.9' and i.cnv_seq_flag.classify!='1.3.0' and
                              # i.cnv_seq_flag.classify!='1.3.1']
    
    h = 'chr start end diff tag vaf_flag.classify vaf_flag.show_string vaf_flag.comments'.split()
    res = show(convs, h)
    
    output = out_wk.get_filename(VAF_FLAG_FILE)
    fmt_write_file(res, output, headers=h)



if __name__ == '__main__':
    # logging.basicConfig(level=logging.NOTSET)
    logging.basicConfig(level=logging.DEBUG)
    # logging.basicConfig(level=logging.INFO)
    # logging.basicConfig(level=logging.WARNING)
    # logging.basicConfig(level=logging.ERROR)
    # logging.basicConfig(level=logging.CRITICAL)
    fire.Fire(run)