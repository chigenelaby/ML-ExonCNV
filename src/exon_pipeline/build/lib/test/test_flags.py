#-*-coding:utf-8-*-
"""test标签测试"""

import logging
import fire
from exon_pipeline.core.contrast import *
from exon_pipeline.core.controlled import *
from exon_pipeline.core.conv import *

from exon_pipeline.utils.io import fmt_write_file, show
from exon_pipeline.apps.reliability import *

logger = logging.getLogger(__name__)

EXONS_CONV_FILE = 'exons_conv_file'
VCF_PATH = 'vcf_path'
FLAG_INFO_FILE = 'flag_info_file'
FLAG_INFO_FILE2 = 'flag_info_file2'
FLAG_INFO_FILE3 = 'flag_info_file3'

VAF_DB = 'vaf_db'
REPEAT_OBJ_FILE = 'repeat_obj_file'


def run(wkcode, analyse_path, 
        config_file, wgs_res_path=None, wgs_medi_path=None, family_path1=None, family_path2=None, outdir=None):
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
    db_save_path = wk.get_filename(VAF_DB)
    repeat_obj_file = wk.get_filename(REPEAT_OBJ_FILE)
    
    # 条件
    c1 = GetFlagDiffCondition()
    c2 = GetFlagIsCnv()
    c3 = GetFlagLengthCondition()
    c4 = GetReliableExonCondition()
    c5 = GetCNVFamilyCondition(family_path1, family_path2)
    c6 = GetVafCondition(vcf_path, db_save_path, repeat_obj_file)
    c7 = GetCnvSeqCondition(wgs_res_path, wgs_medi_path)
    
    sf = StatisticsFlagBase(convs)
    sf.add_flag_condition(c1)
    sf.add_flag_condition(c2)
    sf.add_flag_condition(c3)
    sf.add_flag_condition(c4)
    sf.add_flag_condition(c5)
    sf.add_flag_condition(c6)
    sf.add_flag_condition(c7)
    sf.forward()
    
    output = out_wk.get_filename(FLAG_INFO_FILE)
    IO_pickle.write_pickle(convs, output)
    
    output2 = out_wk.get_filename(FLAG_INFO_FILE2)
    IO_pickle.write_pickle(convs, output)
    
    h = 'chr start end diff tag flag_diff.classify flag_is_cnv.classify flag_length.classify flag_reliable_exon.classify flag_cnv_family.classify flag_vaf.classify flag_cnvseq.classify score'.split()
    res = show(convs, h)
    
    output1 = out_wk.get_filename(FLAG_INFO_FILE2)
    fmt_write_file(res, output1, headers=h)
    
    cc = [i for i in convs if i.tag!='NA' and i.flag_diff.classify != False and i.score>0]
    res = show(cc, h)
    
    output2 = out_wk.get_filename(FLAG_INFO_FILE3)
    fmt_write_file(res, output2, headers=h)



if __name__ == '__main__':
    # logging.basicConfig(level=logging.NOTSET)
    logging.basicConfig(level=logging.DEBUG)
    # logging.basicConfig(level=logging.INFO)
    # logging.basicConfig(level=logging.WARNING)
    # logging.basicConfig(level=logging.ERROR)
    # logging.basicConfig(level=logging.CRITICAL)
    fire.Fire(run)