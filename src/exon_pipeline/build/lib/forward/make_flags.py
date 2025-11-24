#-*-coding:utf-8-*-
"""
标签
非整倍体
"""

import logging
import fire
import pathlib
from exon_pipeline import __version__ as vs
from exon_pipeline.core.contrast import *
from exon_pipeline.core.controlled import *
from exon_pipeline.core.conv import *

from exon_pipeline.utils.io import fmt_write_file, show
from exon_pipeline.apps.aneuploid import TolomereConv
from exon_pipeline.apps.reliability import *

logger = logging.getLogger(__name__)

EXONS_INFOS_FILE = 'exons_infos_file'
EXONS_CONV_FILE = 'exons_conv_file'
VCF_PATH = 'vcf_path'
FLAG_INFO_FILE = 'flag_info_file'
TOLOMERE_FILE = 'tolomere_file'

FMT_EXONS_CONV_FILE = 'fmt_exons_conv_file'
VAF_DB = 'vaf_db'
REPEAT_OBJ_FILE = 'repeat_obj_file'
ANEUPLOID_FILE = 'aneuploid_file'


def make_flags(wkcode, analyse_path, 
        config_file, wgs_res_path=None, wgs_medi_path=None, family_path1=None, family_path2=None, outdir=None):
    """
    --wkcode 文库号
    --analyse-path wkcode-分析目录
    --config-file 配置文件
    --outdir 输出分析目录 默认为analyse_path
    --wgs-res-path CNVseq结果文件 compatible.cnv.tt
    --wgs-medi-path CNVseq中间文件 lib.cnv
    --family-path1 家系父母conv文件1
    --family-path2 家系父母conv文件2
    """
    # 配置
    if wgs_res_path != None and (not pathlib.Path(wgs_res_path).exists() or not pathlib.Path(wgs_medi_path).exists()):
        assert False, f"CNV-seq文件 {wgs_res_path} {wgs_medi_path} 不存在，请检查"
    
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
    exons_infos_file = wk.get_filename(EXONS_INFOS_FILE)
    exons = IO_pickle.load_pickle(exons_infos_file)
    
    # 非整倍体校正
    tolomere_file = wk.get_filename(TOLOMERE_FILE)
    aneuploid_file = out_wk.get_filename(ANEUPLOID_FILE)
    fmt_exons_conv_file = out_wk.get_filename(FMT_EXONS_CONV_FILE)
    tc = TolomereConv(tolomere_file, exons, convs)
    tc.forward()

    IO_pickle.write_pickle(convs, fmt_exons_conv_file)
    aneuploid_head = 'id_string chr start end diff fmt_tag show_string'.split()
    aneuploid_data = show(sorted(tc.aneuploid), attrs=aneuploid_head)
    fmt_write_file(aneuploid_data, aneuploid_file, headers=aneuploid_head, infos={'version': vs, 'wkcode': convs.wkcode})
    
    # flags
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


if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    fire.Fire(make_flags)