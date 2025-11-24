#-*-coding:utf-8-*-
"""exon模块gene格式标准输出"""

import logging
from operator import itemgetter
import fire
from exon_pipeline import __version__ as vs
from exon_pipeline.utils.sample_tools import get_info, read_config
from exon_pipeline.utils.sample_info import get_wkcode, WkCodeBase
from exon_pipeline.core.conv import *
from exon_pipeline.apps.get_report_tag import *
from exon_pipeline.apps.family_check import ConvFamilyCheck
from exon_pipeline.apps.reliability.exon_reliability import *
from exon_pipeline.utils.io import fmt_write_file, show
from exon_pipeline.utils.fmt_tools import *

logger = logging.getLogger(__name__)
# 输入
EXONS_INFOS_FILE = 'exons_infos_file'
FREQ_FILE = 'test_freq_file'
REPEAT_OBJ_FILE = 'repeat_obj_file'
QC_FILE = 'qc_output'
GENE_ATTENTIONS_FILE = 'gene_attentions_file'
GENE_FILTERFALSE_FILE = 'gene_filterfalse_file'

GENE_ATTENTIONS_FILE1 = 'gene_attentions_file1' ## by housy
GENE_ATTENTIONS_FILE2 = 'gene_attentions_file2' ## by housy

# 输出
FORMAT_GENE_INFO = 'format_gene_info'
FORMAT_CNV_INFO = 'format_cnv_info'
FORMAT_CNV_INFO_FILTER = 'format_cnv_info_filter'


family_list_key = itemgetter(0, 1)


def format_info(wkcode, analyse_path, config_file, 
                     family_list_file=None, output=None, 
                     only_gene=False, only_cnv=False, ):
    """
    exon_pipeline模块gene格式标准输出
    --wkcode 文库号
    --analyse-path wkcode-分析目录
    --config-file 配置文件
    --family-list 家系列表 每一行为家系成员 wkcode analyse_path
    --output 输出格式化gene文件 默认为analyse_path FORMAT_GENE_INFO
    --only-gene 仅输出gene整合结果
    --only-cnv 仅输出cnv整合结果
    """
    logger.info('配置参数...')
    config = read_config(config_file)
     
    wk = get_wkcode(wkcode, analyse_path, config=config)
    
    convs_path = wk.get_filename(FREQ_FILE)
    convs = IO_pickle.load_pickle(convs_path)
    
    gene_attentions_file = wk.get_filename(GENE_ATTENTIONS_FILE)
    gene_attentions_file1 = wk.get_filename(GENE_ATTENTIONS_FILE1) ## by housy
    gene_attentions_file2 = wk.get_filename(GENE_ATTENTIONS_FILE2) ## by housy

    gene_filterfalse_file = wk.get_filename(GENE_FILTERFALSE_FILE)
    repeat_obj_file = wk.get_filename(REPEAT_OBJ_FILE)
    qc_file = wk.get_filename(QC_FILE)
    report_tag_key1 = get_report_tag_conv_key(gene_attentions_file, repeat_obj_file, qc_file) ## by housy 20230106
    # report_tag_key1 = get_report_tag_conv_key(gene_attentions_file1, gene_attentions_file2, repeat_obj_file, qc_file) ## by housy
    convs.add_attr('flag_report_tag', key=report_tag_key1)
    convs.add_attr('report_tag', key=fmt_report_tag)
    convs.add_attr('infos', key=lambda x: 
                            fmt_infos_json(x, attrs=('reliability_show', 'is_aneuploid', 'report_tag'))
               )
    
    exons_path = wk.get_filename(EXONS_INFOS_FILE)
    exons = IO_pickle.load_pickle(exons_path)
    
    wks = [wk]
    
    convs_list = []
    exons_list = []
    if family_list_file:
        family_data = load_file(family_list_file, '\t', '#')
        for fm_info in family_data:
            fm_wk = get_wkcode(*family_list_key(fm_info), config=config)
            if fm_wk == wk:
                continue
            wks.append(fm_wk)
            fm_convs_path = fm_wk.get_filename(FREQ_FILE)
            fm_convs = IO_pickle.load_pickle(fm_convs_path)
            convs_list.append(fm_convs)
            
            fm_exons_path = fm_wk.get_filename(EXONS_INFOS_FILE)
            fm_exons = IO_pickle.load_pickle(fm_exons_path)
            exons_list.append(fm_exons)

    # gene filter
    gene_filter_key = get_report_tag_gene_key(gene_attentions_file1, gene_attentions_file2, gene_filterfalse_file, qc_file) # by housy 20230106 白名单不再无条件保留
    # gene_filter_key = get_report_tag_gene_key(gene_attentions_file1, gene_attentions_file2, gene_filterfalse_file, qc_file) ## by housy
    cnv_filter_key = get_report_tag_cnv_key(gene_filterfalse_file, qc_file)
    

    cfc = ConvFamilyCheck(convs, exons, convs_list, exons_list, gene_filter_key, cnv_filter_key)
    res = cfc.forward()
    header = cfc.get_header()
    if output is None:
        output = wk.get_filename(FORMAT_GENE_INFO)
    
    if not only_cnv:
        fmt_write_file(res, output, headers=header, infos={'version': vs, 'wkcode': convs.wkcode})
       
    if not only_gene:
        cnv_dict, res_filter_dict = cfc.forward_cnv()
        cnv_header = cnv_dict['header']
        for wkcode in wks:
            res_cnv = cnv_dict[wkcode]
            cnv_path = wkcode.get_filename(FORMAT_CNV_INFO)
            fmt_write_file(res_cnv, cnv_path, headers=cnv_header, 
                           infos={'version': vs, 'wkcode': wkcode, 'is_filter': False})
            
            res_cnv_filter = res_filter_dict[wkcode]
            cnv_filter_path = wkcode.get_filename(FORMAT_CNV_INFO_FILTER)
            fmt_write_file(res_cnv_filter, cnv_filter_path, headers=cnv_header, 
                           infos={'version': vs, 'wkcode': wkcode, 'is_filter': True})

