#-*-coding:utf-8-*-
"""test标签测试"""

import logging
import fire
import copy
from exon_pipeline import __version__ as vs
from exon_pipeline.core.contrast import *
from exon_pipeline.core.controlled import *
from exon_pipeline.core.conv import *
from exon_pipeline.utils.io import fmt_write_file, show
from exon_pipeline.apps.reliability import *
from exon_pipeline.utils.conv_tools import conv_filter
from exon_pipeline.utils.fmt_tools import get_fmt_tag
from exon_pipeline.apps.freq import *
from exon_pipeline.apps.freq.build_freq_data import *
from exon_pipeline.apps.reliability.reliability_value import *
from exon_pipeline.core.format import *
from exon_pipeline.apps.filter_func import *
from exon_pipeline.utils.fmt_tools import *
import exon_pipeline.parameters as opt

logger = logging.getLogger(__name__)

NT_TARGET_FILE = 'nt_target'
FLAG_INFO_FILE = 'flag_info_file'
QC_REPORT_FILE = 'qc_report_file'
FREQ_FILE = 'test_freq_file'
FREQ_FILE2 = 'test_freq_file2'
GENE_ATTENTIONS_FILE = 'gene_attentions_file'

get_sexuality = opt.get_sexuality

def get_reliability_conv(wkcode, analyse_path, 
        config_file, sex,freq_db_file, outdir=None):
    """"""
    # 配置
    logger.info('配置参数...')
    config = read_config(config_file)
    if outdir is None:
        outdir = analyse_path
    wk = get_wkcode(wkcode, analyse_path, config=config)
    out_wk = get_wkcode(wkcode, outdir, config=config)
    
    mode = 2    
    
    # 样本
    logger.info('读取样本数据...')
    path = wk.get_filename(FLAG_INFO_FILE)
    convs = IO_pickle.load_pickle(path)
    
    
    if convs.sexuality == 'XX':
        conv_filter(convs, lambda x:x.chr!='chrY')
    # print(dir(convs))
    qc_report_file = wk.get_filename(QC_REPORT_FILE)
    sexuality = sex
    
    # freq_db_file = wk.get_filename(FREQ_DB).format(sexuality)     改为直接输入 yzh
    # print('频率库:',freq_db_file)     
    freq_db = FreqBase()
    freq_db.load(freq_db_file)
    
    convs.add_attr('freq', key=lambda x:freq_db.get_freq(*x.interval(), mode).get(x.tag))
    convs.add_attr('all_freq', key=lambda x:sum(freq_db.get_freq(*x.interval(), mode).values()))
    convs.check_attr('freq',attr_default=0)
    convs.check_attr('all_freq',attr_default=0)
    
    # for a in convs:
    #     if a.chr == "chrX":
    #         print(a.chr, a.start, a.end, a.freq, a.all_freq)

    # 过滤
    cs = ConditionSet()
    gene_attentions_file = wk.get_filename(GENE_ATTENTIONS_FILE)
    fca = FilterConditionAttentions(gene_attentions_file)
    cs.add_filter_condition(filter_NA_key, 'function')
    cs.add_attention_condition(fca)
    cs.add_filter_condition(filter_key, 'function')
    
    conv_filter(convs, cs.check)

    if hasattr(convs[0], "is_aneuploid"):
        tolomere_convs = copy.copy(convs)
        conv_filter(tolomere_convs, key=lambda x:x.is_aneuploid==True)
        conv_filter(convs, key=lambda x: (not tolomere_convs.find_one(x, action='include')) or tolomere_convs.find_one(x))

    #reliability
    convs.add_attr('reliability_value', key=get_reliability)
    convs.add_attr('reliability_show', key=get_reliability_show)
    convs.add_attr('infos', key=lambda x: 
                                fmt_infos_json(x, attrs=('reliability_show', 'is_aneuploid'))
                   )
    convs.add_attr('fmt_tag', key=get_fmt_tag)
    
        
    # format
    gene_intervals_path = wk.get_filename(NT_TARGET_FILE)
    gene_intervals = make_gene_intervals(gene_intervals_path)
    gf = GeneFormat(convs, gene_intervals)
    gf.forward()
    
    output = out_wk.get_filename(FREQ_FILE)
    IO_pickle.write_pickle(convs, output)
    
    output2 = out_wk.get_filename(FREQ_FILE2)
    attrs = 'chr start end fmt_tag diff reliability_value freq all_freq infos'.split()
    r = show(convs, attrs=attrs)
    fmt_write_file(r, output2, headers=attrs, infos={'version': vs, 'wkcode': convs.wkcode})

if __name__ == '__main__':
    fire.Fire(get_reliability_conv)