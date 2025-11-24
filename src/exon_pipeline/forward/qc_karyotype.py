#-*-coding:utf-8-*-
"""核型标签"""
import logging
import fire
from exon_pipeline.core.contrast import *
from exon_pipeline.core.controlled import *
from exon_pipeline.apps.qc.get_karyotype import get_karyotype

logger = logging.getLogger(__name__)

#输入
ANEUPLOID_FILE = 'aneuploid_file'
#输出
KARYOTYPE_FILE = 'karyotype_file'

def qc_karyotype(wkcode, analyse_path, 
                config_file, outdir=None):
    """质控，核型"""
    # 配置
    logger.info('配置参数...')
    config = read_config(config_file)
    if outdir is None:
        outdir = analyse_path
    wk = get_wkcode(wkcode, analyse_path, config=config)
    out_wk = get_wkcode(wkcode, outdir, config=config)
    
    path = wk.get_filename(ANEUPLOID_FILE)
    output = out_wk.get_filename(KARYOTYPE_FILE)
    get_karyotype(path, output)
    return 