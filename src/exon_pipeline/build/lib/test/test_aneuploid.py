#-*-coding:utf-8-*-
"""非整倍体"""
from exon_pipeline.core.conv import *
from exon_pipeline.utils.load_qc import get_sexuality_from_qc_report
from exon_pipeline.utils.io import show





def get_aneuploid_cnvs(aneuploid_file, exons, sexuality):
    """获得非整倍体cnvs"""
    data = load_file(aneuploid_file, '\t', '#')
    cnvs = BucketIndexIntervalList(data, default_base=IntervalBase, key=itemgetter(0,1,2))
    cnvs.make_index()
    set_conv(exons, cnvs, sexuality)
    return cnvs

def run(aneuploid_file, exons_info_file, sexuality_file, output):
    exons = IO_pickle.load_pickle(exons_info_file)
    sexuality = get_sexuality_from_qc_report(sexuality_file)
    cnvs = get_aneuploid_cnvs(aneuploid_file, exons, sexuality)
    attrs = 'chr start end diff tag'.split()
    r = show(cnvs, attrs)
    fmt_write_file(r, output, headers=attrs, )


def merge_chr_cnvs(cnvs):
    """按chr合并"""
    cnvs_chr = {}
    for cnv in cnvs:
        cnvs_chr.setdefault(cnv.chr, []).append(cnv)
    
    res = []
    