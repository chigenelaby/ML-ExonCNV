#!/usr/bin/env python3
#-*-coding:utf-8-*-
"""exon_pipeline外显子流程"""
import os
import sys
import logging
import fire

_self_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(_self_path)

from exon_pipeline.core.contrast import build_contrast
from exon_pipeline.forward.make_exons import make_exons
from exon_pipeline.forward.make_flags import make_flags
from exon_pipeline.forward.get_freq_reliability import get_reliability_conv
from exon_pipeline.forward.format_info import format_info
from exon_pipeline.forward.make_qc import make_qc
from exon_pipeline.forward.make_exons_details import make_exons_details
from exon_pipeline.forward.qc_karyotype import qc_karyotype
from exon_pipeline.forward.format_merge_gene import format_merge_gene_file
from exon_pipeline.forward.format_merge_sv import format_merge_sv

if __name__ == '__main__':
    # logging.basicConfig(level=logging.NOTSET)
    # logging.basicConfig(level=logging.DEBUG)
    # logging.basicConfig(level=logging.INFO)
    logging.basicConfig(level=logging.WARNING)
    # logging.basicConfig(level=logging.ERROR)
    # logging.basicConfig(level=logging.CRITICAL)
    fire.Fire({'build_contrast': build_contrast, 
               'make_exons': make_exons, 
               'make_flags': make_flags, 
               'get_reliability_conv': get_reliability_conv,
               'qc': make_qc, 
               'qc_karyotype': qc_karyotype,
               'make_exons_result': make_exons_details,
               'format': format_info,
               'merge_gene_file': format_merge_gene_file, 
               'format_merge_sv': format_merge_sv, 
               })
