#!/usr/bin/env python3
#-*-coding:utf-8-*-
"""
vcf转table
PyVCF是一个第三方python包,用于读取处理VCF文件
"""
__date__    = '2020/09/02'
__author__  = 'wangxf'
__Version__ = '1.4.0'

import vcf
import os
import sys
import argparse
import textwrap

def parseargs():
    """参数调取"""
    parser = argparse.ArgumentParser(
        prog='vcf2table',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent("""\
            version  : 1.0.0
            Function : vcf2table
            """))
    parser.add_argument('-vcf', dest='vcf', metavar='<input>', help=r"vcf_file", default=None)
    parser.add_argument('-table', dest='table', metavar='<output>', help=r"vcf_table_file", default=None)
    parser.add_argument('-S', '--is_single', help='是否单人', action="store_true")
    argv = parser.parse_args()
    return argv, parser.print_help


def _mut_ratio(ad, dp):
    """计算突变率"""
    if ad == 0 or dp == 0:
        mut_ratio = 0
    else:
        mut_ratio = round(ad/dp, 4)
    return mut_ratio

def _mut_revise(pos, ref, alt):
    """突变结果校正"""
    ref_seq = list(ref)
    alt_seq = list(str(alt))
    while(len(ref_seq) > 0 and len(alt_seq) > 0):
#        if(ref_seq[0] == alt_seq[0]):
#            pos += 1
#            del ref_seq[0]
#            del alt_seq[0]
#        elif(ref_seq[-1] == alt_seq[-1]):
#            del ref_seq[-1]
#            del alt_seq[-1]
#       20200902 by wangxf 修改indel处理策略，先pop再shift，尽量保持突变位置靠近基因组5’端，后续再做暴力穷举
        if(ref_seq[-1] == alt_seq[-1]):
            del ref_seq[-1]
            del alt_seq[-1]
        elif(ref_seq[0] == alt_seq[0]):
            pos += 1
            del ref_seq[0]
            del alt_seq[0]
        else:
            break
    if len(ref_seq) == 0:
        pos -= 1
        end = pos + 1
        alt = ''.join(alt_seq)
        ref = '-'
    elif len(alt_seq) == 0:
        end = pos + len(ref_seq) - 1
        alt = '-'
        ref = ''.join(ref_seq)
    else:
        end = pos + len(ref_seq) - 1
        alt = ''.join(alt_seq)
        ref = ''.join(ref_seq)
    mut_info = str(pos) + '\t' + str(end) + '\t' + ref + '\t' + alt
    return mut_info

if __name__ == "__main__":
    """vcf转table"""
    argv, help_doc = parseargs()
    if not argv.vcf or not argv.table:
        help_doc()
        sys.exit(1)
    with open(argv.table,'wt') as table:
        vcf_reader = vcf.Reader(filename = argv.vcf)
        sams = '\t'.join(vcf_reader.samples)
        # 兼容单人模式 by chenw 2020-4-27
        info_single = "\tinfo" if argv.is_single else ""
        # print("#chr\tstart\tend\tref\talt", sams, sep = '\t', file = table)
        print(f"#chr\tstart\tend\tref\talt{info_single}\t{sams}", file = table)
        for record in vcf_reader:
            alt_count = 1
            for alt in record.ALT:
                mut_revise_info = _mut_revise(record.POS, record.REF, alt)
                mut_info_list = [record.CHROM, mut_revise_info]
                # 兼容单人模式 by chenw 2020-4-27
                if argv.is_single:
                    mut_qc_info_list = []
                    for qc in record.INFO:
                        if isinstance(record.INFO[qc], list):
                            mut_qc_info_list.append(qc+'='+str(record.INFO[qc][alt_count - 1]))
                        else:
                            mut_qc_info_list.append(qc+'='+str(record.INFO[qc]))
                    mut_qc_info = ';'.join(mut_qc_info_list)
                    mut_info_list = [record.CHROM, mut_revise_info, mut_qc_info]

                flag_dp = 1
                for sam in record.samples:
                    """添加判断是否有DP,如无则跳过此行,20200424 修改by wangxf"""
                    try:
                        sam['DP']
                    except:
                        flag_dp = 0
                        continue
                    if str(sam['DP']) == 'None': #兼容某成员dp为None，20200809 修改by wangxf
                        ratio = 0
                        sample_info = '0' + '/' + '0' + '=' + str(ratio)
                    else:
                        if sam['AD']: # when sample's 'AD' value not empty(!=".") 202004 修改 by liyt
                            ratio = _mut_ratio(sam['AD'][alt_count], sam['DP'])
                            sample_info = str(sam['AD'][alt_count]) + '/' + str(sam['DP']) + '=' + str(ratio)
                        else: #when sample's 'AD' value is empty(".")
                            ratio = 0
                            sample_info = '0' + '/' + '0' + '=' + str(ratio)
                    mut_info_list.append(sample_info)
                if flag_dp == 1:
                    alt_count += 1
                    mut_info_out = '\t'.join(mut_info_list)
                    print(mut_info_out, file = table)