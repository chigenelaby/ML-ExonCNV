#-*-coding:utf-8-*-

"""
家系类校正
"""
import copy
import random
import json
import numpy as np

from exon_pipeline.utils.conv_tools import conv_filter
from exon_pipeline.utils.io import show
from exon_pipeline.utils.sample_tools import fmt_id_string
from exon_pipeline.utils.fmt_tools import get_fmt_tag
from exon_pipeline.apps.family_info import *
from exon_pipeline.apps.reliability.exon_reliability import *
from exon_pipeline.utils.fmt_tools import *
from exon_pipeline.apps.get_report_tag import *

# cf = ConvFamilyBase()

class ConvFamilyCheck(ConvFamilyBase):
    """
    self.condition.forward_one(tag, ex) 返回值 0 1 2
    """
    def __init__(self, propositus_convs, propositus_exons, 
                 convs_list, exons_list, gene_filter_key=None,
                 cnv_filter_key=None):
        self.condition = ExonReliabilityCondition()
        self.convs_list = convs_list
        self.exons_list = exons_list
        self.propositus_convs = propositus_convs
        self.propositus_exons = propositus_exons
        self.gene_filter_key = gene_filter_key
        self.cnv_filter_key = cnv_filter_key
    
    def intersect_exons(self, gene, conv):
        """相交的外显子"""
        res = []
        for ex in conv.exons:
            if ex in gene.exons:
                res.append(ex)
        return res
    
    def get_condition_exons(self, exons_interval, tag):
        """统计exons可靠性>1的exon个数"""
        res = []
        for ex in exons_interval:
            if self.condition.forward_one(tag, ex) > 0 and not ex.contrast.flag_bad_contrast:
                res.append(ex)
        return res
    
    def find_best_exon(self, reliability_exons, tag):
        """最佳外显子"""
        if tag not in self.condition.diff_condition.tags:
            return None
        if not len(reliability_exons):
            return None
        base_diff = getattr(self.condition.diff_condition, '%s_base'%tag)
        exons_sort = sorted(reliability_exons, 
                key=lambda x:(self.condition.forward_one(tag, x), 
                              -abs(x.diff-base_diff)),
                reverse=True)
        best_exon = exons_sort[0]
        return best_exon
    
    def fmt_best_exon(self, best_exon, gene):
        """格式化最佳外显子, EXON:4(chr1:104294554-104294752)"""
        if not best_exon:
            return 'NA'
        idx = gene.exons.index(best_exon)
        ex = gene.exons[idx]
        # string = 'EXON:%s(%s:%s-%s)'%(ex.exon_num, ex.chr, ex.start, ex.end)
        string = 'REG:%s(%s:%s-%s)'%(ex.exon_num, ex.chr, ex.start, ex.end)
        return string
    
    
    def get_propositus_gene_items(self):
        """获取先证者的gene_items"""
        conv1 = self.propositus_convs
        exons1 = self.propositus_exons
        for cv in conv1:
            for gene in cv.gene_divides:
                exs = exons1.find_all(gene, action='be_included')
                if hasattr(cv, 'reliability_value'):
                    reli_value = cv.reliability_value
                    reli_show = cv.reliability_show
                    # print(reli_show)
                else:
                    reli_value = "未知"
                    reli_show = []
                setattr(gene, 'reliability_value', reli_value)
                setattr(gene, 'reliability_show', reli_show) ##by wangxf 增加reliability_show属性
                setattr(gene, 'tag', cv.tag) ##by wangxf 增加tag属性 20230414
                yield (gene, cv, exs)
    
    def get_gene_item(self, gene, tag, exons, convs):
        """输入先证者conv以及对应wk的外显子信息和convs"""
        exs = exons.find_all(gene, action='be_included')
        return (gene, self.get_intersect_conv(gene, tag, convs), exs)
    
    def fmt_gene_item(self, item):
        """格式化item (gene, conv, wkcode)"""
        gene, conv, exs = item
        # chr = gene.chr
        # start = gene.start
        # end = gene.end
        # gene_name = gene.gene
        # exon_info = gene.exon_info
        gene_len = len(gene.exons)
        diff = round(np.mean([i.diff for i in exs]), 4)  ## 20220727 by housy 保留4位
        gene_all_str = ['', '(all)'][gene.is_all_gene]
        # gene_info_str = '(%s)%s'%(exon_info, gene_all_str)
        
        if conv:
            tag = conv.tag
            fmt_tag = get_fmt_tag(conv)
            exons_interval = self.intersect_exons(gene, conv)
            reliability_exons = self.get_condition_exons(exons_interval, tag)
            reliability_count = len(reliability_exons)
            # best_exon = self.find_best_exon(reliability_exons, tag)
            # freq = conv.freq
            # all_freq = conv.all_freq
            if hasattr(conv, 'reliability_value'):
                reliability_value = conv.reliability_value
            else:
                reliability_value = '未知'
        else:
            fmt_tag = tag = 'wild'
            reliability_count = 0
            # best_exon = None
            # freq = 'NA'
            # all_freq = 'NA'
            reliability_value = '相对可靠'
        count_str = '%s/%s=%s'%(reliability_count, gene_len, diff)
        # best_exon_str = self.fmt_best_exon(best_exon, gene)
        
        return (count_str, fmt_tag, 
                reliability_value, )
    
    def fmt_base(self, p_item):
        """格式化基础信息"""
        gene, conv, exs = p_item
        chr = gene.chr
        start = gene.start
        end = gene.end
        gene_name = gene.gene
        exon_info = gene.exon_info
        gene_all_str = ['', '(all)'][gene.is_all_gene]
        gene_info_str = '(%s)%s'%(exon_info, gene_all_str)
        gene_info_str = gene_info_str.replace("EXON", "REG")

        tag = conv.tag
        exons_interval = self.intersect_exons(gene, conv)
        reliability_exons = self.get_condition_exons(exons_interval, tag)
        best_exon = self.find_best_exon(reliability_exons, tag)
        freq = conv.freq
        all_freq = conv.all_freq
        best_exon_str = self.fmt_best_exon(best_exon, gene)
        
        fmt_info = conv.infos
        # gene filter
        if self.gene_filter_key is None:
            pass

        elif not hasattr(conv, 'flag_report_tag') or not conv.flag_report_tag:
            pass
        else:
            #if gene == "AMY2A":
            #    print(self.gene_filter_key(gene))
            flag_report_tag = self.gene_filter_key(gene)
            report_tag = fmt_report_tag_base(flag_report_tag)
            fmt_info = update_infos_json_dict(fmt_info, report_tag=report_tag)
        
        return (chr, start, end, 
                gene_name, gene_info_str, 
                best_exon_str, freq, all_freq, fmt_info)
    
    def get_wkcodes(self):
        """获取文库号列表"""
        res = []
        res.append(self.propositus_convs.wkcode)
        for convs in self.convs_list:
            res.append(convs.wkcode)
        return res
    
    def get_header(self):
        """获取header"""
        h = "chr start end gene_name gene_info_str best_exon_str freq all_freq infos".split()
        wkcode_head_fmt = '{0}_exon {0}_tag {0}_result'
        wks = self.get_wkcodes()
        for wk in wks:
            h.extend(wkcode_head_fmt.format(wk).split())
        return h
    
    
    def forward(self):
        """输入"""
        res = []
        for p_item in self.get_propositus_gene_items():
            line_list = []
            res.append(line_list)
            gene, cv, exs = p_item
            line_list.extend(self.fmt_base(p_item))
            line_list.extend(self.fmt_gene_item(p_item))
            
            for convs, exons in zip(self.convs_list, self.exons_list):
                item = self.get_gene_item(gene, cv.tag, exons, convs)
                line_list.extend(self.fmt_gene_item(item))
        
        return res
    
    ## cnv文件格式化
    def is_exact_cnv(self, conv):
        """判断是否为准确的cnv
        conv type: conv
        rtype: bool 
        """
        tolerance_num = 1
        normal_min = 10
        tag = conv.tag
        if tag == 'loss2':
            tolerance_min_exon = 5
        else:
            tolerance_min_exon = 7
        length = len(conv.exons)
        
        condition_exons = self.get_condition_exons(conv.exons, tag)
        condition_exons_count = len(condition_exons)
        
        badcapture_exons = [ex for ex in conv.exons if ex.contrast.flag_bad_capture]
        badcapture_exons_count = len(badcapture_exons)
        
        # normal_min或以上直接比较， [tolerance_min_exon, normal_min]容忍个数为tolerance_num
        if condition_exons_count >= normal_min:
            flag = True
        elif condition_exons_count >= tolerance_min_exon:
            # 错误长度 = 总长 - 同类tag - badcapture <= 容忍长度
            flag = (length - condition_exons_count - badcapture_exons_count <= tolerance_num)
        else:
            flag = False
        return flag
    
    def get_conv_fmt(self, convs, exons, filter_p_convs):
        """
        输入一个wk对应的convs和exons信息，以及校正后的p_convs
        返回校正后的cnv
        'id_string chr start end diff fmt_tag reliability_value infos freq all_freq'
        """
        res = []
        for conv in filter_p_convs:
            chr = conv.chr
            start = conv.start
            end = conv.end
            diff = round(np.mean([i.diff for i in exons.find_all(conv, action='be_included')]), 4) ## 20220727 by housy 保留4位
            intersect_conv = self.get_intersect_conv(conv, conv.tag, convs)
            if intersect_conv is None:
                tag = 'wild'
                reliability_value = '相对可靠'
                fmt_tag = tag
            else:
                tag = intersect_conv.tag
                reliability_value = intersect_conv.reliability_value
                fmt_tag = get_fmt_tag2(tag, intersect_conv.is_male_sex_chr)
            freq = conv.freq
            all_freq = conv.all_freq
            id_string = fmt_id_string(chr, start, end, fmt_tag, fmt_tag)
            infos = fmt_infos_json(conv)
            res.append((id_string, chr, start, end, diff, fmt_tag, reliability_value, freq, all_freq, infos))
        return res
    
    def forward_cnv(self, ):
        """
        输出为 cnv数据 dict wkcode: list[list]
        'chr start end diff tag reliability_value freq all_freq infos'
        """
        attrs = 'id_string chr start end diff fmt_tag reliability_value freq all_freq infos'.split()
        p_convs = copy.copy(self.propositus_convs)
        conv_filter(p_convs, self.is_exact_cnv)
        p_convs.add_attr('id_string', key=lambda x: fmt_id_string(x.chr, x.start, x.end, x.fmt_tag, x.fmt_tag))
        res = {}
        # 先证者
        res[p_convs.wkcode] = show(p_convs, attrs=attrs)
        # 家系成员
        for fconvs, fexons in zip(self.convs_list, self.exons_list):
            wkcode = fconvs.wkcode
            fmt_res = self.get_conv_fmt(fconvs, fexons, p_convs)
            res[wkcode] = fmt_res
        res['header'] = attrs
        
        ## cnv_filter:
        cnv_filter_key = self.cnv_filter_key
        if cnv_filter_key is None:
            cnv_filter_key = lambda x: x
        
#        for a in p_convs:
#            print("====")
#            print(a)
#            for attr in dir(a):
#                if not attr.startswith("__"):
#                    print(attr, getattr(a, attr))
#            print(a.flag_report_tag)

        conv_filter(p_convs, lambda x: (hasattr(x, "flag_report_tag") and x.flag_report_tag) and self.cnv_filter_key(x))
        p_convs.add_attr('id_string', key=lambda x: fmt_id_string(x.chr, x.start, x.end, x.fmt_tag, x.fmt_tag))
        res_filter = {}
        # 先证者
        res_filter[p_convs.wkcode] = show(p_convs, attrs=attrs)
        # 家系成员
        for fconvs, fexons in zip(self.convs_list, self.exons_list):
            wkcode = fconvs.wkcode
            fmt_res = self.get_conv_fmt(fconvs, fexons, p_convs)
            res_filter[wkcode] = fmt_res
        res_filter['header'] = attrs
        return res, res_filter


