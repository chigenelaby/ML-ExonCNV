#-*-coding:utf-8-*-
"""报出标签"""
import numpy as np
from operator import attrgetter
from exon_pipeline.core.conv import *
from exon_pipeline.apps.filter_func import *


REPORT_TAG_STRING = "report_tag"

def fmt_report_tag(conv, attr='flag_report_tag'):
    flag = getattr(conv, attr)
    if flag:
        return "Report_loci"
    else:
        return "NA"

def fmt_report_tag_base(flag):
    if flag:
        return "Report_loci"
    else:
        return "NA"

# class FilterConditionOther(FilterConditionBase):  ## 20211009 by housy
#     """loss1/gain1过滤条件：外显子个数<=2 and 不可靠/可能不可靠/未知 """
#     def check(self, conv):
#         length = len(conv.exons)
#         if conv.tag == "loss1" or conv.tag == "gain1":
#             if length <= 2 and (conv.reliability_value == "可能不可靠" or conv.reliability_value == "不可靠" or conv.reliability_value == "未知"):
#                 return False
#         return True

class FilterConditionOther(FilterConditionBase):
    """loss1/gain1，外显子个数<=2的过滤条件"""
    def __init__(self, unbalance_file):
        self.unbalance_file = unbalance_file
        self.unbalance_string = self.get_unbalance_info(unbalance_file)
        self.is_normal_cutoff = (self.get_unbalance_value(unbalance_file) <= 0.01)
        self.is_normal = (self.unbalance_string == "样本正常")
    
    def get_unbalance_info(self, unbalance_file):
        data = load_file(unbalance_file, '\t', '#')
        for line in data:
            try:
                if 'degradation' in line[0]:
                    return line[1].split('|')[0]
            except:
                pass
        return ''

    def get_unbalance_value(self, unbalance_file):
        data = load_file(unbalance_file, '\t', '#')
        for line in data:
            try:
                if 'degradation' in line[0]:
                    return float(line[1].split('(')[1].split(')')[0])
            except:
                pass
        return float("1")


    def get_z_score(self, conv):
        z_score_list = []
        for ex in conv.exons:
            if not ex.contrast.std:
                z_score_list.append(0.0)
            else:
                z_score_list.append((ex.sample_value - ex.contrast.mean) / ex.contrast.std)
        return z_score_list
    
    def check(self, conv):
        length = len(conv.exons)
        cv_list = [ex.contrast.cv for ex in conv.exons]
        mean_cv = np.mean(cv_list)
        
        
        if length == 2 and (conv.tag == "loss1" or conv.tag == "gain1"):
            if conv.reliability_value == "特别可靠" or conv.reliability_value == "相对可靠":
                return True
            if self.is_normal and conv.tag == "loss1" and conv.diff >= 0.4 and conv.diff <= 0.6 and mean_cv <= 0.15:
                return True
            if self.is_normal and conv.tag == "gain1" and conv.diff >= 1.4 and conv.diff <= 1.6 and mean_cv <= 0.15:
                return True
            return False
        if length == 1 and (conv.tag == "loss1" or conv.tag == "gain1"):
            if conv.reliability_value == "特别可靠" or conv.reliability_value == "相对可靠":
                return True
            
            z_score_list = self.get_z_score(conv)

            if self.is_normal_cutoff and conv.tag == "loss1" and conv.diff >= 0.4 and conv.diff <= 0.6 and mean_cv <= 0.15 and abs(float(z_score_list[0])) > 3 and conv.freq <= 0.01:
                return True
            
            if self.is_normal_cutoff and conv.tag == "gain1" and conv.diff >= 1.4 and conv.diff <= 1.6 and mean_cv <= 0.15 and abs(float(z_score_list[0])) > 3 and conv.freq <= 0.01:
                return True 
            return False
        return True

## 白名单不再无条件保留，需要diff值和meanzacore满足条件
class FilterConditionOther2(FilterConditionList):
    def get_mean_z_score(self, conv):
        z_score_list = []
        for ex in conv.exons:
            if ex.contrast.flag_bad_capture == "True":
                pass
            if not ex.contrast.std:
                z_score_list.append(0.0)
            else:
                z_score_list.append(abs((ex.sample_value - ex.contrast.mean) / ex.contrast.std))
        return z_score_list

    def get_value(self, conv):
        num = len(conv.exons)

        intersect_cv = [float(ex.contrast.cv) for ex in conv.exons if ex.contrast.flag_bad_capture != "True"]
        intersect_diff = [float(ex.diff) for ex in conv.exons if ex.contrast.flag_bad_capture != "True" ]
        # intersect_z_score = [abs(float(ex.z_score)) for ex in conv.exons if ex.contrast.flag_bad_capture != "True"]
        intersect_z_score = self.get_mean_z_score(conv)
        intersect_flag_outlier = [ex.flag_outlier for ex in conv.exons if ex.contrast.flag_bad_capture != "True"]

        mean_diff = round(np.mean(intersect_diff), 4)
        mean_z_score = round(np.mean(intersect_z_score), 4)
        cvratio = len(list(filter(lambda x:x<0.2, intersect_cv)))/num
        flag_outlier_ratio = len(list(filter(lambda x:x=="True", intersect_flag_outlier)))/num
        return mean_diff, mean_z_score, cvratio, flag_outlier_ratio


    def flag_value(self, conv, interval_diff=[None, None], interval_zscore=[None, None], interval_cvratio=[None, None], interval_outlierratio=[None, None]):
        
        mean_diff, mean_z_score, cvratio, flag_outlier_ratio = self.get_value(conv)
        
        def compare(value, interval_value):
            interval_value1, interval_value2 = interval_value
            if not interval_value1:
                interval_value1 = -float("inf")
            if not interval_value2:
                interval_value2 = float("inf")
            if value >= interval_value1 and value <= interval_value2:
                return True
            return False

        return compare(mean_diff, interval_diff) and compare(mean_z_score, interval_zscore) and compare(cvratio, interval_cvratio) and compare(flag_outlier_ratio, interval_outlierratio)

    def check(self, conv):
        if conv.reliability_value in ['特别可靠', '相对可靠']:
            return True
        if 'vaf支持' in conv.reliability_show or 'wgs支持' in conv.reliability_show:
            return True
        
        num = len(conv.exons)
        
        flag_white_list = False
        if self.list_convs.find_one(conv, action='intersect'):
            flag_white_list = True
        
        if conv.tag == "loss1":
            f1 = self.flag_value(conv, interval_diff=[0.4, 0.7], interval_zscore=[1.5, None], interval_cvratio=[0.5, None], interval_outlierratio=[None, None])
            f2 = self.flag_value(conv, interval_diff=[0.4, 0.7], interval_zscore=[1.5, None], interval_cvratio=[None, None], interval_outlierratio=[0.5, None])
            f3 = self.flag_value(conv, interval_diff=[0.4, 0.7], interval_zscore=[2, None], interval_cvratio=[None, None], interval_outlierratio=[None, None])
            f4 = self.flag_value(conv, interval_diff=[0.4, 0.7], interval_zscore=[1.5, None], interval_cvratio=[None, None], interval_outlierratio=[None, None])

            if num >= 3 and flag_white_list:
                return f4
            return f1 or f2 or f3


        elif conv.tag == "loss" or conv.tag == "loss2":
            # if flag_white_list:
            #     return self.flag_value(conv, interval_diff=[None, 0.15], interval_zscore=[2, None], interval_cvratio=[None, None], interval_outlierratio=[None, None])
            # return self.flag_value(conv, interval_diff=[None, 0.1], interval_zscore=[2, None], interval_cvratio=[None, None], interval_outlierratio=[None, None])  ## 暂时不考虑离群比例
            return True

        elif conv.tag == "gain1":
            f1 = self.flag_value(conv, interval_diff=[1.3, 1.7], interval_zscore=[2, None], interval_cvratio=[None, None], interval_outlierratio=[None, None])
            f2 = self.flag_value(conv, interval_diff=[1.3, 1.7], interval_zscore=[1.5, None], interval_cvratio=[0.5, None], interval_outlierratio=[None, None])
            f3 = self.flag_value(conv, interval_diff=[1.3, 1.7], interval_zscore=[1.5, None], interval_cvratio=[None, None], interval_outlierratio=[0.5, None])
            
            if num >= 3:
                return f1 or f2 or f3
            return f1

        else:
            if flag_white_list:
                return self.flag_value(conv, interval_diff=[1.7, None], interval_zscore=[2, None], interval_cvratio=[None, None], interval_outlierratio=[None, None])
            return self.flag_value(conv, interval_diff=[1.7, None], interval_zscore=[2, None], interval_cvratio=[None, None], interval_outlierratio=[None, None]) ## 暂时不考虑离群比例


class FilterConditionRepeat(FilterConditionBase):
    def __init__(self, repeat_file):
        self.rp_conv = IO_pickle.load_pickle(repeat_file)
    
    def check(self, conv):
        if len(conv.exons) > 2:
            # if conv.start == 104160062 :
            #     print("Repeatyes")
            return True
        if conv.tag != 'gain1':
            return True
        intervals = self.rp_conv.find_all(conv, action='intersect')
        if intervals:
            max_rp = max(map(attrgetter('repeat_num'), intervals))
            if max_rp >= 3:
                return False
        return True

class FilterConditionFreq(FilterConditionBase):
    """频率过滤条件"""
    def check(self, conv):
        length = len(conv.exons)
        if length > 3:  ## 20211009 by housy
            if conv.freq > 0.1:
                return False
        if length <= 3:
            if conv.all_freq > 0.1:
                return False
        # if conv.start == 104160062 :
        #     print("Freqyes")
        return True

class FilterConditionCnvLength(FilterConditionBase):
    """2、大片段CNV过滤 ：
        >=100k  <500k ，报特别可靠，相对可靠
        >=500k，全部报出"""
    def check(self, conv):
        length = len(conv)
        if length >= 500000:
            return True
        elif length >= 100000:
            if conv.reliability_value == '特别可靠' or conv.reliability_value == '相对可靠':
                return True
        return False

class FilterConditionUnbalance(FilterConditionBase):
    """不均衡过滤
    严重不均衡  is_cnv=True: length >= 500000 and reliability_value == '特别可靠' True else False
                is_cnv=False: False
    轻微不均衡 length>=3 and 可靠 True
    """
    def __init__(self, unbalance_file, is_cnv=False):
        self.unbalance_file = unbalance_file
        self.unbalance_string = self.get_unbalance_info(unbalance_file)
        self.is_severity = (self.unbalance_string == '严重不均衡')
        self.is_slight = (self.unbalance_string == '轻微不均衡')
        self.is_cnv = is_cnv
    
    def get_unbalance_info(self, unbalance_file):
        data = load_file(unbalance_file, '\t', '#')
        for line in data:
            try:
                if 'degradation' in line[0]:
                    return line[1].split('|')[0]
            except:
                pass
        return ''
    
    def check(self, conv):
        if self.is_severity:
            if conv.tag == "loss" or conv.tag == "loss2":  ###### by wangxf 20230412 轻微不均衡报出可靠的loss2
               if conv.reliability_value == '特别可靠' or conv.reliability_value == '相对可靠':
                    return True
            if not self.is_cnv:
                return False
            length = len(conv)
            if length >= 500000 and conv.reliability_value == '特别可靠':
                return True
            else:
                return False
        if self.is_slight:
            # if conv.start == 104160062 :
            #     print(self.is_slight, conv.reliability_value, len(conv.exons))
            if conv.tag == "loss" or conv.tag == "loss2":  ###### by wangxf 20230412 轻微不均衡报出可靠的loss2
               if conv.reliability_value == '特别可靠' or conv.reliability_value == '相对可靠':
                    return True
            if len(conv.exons) >= 3:
                if conv.reliability_value == '特别可靠':
                    return True
                if conv.reliability_value == '相对可靠':
                    if '满足外显子个数条件' in conv.reliability_show:  ##by wangxf 对于轻微不均衡样本,报出的结果需同时满足外显子个数条件和diff值条件
                        if '满足diff值条件' in conv.reliability_show:
                            #print(conv.start,conv.reliability_show)
                            return True
                        else:
                            return False
                    else:
                        return False                                 ##by wangxf 20230412 发现不满足个数条件反而给了True
                        #return True
                    #if conv.start == 9910775 :
                        #print(conv.reliability_show)
                    #return True
            return False
        return True


def get_report_tag_conv_key(gene_attentions_file, repeat_file, unbalance_file):  ## by housy
    """获取report_tag"""
    cs = ConditionSet()

    ## 基于白名单和conv长度，过滤diff值和meanzscore不满足条件的conv ## 20230105 by housy
    fc_ot2 = FilterConditionOther2(gene_attentions_file)
    cs.add_filter_condition(fc_ot2)

    # # gene白名单1档
    # fc_at = FilterConditionAttentions(gene_attentions_file1)
    # cs.add_attention_condition(fc_at)

    # # gene白名单2档 by housy
    # fc_at_c = FilterConditionAttentions_condition(gene_attentions_file2)
    # cs.add_attention_condition(fc_at_c)

    # repeat
    fc_rp = FilterConditionRepeat(repeat_file)
    cs.add_filter_condition(fc_rp)
    # 不均衡
    # fc_ub = FilterConditionUnbalance(unbalance_file, is_cnv=False)
    # cs.add_filter_condition(fc_ub)
    # 频率
    fc_freq = FilterConditionFreq()
    cs.add_filter_condition(fc_freq)

    # 其它条件（loss1/gain1，过滤可靠性为 不可靠/可能不可靠/未知）## 20211009 by housy
    fc_ot = FilterConditionOther(unbalance_file)
    cs.add_filter_condition(fc_ot)
    
    return cs.check

def get_report_tag_gene_key(gene_attentions_file1, gene_attentions_file2, gene_filterfalse_file, unbalance_file):
    """gene额外过滤条件"""
    cs = ConditionSet()
    
    fc_at = FilterConditionAttentions(gene_attentions_file1)
    cs.add_attention_condition(fc_at)
    # gene白名单2档 by housy
    fc_at_c = FilterConditionAttentions_condition(gene_attentions_file2)
    cs.add_attention_condition(fc_at_c)

    # gene黑名单
    fc_ff = FilterConditionListFilterFalse(gene_filterfalse_file)
    cs.add_filter_condition(fc_ff.check_gene, model='function')
    
    # 不均衡
    fc_ub = FilterConditionUnbalance(unbalance_file, is_cnv=False)
    cs.add_filter_condition(fc_ub)
    
    return cs.check

def get_report_tag_cnv_key(gene_filterfalse_file, unbalance_file):
    """gene额外过滤条件"""
    cs = ConditionSet()
    # 
    # gene黑名单
    fc_ff = FilterConditionListFilterFalse(gene_filterfalse_file)
    cs.add_filter_condition(fc_ff)
    # cnv长度
    fc_cl = FilterConditionCnvLength()
    cs.add_filter_condition(fc_cl)
    # 不均衡
    fc_ub = FilterConditionUnbalance(unbalance_file, is_cnv=True)
    cs.add_filter_condition(fc_ub)
    # 频率
    fc_freq = FilterConditionFreq()
    cs.add_filter_condition(fc_freq)
    return cs.check
