#-*-coding:utf-8-*-

"""草考原外显子可靠性模块
与cnvseq数据结果文件和中间文件对比
cnvseq中间文件 wgs_res_path 2.STAT/compatible.cnv.tt 取0,1,2,3列
cnvseq中间文件 wgs_medi 2.STAT/lib.cnv 取0,1,2,5列
"""
import re
from operator import *

import numpy as np

from mod_tools.mod_data_IO import *
from intervaltools.core import *
from exon_pipeline.core.contrast import *
from exon_pipeline.core.controlled import *
from exon_pipeline.core.conv import *
from exon_pipeline.utils.io import fmt_write_file, show
from exon_pipeline.utils.flags import FlagInfoBase

logger = logging.getLogger(__name__)


# 最短校正限制
LEN_LIMIT = 10000

# 结果文件tag限制比例，交集区段tag小于该比例的不计入相交区段
# 排除同类tag总长占比小于限制比例的类，该参数为类型占比下限
RESULE_FILE_TAG_LIMIT_PROPORTION = 0.1

# wgs中间文件与wes结果区段相交，且满足该相交区域占wgs窗口长度比例θ>=0.7。
INTERMEDIATE_FILE_FILTER_VALUE = 0.7

reg = re.compile('([1-9xyXY]\d*)(?=(\.\d+)*$)')
_chromosome_map = {'23': 'X', '24': 'Y', 'x': 'X', 'y': 'Y'}
wgs_tag_dict = {'1': 'gain1', '2': 'gain2', '-1':'loss1', '-2': 'loss2'}

class CnvSeqFlagInfo(FlagInfoBase):
    def __init__(self, classify, show_string):
        flag_type = 'cnv_seq_flag'
        super().__init__(flag_type, classify, show_string)

##flags
seq_flag_1_0 = CnvSeqFlagInfo('1.1.0', 'CnvSeq结果文件存在，有交集，且tag不同类')
seq_flag_1_1 = CnvSeqFlagInfo('1.1.1', 'CnvSeq结果文件存在，有交集，且tag同类')
seq_flag_1_2 = CnvSeqFlagInfo('1.1.2', 'CnvSeq结果文件存在，有交集，且存在同类tag和不同类tag')
seq_flag_2_1_0 = CnvSeqFlagInfo('1.2.1.0',
    'CnvSeq结果文件无交集，中间文件区段diff值满足强条件，且tag不同类')
seq_flag_2_1_1 = CnvSeqFlagInfo('1.2.1.1',
    'CnvSeq结果文件无交集，中间文件区段diff值满足强条件，且tag同类')
seq_flag_2_2_0 = CnvSeqFlagInfo('1.2.2.0',
    'CnvSeq结果文件无交集，中间文件区段diff值仅满足弱条件，且tag不同类')
seq_flag_2_2_1 = CnvSeqFlagInfo('1.2.2.1',
    'CnvSeq结果文件无交集，中间文件区段diff值满足弱条件，且tag同类')
seq_flag_3_0 = CnvSeqFlagInfo('1.3.0',
    'CnvSeq结果文件无交集，中间文件区段无交集')
seq_flag_3_1 = CnvSeqFlagInfo('1.3.1',
    'CnvSeq结果文件无交集，中间文件区段存在正常位点区间的')
seq_flag_9 = CnvSeqFlagInfo('1.9',
    '其他, CnvSeq结果无证据')


def chr_conv(chr):
    """MC格式转化chr"""
    num = reg.search(chr).group(1)
    fmt_num = _chromosome_map.get(num, num)
    return 'chr%s'%fmt_num


def tag_conv(tag):
    """tag标记转化"""
    return wgs_tag_dict.get(tag, tag)

def get_tag_type(tag):
    """获取tag类型"""
    return tag[:4]



class WgsDiff(IntervalBase):
    """cnv seq 中间文件类"""
    def __init__(self, chr, start, end, diff):
        super().__init__(chr, start, end)
        self.diff = self.format_factory(float)(diff)
        self.chr = chr_conv(self.chr)


class WgsTag(IntervalBase):
    """cnv seq 结果文件类"""
    def __init__(self, chr, start, end, tag):
        super().__init__(chr, start, end)
        self.tag = tag_conv(tag)
        self.chr = chr_conv(self.chr)


class ConditionDiffBase:
    """
    # 5. diff copy值条件
    # 5.1. 强条件
    # 5.1.1. >1.75 gain2
    # 5.1.2. <=1.75 且 >1.4 gain1
    # 5.1.3. <0.15 loss2
    # 5.1.4. >=0.15 且 <0.65 loss1
    # 5.2. 弱条件
    # 5.2.1. <0.8 loss1
    # 5.2.2. >1.2 gain1
    
    # 存在接口strong_condition, weak_condition
    """
    def __init__(self):
        pass
    
    def strong_condition(self, diff_copy_value):
        """强条件标签"""
        if diff_copy_value > 1.75:
            return 'gain2'
        elif diff_copy_value > 1.4:
            return 'gain1'
        elif diff_copy_value < 0.15:
            return 'loss2'
        elif diff_copy_value < 0.65:
            return 'loss1'
        else:
            return ''
    
    def weak_condition(self, diff_copy_value):
        """弱条件标签"""
        if diff_copy_value < 0.8:
            return 'loss1'
        elif diff_copy_value > 1.2:
            return 'gain1'
        else:
            return ''


class GetFlagCnvSeq:
    """
    wgs_res_path snvseq结果文件
    wgs_medi_path snvseq中间文件
    get_flag conv应该的外显子连接的结果，返回flag
    """
    def __init__(self, wgs_res_path, wgs_medi_path, condition_diff_obj=None):
        self.condition_diff_obj = condition_diff_obj or ConditionDiffBase()
        self.wgs_res_obj = self.load_cnv_seq_data(wgs_res_path, base=WgsTag, key=itemgetter(0, 1, 2, 3))
        self.wgs_medi_obj = self.load_cnv_seq_data(wgs_medi_path, base=WgsDiff, key=itemgetter(0, 1, 2, 5))
        self.limit_proportion = RESULE_FILE_TAG_LIMIT_PROPORTION
        self.filter_value = INTERMEDIATE_FILE_FILTER_VALUE
        self.len_limit = LEN_LIMIT
    
    def load_cnv_seq_data(self, path, base=IntervalBase, key=None):
        """读取数据"""
        data = load_file(path, '\t', '#')
        obj = BucketIndexIntervalList(data, default_base=base, key=key)
        obj.make_index(10000)
        return obj
    
    def get_tag_type(self, tag):
        """返回tag类型"""
        return get_tag_type(tag)
    
    def get_seqs_proportion(self, conv, seqs):
        """输入目标区段, 以及相交的seqs，返回相交seqs总相交长度/ 目标总长"""
        total_len = len(conv)
        intersect_len = 0
        for seq in seqs:
            intersect_len += max(min(seq.end, conv.end) - max(seq.start, conv.start), 0)
        return intersect_len / total_len   ## 20211009 by housy
    
    def get_filter_seqs(self, conv, seqs):
        """输入目标区段, 以及相交的中间文件seqs
        取交集比例占seq一定区间的seq
        """
        # 过滤
        filter_seqs = []
        for seq in seqs:
            intersect_len =  max(min(seq.end, conv.end) - max(seq.start, conv.start), 0)
            intersect_proportion = intersect_len / len(seq)
            if intersect_proportion >= self.filter_value:
                filter_seqs.append(seq)
        return filter_seqs
    
    def get_flag(self, conv):
        """一个外显子连接区段返回flag"""
        if len(conv) < self.len_limit:
            return seq_flag_9
        
        tag = conv.tag
        tag_type = self.get_tag_type(tag)
        
        res_seqs = self.wgs_res_obj.find_all(conv, action='intersect')
        # 按tag类分类
        tag_type_seqs_dict = {}
        for seq in res_seqs:
            tag_type_seqs_dict.setdefault(self.get_tag_type(seq.tag), []).append(seq)
        
        # 结果文件比较
        flag_res_support = False
        flag_res_against = False
        for seq_tag_type, seqs in tag_type_seqs_dict.items():
            intersect_proportion = self.get_seqs_proportion(conv, seqs)
            if intersect_proportion > self.limit_proportion:
                if seq_tag_type == tag_type:
                    flag_res_support = True
                elif seq_tag_type != tag_type:  ##20211009 by housy
                    flag_res_against = True
        
        if not flag_res_support and flag_res_against:
            return seq_flag_1_0
        elif flag_res_support and not flag_res_against:
            return seq_flag_1_1
        elif flag_res_support and flag_res_against:
            return seq_flag_1_2
        
        ## 结果文件无结果的，进行中间文件比较，取交集比例占seq一定区间的seq计算diff值均值
        medi_seqs = self.wgs_medi_obj.find_all(conv, action='intersect')
        filter_seqs = self.get_filter_seqs(conv, medi_seqs)
        if not filter_seqs:
            return seq_flag_3_0
        
        mean_diff = np.mean(list(map(attrgetter('diff'), filter_seqs)))
        
        # 强条件配置
        seq_tag1 = self.condition_diff_obj.strong_condition(mean_diff)
        seq_tag_type1 = self.get_tag_type(seq_tag1)
        if seq_tag_type1 == tag_type:
            return seq_flag_2_1_1
        elif seq_tag_type1 and seq_tag_type1 != tag_type:
            return seq_flag_2_1_0
        
        # 弱条件配置
        seq_tag2 = self.condition_diff_obj.weak_condition(mean_diff)
        seq_tag_type2 = self.get_tag_type(seq_tag2)
        if seq_tag_type2 == tag_type:
            return seq_flag_2_2_1
        elif seq_tag_type2 and seq_tag_type1 != tag_type:
            return seq_flag_2_2_0
        elif not seq_tag_type2:
            return seq_flag_3_1
        
        return seq_flag_9

    


