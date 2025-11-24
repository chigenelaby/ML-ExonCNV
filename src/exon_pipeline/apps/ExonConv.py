#-*-coding:utf-8-*-

"""
# ExonConv
外显子连接算法 双指针标记算法
"""
import os
import copy
import logging
import itertools
import functools
from operator import itemgetter
from collections import OrderedDict
import numpy as np
import pandas as pd

from mod_tools.mod_data_IO import *
from mod_tools.tools_time import wrapper_time

from ExonConv import __version__ as version
from ExonConv.utils.load_qc import get_sexuality_from_qc_report, QCFormatError

logger = logging.getLogger(__name__)

#最大显示个数
MAXSIZE = 20
def get_is_autosomal_region(exon):
    """
    拟常染色体区判断
    chrX    60,001  2,699,520
    chrY    10,001  2,649,520
    chrX    154,931,044 155,260,560
    chrY    59,034,050  59,363,566
    """
    if exon.chr == 'chrX' and (exon.end <= 2699520 or exon.start >= 154931044):
        return True
    if exon.chr == 'chrY' and (exon.end <= 2649520 or exon.start >= 59034050):
        return True
    return False

@functools.lru_cache(maxsize=128)
def get_tag_from_diff(diff_copy_value, is_male_sex_chr=False, is_autosomal_region=False):
    """获取标签
    is_male_sex_chr：男性性染色体阈值
    """
    if not is_male_sex_chr or is_autosomal_region:
        if diff_copy_value >= 1.7:
            return 'gain2'
        elif diff_copy_value >= 1.25:
            return 'gain1'
        elif diff_copy_value <= 0.15:
            return 'loss2'
        elif diff_copy_value <= 0.75:
            return 'loss1'
        else:
            return 'NA'
    else:#男性性染色体 且不是拟常染色体区
        # if diff_copy_value >= 1.65:
        if diff_copy_value >= 1.55:
            return 'gain2'
        elif diff_copy_value <= 0.2:
            return 'loss2'
        else:
            return 'NA'

def is_sex_chr(chr):
    """是否是性染色体"""
    sex_chrs = ['NC_000023.10', 'NC_000024.9', 'NC_000023', 'NC_000024', 'chrX', 'chrY']
    return chr in sex_chrs

class Exon:
    """
    一个外显子
    """
    def __init__(self, chr, start, end, diff_value, sexuality='XX'):
        self.chr = chr
        self.start = int(start)
        self.end = int(end)
        self.diff_value = self.format_diff_value(diff_value)
        self.sexuality = sexuality
        assert sexuality in ['XX', 'XY']
    
    def __len__(self):
        return self.end - self.start + 1
    
    def __str__(self):
        return "Exon(%s, %s, %s)"%(self.chr, self.start, self.end)
    
    def __repr__(self):
        return "Exon(%s, %s, %s)"%(self.chr, self.start, self.end)
    
    def format_diff_value(self, value):
        try:
            value = float(value)
        except (ValueError, TypeError):
            pass
        return value

    def _is_sex_chr(self):
        """是否是性染色体的外显子"""
        return is_sex_chr(self.chr)

    def is_male_sex_chr(self):
        """是否是男性性染色体的外显子"""
        return self.sexuality == 'XY' and self._is_sex_chr()

    def get_tag(self):
        """获取标签"""
        is_male_sex_chr = self.is_male_sex_chr()
        return get_tag_from_diff(self.diff_value, is_male_sex_chr=is_male_sex_chr, is_autosomal_region=get_is_autosomal_region(self))

    def exon(self):
        """外显子信息"""
        return (self.chr, self.start, self.end)

    def cmp_tag(self, tag):
        """比较tag"""
        return self.get_tag() == tag
    
    def tag_type(self, tag=None):
        """
        loss: loss1, loss2
        gain: gain1, gain2
        NA:NA
        """
        if tag is None:
            tag = self.get_tag()
        return tag[:4]
    
    def cmp_tag_type(self, tag):
        """比较tag的类别
        loss: loss1, loss2
        gain: gain1, gain2
        NA:NA
        """
        return self.tag_type() == self.tag_type(tag)
    
    def neg_tag_type(self, tag):
        """否定的类"""
        if self.tag_type() == 'NA' or self.tag_type(tag) == 'NA':
            return False
        else:
            return self.tag_type() != self.tag_type(tag)

    # Comparisons of Exon objects with other.
    def _cmperror(x, y):
        raise TypeError("can't compare '%s' to '%s'" % (
                    type(x).__name__, type(y).__name__))

    def __eq__(self, other):
        if isinstance(other, Exon):
            return self.exon() == other.exon()
        else:
            return False

    def __le__(self, other):
        if isinstance(other, Exon):
            return self.exon() <= other.exon()
        else:
            self._cmperror(self, other)

    def __lt__(self, other):
        if isinstance(other, Exon):
            return self.exon() < other.exon()
        else:
            self._cmperror(self, other)

    def __ge__(self, other):
        if isinstance(other, Exon):
            return self.exon() >= other.exon()
        else:
            self._cmperror(self, other)

    def __gt__(self, other):
        if isinstance(other, Exon):
            return self.exon() > other.exon()
        else:
            self._cmperror(self, other)
    
    #常用函数
    def include(self, other):
        """self包含other"""
        if not isinstance(other, Exon):
            self._cmperror(self, other)
        return self.chr == other.chr and \
               self.start <= other.start and \
               self.end >= other.end
    
    def be_included(self, other):
        """self被other包含"""
        if not isinstance(other, Exon):
            self._cmperror(self, other)
        return other.include(self)


class ExonConvBase:
    """
    使用动态标记处理外显子连接
    
    输入：
    exon_data type: 输入数据，应当已排序，且是不存在包含情况的数据。
    key type: None or function 对exon_data进行map处理，处理后应该满足_format_raw_data的输入条件
    方法：
    _format_raw_data(): 输入函数的格式化，包括junk_sign，都是对输入数据的处理
    next(): 核心方法之一。寻找下一个符合条件的可连接的exons。
    run(): list 返回所有已连接的cnv
    add_condition: 添加条件
    """
    def __init__(self, exon_data, key=None, sexuality="XX"):
        self.tags = ['loss2', 'loss1', 'NA', 'gain1', 'gain2']
        self.junk_sign = 'bad_capture'
        self.condition_types = ['dynamic', 'static_end', 'static_conv'] #动态条件，静态终止条件，静态连接条件
        self.sexuality = sexuality
        assert sexuality in ['XX', 'XY']
        
        self.raw_data, self.raw_data_index = self._format_raw_data(exon_data, key=key)
        
        self.next_conv_index = 0 
        self.dynamic_conditions = [] # 动态条件列表
        self.static_end_conditions = [] # 静态终止条件列表
        self.static_conv_conditions = [] # 静态连接条件列表

    def _format_raw_data(self, raw_data, key=None):
        """格式化原始数据, 输入数据，应当已排序，且是不存在包含情况的数据。
        raw_data 五列： chr, start, end, diff_value, res_sign
        raw_data_index: (chr, start, end): index(包含junk_sign)
        """
        junk_sign = self.junk_sign
        sexuality = self.sexuality
        if key is None:
            key = lambda x: x
        res_data = []
        raw_data_index = OrderedDict()
        for index, (chr, start, end, diff_value, res_sign) in enumerate(map(key, raw_data)):
            ex = Exon(chr, start, end, diff_value, sexuality=sexuality)
            raw_data_index[ex.exon()] = index
            if res_sign == junk_sign:
                continue
            res_data.append(ex)
        return res_data, raw_data_index


    def next(self, next_conv_index=None):
        """获取下一个可连接的exons, 如果没有，返回StopIteration
        type next_conv_index: None or int
        rtype exons: list
        """
        next_conv_index = next_conv_index or self.next_conv_index #
        raw_data = self.raw_data
        dc = self.dynamic_conditions
        sec = self.static_end_conditions
        scc = self.static_conv_conditions
        
        #初始化状态
        last_state = 'N' 
        start_conv_index = None #'C'的开始
        last_conv_index = None #'C'最新标记处
        start_M_index = None    #'M'的开始
        tag = 'NA' #整体tag
        
        # 依次标记
        for index, exon in enumerate(raw_data[next_conv_index:], next_conv_index):
            if last_state == 'N':
                #(1)开始条件
                if self.start_cases(exon):
                    start_conv_index = index #开始索引
                    last_conv_index = index
                    state = 'C'
                    tag = self.get_tag(exon)
                else: #(0) N -> N:     未连接
                    state = 'N'
            elif last_state == 'C' or last_state == 'M':
                exon_tag = self.get_tag(exon)
                # 初始化状态
                state = 'M' 
                # (2)终止条件(静态条件列表,静态终止条件相比而言都是全局性的)
                for condition in sec:
                    flag_sec = not condition(raw_data[start_conv_index: index+1], 
                                             tag=tag,
                                             )
                    if flag_sec:
                        state = 'E'
                        break
                # 判断是否被静态终止条件终止
                if state == 'E':
                    pass
                # 当前exon的tag和tag0不一致 
                elif exon_tag != tag:
                    for condition in dc:
                        # (2)终止条件(动态条件列表)
                        flag = not condition(raw_data[start_conv_index: index+1], tag=tag)
                        if flag:
                            state = 'E'
                            break
                # 当前exon的tag和tag0一致
                else:
                    flags = []
                    # 判断动态连接条件
                    M_index = start_M_index or index #动态连接区间
                    for condition in dc:
                        flag = condition(raw_data[M_index: index+1], tag=tag)
                        flags.append(flag)
                    # 判断静态连接条件
                    for condition in scc:
                        flag = condition(raw_data[start_conv_index: index+1], 
                                         tag=tag, 
                                         )
                        flags.append(flag)
                    conv_flag = all(flags) #连接条件逻辑
                    # (3) C/M -> C:   连接处理
                    if conv_flag:
                        state = 'C'
                        start_M_index = None #消除前面的'M'
                        last_conv_index = index #更新最新连接index
                # (4) C/M -> M:   连接中断 不满足终止条件或连接条件
                if state == 'M': 
                    start_M_index = start_M_index or index#初始化
            if state == 'E':
                break
            # 更新状态
            last_state = state
        #(5)终止处理
        if start_conv_index is None:
            # 如未发现连接则引发StopIteration，表示无下一个cnv
            raise StopIteration
        else:
            exons = raw_data[start_conv_index: last_conv_index+1]
            self.next_conv_index = last_conv_index+1 #更新下一个索引位置
        return exons

    @wrapper_time("程序运行")
    def run(self, next_conv_index=None):
        """执行，连接后数据"""
        res = []
        while True:
            try:
                exons = self.next()
                res.append(exons)
            except StopIteration:
                break
        return res

    def get_tag(self, exon):
        """获取标签
        return: in self.tags
        """
        return exon.get_tag()

    def start_cases(self, exon):
        """
        exon: 当前外显子
        return: bool
        """
        return self.get_tag(exon) != 'NA'

    def add_condition(self, condition_obj, condition_type='static_end'):
        """添加条件的快速接口
        condition_obj type: ConditionBase, 或者是func
        condition_type type: str
        """
        assert condition_type in self.condition_types
        condition_store = getattr(self, '%s_conditions'%condition_type)
        condition_store.append(condition_obj)
        logger.info('add %s condition: %s'%(condition_type, condition_obj))



class ConditionBase:
    """条件基础类
    exons应该是Exon的实例"""
    def synopsis(self, maxsize=MAXSIZE):
        info = self.__doc__.split('\n')[0]
        if len(info) > maxsize:
            info = '%s...'%(info[:maxsize])
        return info
    
    def condition(self, exons, *args, **keyword):
        """条件判断
        type exons: list(Exon实例) 
        rtype: bool
        """
    
    def __call__(self, exons, *args, **keyword):
        """调用函数
        >>> cd = ConditionBase()
        >>> cd(exons)
        """
        return self.condition(exons, *args, **keyword)
    
    def __str__(self):
        return 'Condition(def="%s")'%(self.synopsis())
    
    def __repr__(self):
        return self.__str__()

class DynamicCondition(ConditionBase):
    """动态条件"""
    _theta_correction_ratio = 0.0001 #修正率
    def __init__(self, dynamic_ratio):
        self.dynamic_ratio = dynamic_ratio
    
    def __str__(self):
        return 'Condition(def="%s", ratio="%.2f")'%(self.synopsis(), self.dynamic_ratio)

class StaticCondition1(ConditionBase):
    """静态条件
    仅使用最后的数据进行判断
    """
    def condition(self, exon, *args, **keyword):
        """条件判断
        type exon: Exon实例 
        rtype: bool
        """

    def __call__(self, exons, *args, **keyword):
        """调用函数
        >>> cd = ConditionBase()
        >>> cd(exons)
        """
        return self.condition(exons[-1], *args, **keyword)

class StaticCondition2(ConditionBase):
    """静态条件
    仅使用最后的两数据进行判断
    """
    def condition(self, exon1, exon2, *args, **keyword):
        """条件判断
        type exon: Exon实例 
        rtype: bool
        """

    def __call__(self, exons, *args, **keyword):
        """调用函数
        >>> cd = ConditionBase()
        >>> cd(exons)
        """
        try:
            ex1 = exons[-2]
            ex2 = exons[-1]
        except IndexError:
            return True
        return self.condition(ex1, ex2, *args, **keyword)

class ChrCondition(StaticCondition2):
    """染色体号条件，保持染色体一致"""
    def condition(self, exon1, exon2, *args, **keyword):
        """条件判断
        type exon: Exon实例 
        rtype: bool
        """
        return exon1.chr == exon2.chr


class ExonNumCondition(DynamicCondition):
    """外显子个数比例条件
    主|间    loss2  loss1  NA  gain1  gian2
    loss2    1.0    0.6   0   -0.1   -0.1
    loss1    1.0    1.0   0   -0.1   -0.1
    NA       0.0    0.0   0    0.0    0.0
    gain1   -0.1   -0.1   0    1.0    1.0
    gian2   -0.1   -0.1   0    0.6    1.0
    参照矩阵 外显子个数分值个数比要>=0.5
    """
    
    def __init__(self, dynamic_ratio=0.5):
        self.dynamic_ratio = dynamic_ratio
        data = [[1, 0.6, 0, -0.1, -0.1],
                [1, 1, 0, -0.1, -0.1],
                [0, 0, 0, 0, 0],
                [-0.1, -0.1, 0, 1, 1],
                [-0.1, -0.1, 0, 0.6, 1]]
        index = ['loss2', 'loss1', 'NA', 'gain1', 'gain2']
        self.score_matrix = pd.DataFrame(data=data, index=index, columns=index)

    def get_score(self, ex_tag, tag):
        """ex_tag: 目标exon的tag
        tag: 整体的tag
        rtype: float"""
        return self.score_matrix[ex_tag][tag]

    def condition(self, exons, tag=None, **keyword):
        """条件判断
        type exons: list(Exon实例) 
        rtype: bool
        """
        dynamic_ratio = self.dynamic_ratio
        try:
            all_count = len(exons)
            tag_score = sum(map(lambda ex:self.get_score(ex.get_tag(), tag), exons))#总分
            ratio = tag_score / all_count
        except ZeroDivisionError:
            ratio = 0
        return ratio >= self.dynamic_ratio

def diff_conv_mean(exons):
    """diff_value均值计算"""
    diff_values = list(map(lambda x:x.diff_value, exons))
    mean_diff_value = np.mean(diff_values)
    return mean_diff_value

class DiffConvCondition(ConditionBase):
    """连接diff值条件，连接后整体tag与标记tag一致"""
    def condition(self, exons, tag=None, **keyword):
        """条件判断
        type exons: list(Exon实例) 
        rtype: bool
        """
        if len(exons):
            ex0 = exons[0]
            is_male_sex_chr = ex0.is_male_sex_chr()
            is_autosomal_region = get_is_autosomal_region(ex0)
        else:
            return False
        mean_diff_value = diff_conv_mean(exons)
        dff_tag = get_tag_from_diff(mean_diff_value, is_male_sex_chr=is_male_sex_chr, is_autosomal_region=is_autosomal_region)
        return dff_tag == tag


class ExonConvCNV(ExonConvBase):
    """将exon连接成cnv"""
    def filter_ex_data(self, ex_data):
        """去掉前后外显子包含部分
        ex_data type: list((exon, res_sign))
        """
        res_ex_data = []
        res_ex_data.append(ex_data[0])
        for (ex1, res_sign1), (ex2, res_sign2) in zip(ex_data, ex_data[1:]):
            # ex1包含ex2
            if ex1.include(ex2):
                continue
            # ex2包含ex1
            if ex1.be_included(ex2):
                res_ex_data.pop(-1)
            res_ex_data.append((ex2, res_sign2))
        return res_ex_data
    
    def _format_raw_data(self, raw_data, key=None):
        """
        输入数据，去掉前后外显子包含部分,去掉部分不计入总索引
        格式化原始数据, 输入数据，应当已排序，且是不存在包含情况的数据。
        raw_data 五列： chr, start, end, diff_value, res_sign
        raw_data_index: (chr, start, end): index(包含junk_sign)
        """
        junk_sign = self.junk_sign
        sexuality = self.sexuality
        if key is None:
            key = lambda x: x

        ex_data = []
        for chr, start, end, diff_value, res_sign in map(key, raw_data):
            # chr, start, end, diff_value, res_sign = line
            ex = Exon(chr, start, end, diff_value, sexuality=sexuality)
            ex_data.append((ex, res_sign))
        
        #去掉前后外显子包含部分,去掉部分不计入总索引
        ex_data = self.filter_ex_data(ex_data)
        
        res_data = [] #不包含junk_sign数据
        raw_data_index = OrderedDict() # 包含所有去重数据的index
        for index, (ex, res_sign) in enumerate(ex_data):
            raw_data_index[ex.exon()] = index
            if res_sign == junk_sign:
                continue
            res_data.append(ex)
        return res_data, raw_data_index

    def get_one_cnv_info(self, cnv):
        """
        cnv type: exons, run()的输出list的单元
        rtype: ["#chr", "start", "end", "length", "type", "diff_ave", "exon_num", "exon_num_impact", "exon_num_badcapture"]
        """
        ex_index = self.raw_data_index
        
        chr = cnv[0].chr
        start = cnv[0].start
        end = cnv[-1].end
        length = end - start + 1
        tag = cnv[0].get_tag()
        diff_ave = diff_conv_mean(cnv)
        exon_num_no_junk = len(cnv)
        exon_num = ex_index[cnv[-1].exon()] - ex_index[cnv[0].exon()] + 1
        exon_num_impact = len(list(filter(lambda x: x.cmp_tag_type(tag), cnv)))
        exon_num_badcapture = exon_num - exon_num_no_junk # 总index长度-过滤后长度
        return (chr, start, end, length, tag, diff_ave, exon_num, exon_num_impact, exon_num_badcapture)

    def is_exact_cnv(self, cnv_info):
        """判断是否为准确的cnv
        cnv_info type: get_one_cnv_info()的输出list中的单元
        rtype: bool 
        """
        tolerance_num = 1
        normal_min = 10
        tag = cnv_info[4]
        if tag == 'loss2':
            tolerance_min_exon = 5
        else:
            tolerance_min_exon = 7 
        
        # normal_min或以上直接比较， [tolerance_min_exon, normal_min]容忍个数为tolerance_num
        if cnv_info[7] >= normal_min:
            flag = True
        elif cnv_info[7] >= tolerance_min_exon:
            # 错误长度 = 总长 - 同类tag - badcapture <= 容忍长度
            flag = (cnv_info[6] - cnv_info[7] - cnv_info[8] <= tolerance_num)
        else:
            flag = False
        return flag

    def get_cnvs(self, headers=True):
        """获取所有可拼接cnv"""
        headers_info = ["#chr", "start", "end", "length", "type", "diff_ave", "exon_num", "exon_num_impact", "exon_num_badcapture"]
        version_header = ['## version %s'%version]
        cnvs_data = self.run()
        cnv_infos = map(self.get_one_cnv_info, cnvs_data)
        # exon_num_impact>=10
        cnvs = list(filter(self.is_exact_cnv, cnv_infos))
        if headers:
            cnvs.insert(0, headers_info)
            cnvs.insert(0, version_header)
        return cnvs

class ExonConvCNVInner(ExonConvCNV):
    """用于递归，输入变成了单段的exons和index_map"""
    def __init__(self, exons, exon_index=None):
        super().__init__((exons, exon_index), key=None)
    
    def _format_raw_data(self, raw_data, key=None):
        """直接返回即可"""
        exons, exon_index = raw_data
        return exons, exon_index

def view_conv_tag(cnv):
    """将cnv的exons的tag信息进行展示，按tag分类"""
    view_info = []
    last_tag = None
    for ex in cnv:
        ex_tag = ex.get_tag()
        if ex_tag != last_tag:
            view_info.append([1, ex_tag])
        else:
            view_info[-1][0] += 1
        last_tag = ex_tag
    return view_info

def view_same_tag(cnv, tag):
    """将cnv的exons的tag信息进行展示,按tag相同和不同进行分类"""
    assert tag in ['loss1', 'loss2', 'gain1', 'gain2']
    view_info = []
    last_tag = None
    for ex in cnv:
        ex_tag = ex.get_tag()
        if ex_tag != tag:
            ex_tag = 'NA'
        if ex_tag != last_tag:
            ex_tag = 'NA'
            view_info.append([1, ex_tag])
        else:
            view_info[-1][0] += 1
        last_tag = ex_tag
    return view_info

# 终止条件：连接内部不应该包含准确的cnv
class InnerCnvCondition(ConditionBase):
    """内部cnv条件，内部不应该出现准确的cnv，静止终止条件"""
    def is_exact_cnv(self, num, tag):
        """判断是否为准确的cnv
        rtype: bool 
        """
        if tag == 'loss2':
            min_num = 5
        else:
            min_num = 10
        return num >= min_num
    
    def condition(self, exons, tag=None, **keyword):
        """条件判断
        type exons: list(Exon实例) 
        M_index: int
        rtype: bool
        """
        view_info = view_conv_tag(exons)
        condition_flag = True
        filter_tags = ['loss1', 'loss2', 'gain1', 'gain2']
        try:
            filter_tags.remove(tag)
        except ValueError:
            pass
        for count, ex_tag in view_info:
            if ex_tag in filter_tags and self.is_exact_cnv(count, ex_tag):
                condition_flag = False
        return condition_flag

# 终止条件：连续10个非tag外显子即中断，连续否定tag都为loss2，则为5个。
class NegativeNumCondition(ConditionBase):
    """否定个数条件。连续10个非tag外显子即中断。"""
    def is_negative(self, num):
        """判断是否满足否定个数
        rtype: bool 
        """
        min_num = 10
        return num >= min_num
    
    def condition(self, exons, tag=None, **keyword):
        """条件判断
        type exons: list(Exon实例) 
        M_index: int
        rtype: bool
        """
        view_info_neg = view_same_tag(exons, tag=tag)
        condition_flag = True
        #连续10个非tag外显子即中断
        for count, ex_tag in view_info_neg:
            if ex_tag != tag and self.is_negative(count):
                condition_flag = False
        return condition_flag

# 连接条件：外显子间距大于500k时，单侧长度必须大于间距0.33或单侧外显子个数不小于30连接
class ExonSpacingCnvCondition(ConditionBase):
    """外显子间距条件：外显子间距大于500k时，单侧长度必须大于间距0.33或单侧外显子个数不小于10才能连接或单侧连接exon连续4个tag一致。"""
    continuous_tag_num = 4
    
    def is_continuous_tag(self, exons, tag, is_right_side=True):
        """
        单侧连接exon连续tag数判断
        exons type: 单侧exons
        rtype: bool
        """
        ctn = self.continuous_tag_num
        #右侧判断还是左侧判断
        if is_right_side:
            exons = exons[:ctn]
        else:
            exons = exons[-ctn:]
        
        if len(exons) != ctn:
            flag = False
        else:
            flag = all(map(lambda ex:ex.cmp_tag(tag), exons))
        return flag
    
    def is_ignore(self, exons, spacing_lenght, tag, is_right_side=True):
        """可以忽略的情况 单侧长度必须大于间距0.33或单侧外显子个数不小于10或单侧连接exon连续4个tag
        exons type: 单侧exons
        spacing_lenght type: int 区间长度
        rtype: bool 
        """
        exons_lenght = exons[-1].end - exons[0].start + 1
        continuous_flag = self.is_continuous_tag(exons, tag, is_right_side=True)
        return exons_lenght/spacing_lenght > 0.33 or len(exons) >= 10 or continuous_flag
    
    def condition(self, exons, tag=None, **keyword):
        """条件判断
        type exons: list(Exon实例) 
        M_index: int
        rtype: bool
        """
        condition_flag = True
        for index, (ex1, ex2) in enumerate(zip(exons, exons[1:]), 1):
            spacing_lenght = ex2.start - ex1.end - 1
            if ex1.chr == ex2.chr and spacing_lenght > 500000:
                #间距大于指定值且单侧不能达到忽略条件
                if not self.is_ignore(exons[:index], spacing_lenght, tag, is_right_side=False) or \
                    not self.is_ignore(exons[index:], spacing_lenght, tag, is_right_side=True):
                    condition_flag = False
                    break
        return condition_flag


def exon_conv(exon_file, output, qc_report_path=None):
    """
    通过标记状态和状态的变化，判断外显子是否连接。
    -- --help           帮助
    --exon-file         外显子数据
    --output            输出路径
    --qc-report-path    qc文件，提供性别信息，省略默认为"XX"，男性性染色体不做额外处理
    
    示例：
    $ python -m ExonConv --exon-file ./5.3_Batch_EXON_deletion/exon_result_1.txt.bak --qc-report-path acmg2.0_4.1_QC/QC_report.txt --output ./5.3_Batch_EXON_deletion/exon_result_1.txt_CNV
    """
    # 读数据
    logger.info('loading...')
    exon_data = load_file(exon_file, '\t', '#')
    key = itemgetter(0, 1, 2, 10, 11)
    
    if qc_report_path is None:
        sexuality = 'XX'
    else:
        sexuality = get_sexuality_from_qc_report(qc_report_path)
    
    conv_obj = ExonConvCNV(exon_data, key=key, sexuality=sexuality)
    
    #染色体号条件
    chr_conditon = ChrCondition()
    conv_obj.add_condition(chr_conditon, condition_type='static_end')
    
    #外显子个数条件
    ex_num_conditon = ExonNumCondition(0.5)
    conv_obj.add_condition(ex_num_conditon, condition_type='dynamic')
    
    #diff值条件
    dff_conditon = DiffConvCondition()
    conv_obj.add_condition(dff_conditon, condition_type='static_conv')
    
    # 内部cnv条件，内部不应该出现准确的cnv，静止终止条件
    neg_num_conditon = NegativeNumCondition()
    conv_obj.add_condition(neg_num_conditon, condition_type='static_end')
    
    # 否定个数条件。连续10个非tag外显子即中断.
    inner_cnv_conditon = InnerCnvCondition()
    conv_obj.add_condition(inner_cnv_conditon, condition_type='static_end')
    
    # 外显子间距条件：外显子间距大于500k时，单侧长度必须大于间距0.33或单侧外显子个数不小于10连接
    exon_spacing_conditon = ExonSpacingCnvCondition()
    conv_obj.add_condition(exon_spacing_conditon, condition_type='static_conv')
    
    cnvs = conv_obj.get_cnvs(headers=True)
    output = os.path.abspath(output)
    IO_file.write_file(cnvs, output)
    logger.info('output is "%s"'%output)
