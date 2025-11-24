#-*-coding:utf-8-*-
"""
CnvHeteCheck
    输入 
        cnv_data    /share/chg2master/prod/sample_batch/b5344/AAU161_NTb01_NT01T_b5344/5.3_Batch_EXON_deletion/exon_result_1.txt_CNV
        repeat_data    /share/chg2master/EXdev/WES_pipe_v4/db/C_mut_anno_database/repeat.txt_chr_sort
        vcf_data    /share/chg2master/prod/sample_batch/b5344/AAU161_NTb01_NT01T_b5344/acmg2.0_3.1_variation_calling/All_variation.vcf.annoformat
    接口
        get_flag    (WesCnvReliability/core/get_flag_vaf.py)
            input cnv_obj
            return Ture?False?None   
        get_flag2
            input *cnv_item, tag
            return check_res, reason, comments
        get_exon_repeat
            input cnv_item/*cnv_item
            return ave_repeat_num
        check output
            return 文件 *line, check_res, reason, comments

flag_show = {
"": '0. other',
"1.1": "1.1. len<=100k",
"1.2": "1.2. 对应标签tag非gain1或loss1",
"1.3.1": "1.3.1. tag=loss1, vcf_num<10",
"1.3.2": "1.3.2. tag=gain1, vcf_num<10", 
"1.4": "1.4. repeat>1 且 tag对应特征值char_value为 False",
"1.5": "1.5. tag对应特征值char_value为 None", 
"2.1": "2.1. 总体偏离标签deviation_flag为True, tag对应特征值char_value flag为 False",
"3.1": "3.1. 总体偏离标签deviation_flag为True, tag对应特征值char_value flag为 True",
}
"""
import os
import re
import logging
import copy
from operator import *
from itertools import repeat
import fire
import numpy as np
from pybedtools import BedTool

from mod_tools.mod_data_IO import *
from intervaltools.core import *
from collections import Counter
import exon_pipeline.parameters as opt
from exon_pipeline.utils.database_store import DatabaseStore
from exon_pipeline.utils.flags import FlagInfoBase

logger = logging.getLogger(__name__)

vcf_itemgetter = opt.vcf_itemgetter
vcf_regex = opt.vcf_regex

# flag info
class VafFlagInfo(FlagInfoBase):
    def __init__(self, classify, show_string):
        flag_type = 'vaf_flag'
        super().__init__(flag_type, classify, show_string)
    
    def add_comments(self, comments):
        new_obj = copy.deepcopy(self)
        new_obj.comments = comments
        return new_obj

vaf_flag_0 = VafFlagInfo('0', 'other')
vaf_flag_1_1 = VafFlagInfo('1.1', 'len<=10k')
vaf_flag_1_2 = VafFlagInfo('1.2', '对应标签tag非gain1或loss1')
vaf_flag_1_3_1 = VafFlagInfo('1.3.1', 'tag=loss1, vcf_num<10')
vaf_flag_1_3_2 = VafFlagInfo('1.3.2', 'tag=gain1, vcf_num<10')
vaf_flag_1_4 = VafFlagInfo('1.4', 'repeat>1 且 tag对应特征值char_value flag为 False')
vaf_flag_1_5 = VafFlagInfo('1.5', 'tag对应特征值char_value为 None')
vaf_flag_2_1 = VafFlagInfo('2.1', '总体偏离标签deviation_flag为True, tag对应特征值char_value flag为 False')
vaf_flag_3_1 = VafFlagInfo('3.1', '总体偏离标签deviation_flag为True, tag对应特征值char_value flag为 True')


class IntervalTag(IntervalBase):
    def __init__(self, chr, start, end, tag):
        super().__init__(chr, start, end)
        self.tag = tag
    
    def not_include(self, other):
        """self 不包含 other"""
        return not self.include(other)

class IntervalRepeat(IntervalBase):
    def __init__(self, chr, start, end, repeat_num):
        super().__init__(chr, start, end)
        self.repeat_num = self.format_factory(int)(repeat_num)

class CNVInterval(IntervalBase):
    """
    输入 cnvs 类
    """
    def __init__(self, exons, adapted_methods=None):
        """We set the adapted methods in the object's dict"""
        self._obj = exons
        if adapted_methods:
            self.__dict__.update(adapted_methods)
        self.lenght = len(exons)
    
    def __getattr__(self, attr):
        """All non-adapted calls are passed to the object"""
        return getattr(self._obj, attr)
    
    def flag_len(self):
        """1.1"""
        lenght = self.lenght
        return lenght > 10000
    
    def get_type_value(self):
        """获得类型"""
        cnv_type = self.tag
        return cnv_type
    
    def flag_type(self):
        """1.2"""
        cnv_type = self.get_type_value()
        return cnv_type in ['loss1', 'gain1']


class VcfBase(IntervalBase):
    def __init__(self, chr, start, end, vcf_info):
        super().__init__(chr, start, end)
        #self.regex = re.compile('.+/.+=(.+?);(.+?);')
        #self.regex = re.compile('.+/.+=(.+)()')
        self.regex = vcf_regex
        self.vcf_info = vcf_info
        mut_rate, vcf_tag = self.format_vcf_info(vcf_info)
        self.mut_rate = mut_rate
        self.vcf_tag = vcf_tag 
    
    def format_vcf_info(self, vcf_info):
        """从vcf获得突变率"""
        try:
            mut_rate, vcf_tag = self.regex.search(vcf_info).groups()
            mut_rate = float(mut_rate)
        except (AttributeError, IndexError):
            raise TypeError(vcf_info)
        return mut_rate, vcf_tag



def load_sample_vcf(wkcode, vcf_path, interval_path):
    """获取样本数据"""
    vcf_data = load_file(vcf_path, '\t', '#')
    vcf = BucketIndexIntervalList(vcf_data, default_base=VcfBase, key=itemgetter(0,1,2,8))
    
    ## 排除指定区间vcf
    interval_data = load_file(interval_path, '\t', '#')
    interval_obj = BucketIndexIntervalList(interval_data, default_base=CNVInterval, key=itemgetter(0,1,2,4))
    interval_obj.make_index(10000)
    filter_vcf = [mut for mut in vcf._data if not interval_obj.find_one(mut, action='include')]
    vcf._data = filter_vcf
    vcf.make_index(10000)

    ## 添加wkcode信息
    vcf.wkcode = wkcode
    vcf.add_attr('wkcode', repeat(wkcode))
    return vcf

class ControlStatistics:
    """对照库统计
    输入 
        vcfs        样本vcf
        db_vcf_info 数据库样本信息
    """
    def __init__(self, vcfs, db_vcf_info, ):
        self.vcfs = vcfs
        self.db_vcf_info = db_vcf_info
        self.high_freq_homo_interval = self.get_high_freq_homo()
    
    def get_high_freq_homo(self):
        """获取高频纯合mut"""
        c = Counter()
        for muts in self.db_vcf_info:
            #homo_muts = {mut.interval() for mut in muts if 'Homo' in mut.vcf_info}
            homo_muts = {mut.interval() for mut in muts if mut.mut_rate >= 0.85}
            c.update(homo_muts)
        high_freq_ratio = 0.5
        limit = len(self.db_vcf_info) * high_freq_ratio
        junk = {intv for intv, num in c.items() if num >= limit}
        return junk
        
    def get_characteristic_value_loss1(self, mut_rates, mut_counts=10, cutoff1=0.15, cutoff2=None, mut_predict = 1):
        """获取loss1的特征值
        mut_rates type: list(float) or np.array
        return type: np.float64
        """
        mut_rates = np.array(mut_rates)
        mut_rates = mut_rates[mut_rates > cutoff1]
        
        if len(mut_rates) >= mut_counts and len(mut_rates) > 0:
            mse = sum([(x - mut_predict) ** 2 for x in mut_rates])/len(mut_rates)
            return mse
        else:
            return None
    
    def get_characteristic_value_gain1(self, mut_rates, mut_counts=10, cutoff1=0.15, cutoff2=0.85, mut_predict1=0.33, mut_predict2=0.67):
        """获取gain1的特征值
        mut_rates type: list(float) or np.array
        return type: np.float64
        """
        mut_rates = np.array(mut_rates)
        mut_rates = mut_rates[(mut_rates > cutoff1) & (mut_rates < cutoff2)]
        
        if len(mut_rates) >= mut_counts and len(mut_rates) > 0:
            mse = (sum([(x - mut_predict1) ** 2 for x in mut_rates[mut_rates <= 0.5]]) + sum([(x - mut_predict2) ** 2 for x in mut_rates[mut_rates > 0.5]]))/len(mut_rates)
            return mse
        else:
            return None
    
    def get_char(self, vcfs, tag='loss1'):
        """获取特征值"""
        junk = self.high_freq_homo_interval
        mut_rates = [i.mut_rate for i in vcfs if i.interval() not in junk]
        func = getattr(self, "get_characteristic_value_%s"%tag)
        return func(mut_rates)
    
    @property
    def loss_char(self,):
        if not hasattr(self, '_loss_char'):
            self._loss_char = self.get_char(self.vcfs, 'loss1')
        return self._loss_char
    
    @property
    def gain_char(self,):
        if not hasattr(self, '_gian_char'):
            self._gian_char = self.get_char(self.vcfs, 'gain1')
        return self._gian_char
    
    def condition_loss_hete(self, vcfs):
        """loss1的杂合突变比例应当小于50%"""
        junk = self.high_freq_homo_interval
        mut_rates = [i.mut_rate for i in vcfs if i.interval() not in junk]
        total = [i for i in mut_rates if i>0.15]
        hetes = [i for i in total if i<0.8]
        try:
            ratio = len(hetes) / len(total)
            flag = ratio < 0.5
        except:
            flag = True
        return flag
    
    def flag_loss_control(self):
        """正常人q%置信区间q%
        如果杂合比例过高，则返回False
        小于 100 置信区间 支持  True
        大于 90 置信区间 否定   False
        对照库数据不够/其他 不支持不反对   None
        """
        if not self.loss_char:
            return None
        if not self.condition_loss_hete(self.vcfs):
            return False
        r = [self.get_char(vcfs, 'loss1') for vcfs in self.db_vcf_info]
        # 过滤突变个数不够的
        r = [i for i in r if i]
        if len(r) < 10:
            return None
        q100 = np.percentile(r, 0)
        q90 = np.percentile(r, 10)
        if self.loss_char < q100:
            return True
        elif self.loss_char > q90:
            return False
        else:
            None
    
    def flag_gain_control(self):
        """正常人q%置信区间q%
        小于 100 置信区间 支持  True
        大于 90 置信区间 否定   False
        对照库数据不够/其他 不支持不反对   None
        """
        if self.flag_loss_control():
            return False
        if not self.gain_char:
            return None
        r = [self.get_char(vcfs, 'gain1') for vcfs in self.db_vcf_info]
        # 过滤突变个数不够的
        r = [i for i in r if i]
        if len(r) < 10:
            return None
        q100 = np.percentile(r, 0)
        q90 = np.percentile(r, 10)
        if self.gain_char < q100:
            return True
        elif self.gain_char > q90:  ## 20211009 by housy
            return False
        else:
            None
    
    def get_characteristic_comments(self, tag='loss1'):
        """特征值comment"""
        char_value = getattr(self, '%s_char'%tag[:4])
        return '%s_char_value=%s'%(tag, char_value)
    
    def get_flag(self, tag='loss1'):
        """获取特征判定"""
        return getattr(self, 'flag_%s_control'%tag[:4])()


class FlagExonVafBase:
    """
    通过突变分布验证cnv的loss1和gain1信息
    
    对于杂合缺失重复，区段内点突变呈现特殊分布，使用该分布检测cnv杂合缺失重复是否可靠
    
    结果矫正判定 见 CnvHeteCheck.core.__doc__
    
    输入 
        repeat_data    /share/chg2master/EXdev/WES_pipe_v4/db/C_mut_anno_database/repeat.txt_chr_sort
        vcf_data    acmg2.0_3.1_variation_calling/All_variation.vcf.annoformat
    接口
        get_flag    (WesCnvReliability/core/get_flag_vaf.py)
            input cnv_obj
            return Ture?False?None   
        get_flag2
            input *cnv_item, tag
            return check_res, reason, comments
        get_exon_repeat
            input cnv_item/*cnv_item
            return ave_repeat_num
        check output
            return 文件 *line, check_res, reason, comments
    """
    def __init__(self, repeat_obj=None, vcf_obj=None, vcf_db=None):
        self.repeat_obj = repeat_obj
        self.vcf_obj = vcf_obj
        self.vcf_db = vcf_db
        
        self.check_show_str = {None: '-',
                               False: 'no',
                               True: 'yes',
                               }
        self.flag_show = {
                        "": '0. other',
                        "1.1": "1.1. len<=10k",
                        "1.2": "1.2. 对应标签tag非gain1或loss1",
                        "1.3.1": "1.3.1. tag=loss1, vcf_num<10",
                        "1.3.2": "1.3.2. tag=gain1, vcf_num<10", 
                        "1.4": "1.4. repeat>1 且 tag对应特征值char_value flag为 False",
                        "1.5": "1.5. tag对应特征值char_value为 None", 
                        "2.1": "2.1. 总体偏离标签deviation_flag为True, tag对应特征值char_value flag为 False",
                        "3.1": "3.1. 总体偏离标签deviation_flag为True, tag对应特征值char_value flag为 True",
                        }
    
    def get_vcf_infos(self, cnv):
        """输入cnv, 返回list(vcf)"""
        vcfs = self.vcf_obj.find_all(cnv, action='be_included')
        return vcfs
    
    def get_db_vcf_info(self, cnv):
        """输入cnv, 返回对照库的vcf info :list(list(vcf))"""
        vcf_db = self.vcf_db
        return vcf_db.get_data(*cnv.interval())
    
    def get_flag(self, conv):
        """判定逻辑, 返回判断结果和理由"""
        cnv = CNVInterval(conv)
        cnv_vcfs = self.get_vcf_infos(cnv)
        db_vcfs = self.get_db_vcf_info(cnv)
        cs_obj = ControlStatistics(cnv_vcfs, db_vcfs)
        
        flag_len = cnv.flag_len()
        tag = cnv.get_type_value()
        flag_type = cnv.flag_type()
        
        # 长度
        if not flag_len: #
            return vaf_flag_1_1
        
        #类型检测
        if not flag_type:
            return vaf_flag_1_2
        
        # vaf
        # 突变个数
        mut_num = len(cnv_vcfs)
        if mut_num < 10:
            if tag == 'loss1':
                return vaf_flag_1_3_1
            else:
                return vaf_flag_1_3_2
        
        # vaf
        vaf_flag = cs_obj.get_flag(tag=tag)
        comments = cs_obj.get_characteristic_comments(tag=tag)
        if vaf_flag is True:
            return vaf_flag_3_1.add_comments(comments)
        elif vaf_flag is False:
            repeat_flag = self.get_repeat_flag(cnv)
            if repeat_flag:
                return vaf_flag_1_4.add_comments(comments)
            else:
                return vaf_flag_2_1.add_comments(comments)
        else:
            return vaf_flag_1_5.add_comments(comments)

    
    def get_repeat_flag(self, cnv):
        """计算repeat max"""
        repeat_obj = self.repeat_obj
        rp_intervals = repeat_obj.find_all(cnv, action='intersect')
        if rp_intervals:
            rp_intervals = max([i_obj.repeat_num for i_obj in rp_intervals])
        else:
            rp_intervals = 1
        return rp_intervals





def filter_vcf(vcf_path, bed_path, output):
    """指定bed文件 过滤vcf突变位点，保留在该区间内的位点"""
    bed = BedTool(bed_path)
    vcf = BedTool(vcf_path)
    bed_in_vcf = vcf.intersect(bed, output=output, header=True)



def fmt_load_vcf_obj(vcf_path, bed_path=None, output_fmt_vcf=None):
    """读取vcf_data，如果存在bed_path和output_fmt_vcf，则将按bed过滤"""
    if bed_path is None or output_fmt_vcf is None:
        vcf_data = load_file(vcf_path, '\t', '#')
    else:
        # 输出过滤后的vcf
        try:
            filter_vcf(vcf_path, bed_path, output_fmt_vcf)
            vcf_data = load_file(output_fmt_vcf, '\t', '#')
        except:
            vcf_data = load_file(vcf_path, '\t', '#')
    vcf_obj = BucketIndexIntervalList(vcf_data, default_base=VcfBase, key=vcf_itemgetter)
    vcf_obj.make_index(10000)
    return vcf_obj



def get_cnv_check_model(exons_conv_path, vcf_path, repeat_path, db_save_path, output, bed_path=None, output_fmt_vcf=None):
    """构建cnv矫正模型
    ; --exons-conv-path  
    ; --vcf-path        vcf path: All_variation.vcf.annoformat
    ; --repeat-path     repeat_path: repeat.pic
    ; --output          output result file
    ; --bed-path        bed file (filter vcf)
    ; --output-fmt-vcf  format filter vcf output
    ; -- --help         help
    """
    logger.info('loading...')
    convs = IO_pickle.load_pickle(exons_conv_path)
    logger.info('loading2...')
    repeat_data = IO_pickle.load_pickle(repeat_path)
    logger.info('run...')
    vcf_data = fmt_load_vcf_obj(vcf_path, bed_path=bed_path, output_fmt_vcf=output_fmt_vcf)
    ds = ExonDatabaseStore()
    ds.load(db_save_path)
    cnv_check = FlagExonVafBase(repeat_obj=repeat_data, vcf_obj=vcf_data, vcf_db=ds)
    return cnv_check


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    fire.Fire(get_cnv_check_model)
