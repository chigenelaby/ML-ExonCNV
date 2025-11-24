#-*-coding:utf-8-*-
"""对照库构建"""


import os
import re
import logging
from operator import *
from collections import *
from itertools import repeat

import numpy as np

from mod_tools import ProgressBar
from mod_tools.mod_data_IO import *
from intervaltools.core import *

import exon_pipeline.parameters as opt
from exon_pipeline.utils.sample_info import get_wkcode, WkCodeBase
from exon_pipeline.utils.database_store import get_index
from exon_pipeline.utils.exon_cls import IntervalDepth, ExonDatabaseStore, IntervalContrast
from exon_pipeline.utils.sample_tools import get_info, read_config
from exon_pipeline.core.model import models_dict, ContrastModel1, ContrastModel2, ContrastModel3
from exon_pipeline.utils.io import fmt_write_file

logger = logging.getLogger(__name__)

# 配置文件文件名
DEPTH_FILE_STR = 'depth_file'
STAT_FILE_STR = 'stat_file'
NT_TARGET_STR = 'nt_target'

# 文件中字段
STAT_DEPTH_STR = opt.STAT_DEPTH_STR 

# 参数
MIN_CRITICAL = 0 #过滤最低均深临界

fun_contrast = opt.fun_contrast


def get_samlpe_depth(wkcode,):
    """
    输入文库号类 wkcode类
    获取样本深度 float
    """
    path = wkcode.get_filename(STAT_FILE_STR)
    return fun_contrast(path, STAT_DEPTH_STR, type_cls=float)


class ContrastBase:
    """批次对照库
    输入item (wkcode, analyse_path, *_)
    
    """
    def __init__(self, contrast_items, 
                 wkcode_cls=WkCodeBase, config_file=None, 
                 config=None, model=ContrastModel1(), model2=ContrastModel3()):
        self.wkcode_cls = wkcode_cls
        self.config = self._config(config_file, config)
        self.contrast_items = contrast_items
        self.wkcodes = [self.get_wkcode(wk_item) for wk_item in contrast_items]
        self.model = model
        self.model2 = model2
    
    def _config(self, config_file, config):
        """获取config"""
        if config_file:
            config = read_config(config_file)
        else:
            config = config
        return config
    
    def get_wkcode(self, wk_item):
        """获取wkcode类"""
        config = self.config
        wkcode, analyse_path, *_ = wk_item
        return get_wkcode(wkcode, analyse_path, cls=self.wkcode_cls, config=config)
    
    def get_index_info(self, wkcode,):
        """读入文件并构建桶索引"""
        path = wkcode.get_filename(DEPTH_FILE_STR)
        return get_index(path, default_base=IntervalDepth, key=itemgetter(0,1,2,6))
    
    def get_samlpe_depth(self, wkcode):
        """获取样本深度 float"""
        return get_samlpe_depth(wkcode)

    def get_database(self):
        """获取对照库数据"""
        ds = ExonDatabaseStore()
        for wkcode in self.wkcodes:
            idx = self.get_index_info(wkcode)
            # 添加文库信息
            idx.wkcode = wkcode
            wk = str(wkcode)
            idx.add_attr('wkcode', values=repeat(wk))
            # 添加样本值信息
            sample_depth = self.get_samlpe_depth(wkcode)
            idx.sample_depth = sample_depth
            # idx.add_attr('sample_value', key=lambda x: x.depth/idx.sample_depth)
            ds.add_sample(idx)
            logger.debug('加载 %s'%wkcode)
        return ds
    
    # 计算
    def get_filter(self, arr, critical=MIN_CRITICAL):
        """critical"""
        return arr[arr > critical]
    
    def get_normal(self, arr):
        """提取正常点"""
        return self.model.get_normal(arr)
  
    def get_normal2(self, arr):
        return self.model2.get_normal(arr)
    
    def get_norm_arr(self, arr, ):
        """重复检验"""
        c = self.get_filter(arr)
        length = len(c)
        while True:
            next_arr = self.get_normal(c)
            if len(next_arr) == length:
                break
            c = next_arr
            length = len(next_arr)
        return c

    def get_norm_arr2(self, arr, ):
        """重复检验"""
        c = self.get_filter(arr)
        length = len(c)
        while True:
            next_arr = self.get_normal2(c)
            if len(next_arr) == length:
                break
            c = next_arr
            length = len(next_arr)
        return c
 
    def get_arr_sts(self, arr):
        """arr从获取统计值"""
        norm_arr2 = self.get_norm_arr2(arr)
        size = len(norm_arr2)
        if len(norm_arr2) >= 10:
            norm_arr = norm_arr2
        else: 
            norm_arr = self.get_norm_arr(arr)
        size = len(norm_arr)
        if size:
            mean = norm_arr.mean()
        else:
            mean = 0
        if mean == 0:
            std = 0
            cv = 0
            min_critical = 0
            max_critical = 0
        else:
            std = norm_arr.std(ddof=1)
            cv = std / mean
            min_critical = norm_arr.min()
            max_critical = norm_arr.max()
        return mean, cv, std, size, min_critical, max_critical

    def get_contrast(self, chr, start, end):
        """获取对照库信息"""
        ds = self.ds
        exons = ds.get_data(chr, start, end)
        if any(i is None for i in exons):    #yzh
            return None     
        dps = np.array([i.depth for i in exons])
        sample_depths = np.array([idx.sample_depth for idx in ds])
        # print(chr, start, end)
        # print(dps, sample_depths)
        arr = dps / sample_depths
        mean, cv, std, size, min_critical, max_critical = self.get_arr_sts(arr)
        norm_dps = dps[np.where((arr >= min_critical) & (arr <= max_critical))]
        dps_mean = norm_dps.mean() if len(norm_dps) else 0
        exon_contrast =  IntervalContrast(chr, start, end, 
                            mean=mean, cv=cv, std=std, 
                            size=size, 
                            min_critical=min_critical, max_critical=max_critical, dps_mean=dps_mean)
        return exon_contrast
    
    def forwark(self,):
        """执行"""
        logger.info('加载数据')
        ds = self.get_database()
        self.ds = ds
    
    def run(self, is_show=False):
        """run生成对照库数据"""
        self.forwark()
        res = []
        logger.info('计算对照库')
        target_file = self.config['base'][NT_TARGET_STR]
        data = load_file(target_file, '\t', '#')
        if is_show:
            p = ProgressBar(len(data))
        for chr, start, end, exon, gene, *_ in data:
            start= int(start)
            end = int(end)
            exon_contrast = self.get_contrast(chr, start, end)
            if not exon_contrast:   #yzh
                continue
            res.append(exon_contrast)
            if is_show:
                p.read()
        contrast = BucketIndexIntervalList(res)
        contrast.make_index(10000)
        return contrast


def build_contrast(contrast_build_file, config_file, output=None, model='model2', is_show=False):
    """构建对照库
    contrast_build_file 对照库构建文件，wkcode-分析目录
    config_file 配置文件
    output 输出文件
    """
    model_obj = models_dict[model]
    model_choose = models_dict['model3']
    contrast_items = load_file(contrast_build_file, '\t', '#')
    if len(contrast_items) < 10:
        raise Exception('对照库样本少于10个')
    cb = ContrastBase(contrast_items, 
                      wkcode_cls=WkCodeBase, 
                      config_file=config_file, 
                      config=None, 
                      model=model_obj, model2=model_choose)
    contrast = cb.run(is_show)
    logger.info('数据格式化')
    target_file = cb.config['base'][NT_TARGET_STR]
    targets = load_file(target_file, '\t', '#')
    res = [(*item, *attrgetter('mean', 'cv', 'std', 'dps_mean')(obj), 0, 0, 0, 0) for item, obj in zip(targets, contrast)]
    if output:
        output_pic = '%s.%s'%(os.path.splitext(output)[0], 'pic')
        cb.ds.save(output_pic)
        output_pic2 = '%s_contrast.%s'%(os.path.splitext(output)[0], 'pic')
        IO_pickle.write_pickle(contrast, output_pic2)
        headers = ('chr', 'start', 'end', 'cds', 'gene', 
                   'ave', 'cv', 'sd', 'depth_ave',
                   'gain1_rate', 'gain2_rate', 'loss1_rate', 'loss2_rate')
        fmt_write_file(res, output, headers=headers)


import fire

if __name__ == "__main__":
    fire.Fire(build_contrast)




