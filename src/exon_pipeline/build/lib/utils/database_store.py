#-*-coding:utf-8-*-
"""数据存储"""

import os
import re
import logging
from operator import itemgetter

from mod_tools.mod_data_IO import *
from intervaltools.core import IntervalBase, BucketIndexIntervalList




class DatabaseStore:
    """数据库结构"""
    def __init__(self, save_path=None, base=IntervalBase):
        self._store = []
        self.base = base
        self.initialized = False
        if save_path:
            self.load(save_path)
    
    def add_sample(self, intervals):
        """添加一个对照样本"""
        store = self._store
        store.append(intervals)
    
    def initialize(self,):
        """初始化"""
        if self.initialized:
            return
        self.initialized = True
    
    def get_data(self, chr, start, end):
        """获取数据"""
        self.initialize()
        tag_interval = self.base(chr, start, end, )
        res_intervals = []
        for intervals in self._store:
            match_intervals = intervals.find_all(tag_interval, action='be_included')
            res_intervals.append(match_intervals)
        return res_intervals
    
    def save(self, output):
        """存储信息"""
        self.initialize()
        IO_pickle.write_pickle(self._store, output)
    
    def load(self, save_path):
        """读取数据"""
        self._store = IO_pickle.load_pickle(save_path)
        self.initialized = True
    
    def __iter__(self):
        """iter"""
        return iter(self._store)


def get_index(path, default_base=IntervalBase, key=itemgetter(0,1,2,6)):
    """读入文件并构建桶索引"""
    data = load_file(path, '\t', '#')
    indexs = BucketIndexIntervalList(data, default_base=default_base, key=key)
    indexs.make_index(10000)
    return indexs


def get_filter_index(path, low_coverage_path, default_base=IntervalBase, key=itemgetter(0,1,2,6) ):
    """读入target.bed文件并对其进行过滤，再构建桶索引"""
    data = load_file(path, '\t', '#')
    filter_interval = get_index(low_coverage_path, default_base=IntervalBase, key=itemgetter(0,1,2)) 
    _filter_data = []
    for i in data:
        ex = IntervalBase(i[0], i[1], i[2])
        a = filter_interval.find_all(ex, action='intersect')
        if not a:
            _filter_data.append(i)
    indexs = BucketIndexIntervalList(_filter_data, default_base=default_base, key=key)
    indexs.make_index(10000)
    return indexs
