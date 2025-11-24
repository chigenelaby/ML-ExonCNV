#-*-coding:utf-8-*-
"""conv的相关函数"""

import os
import logging
import fire
from itertools import *
import numpy as np

def conv_filter(conv, key=None):
    """过滤相关信息"""
    if key is None:
        key = lambda x: x
    conv._data = [i for i in conv if key(i)]
    conv.make_index()

def get_conv_diff(conv):
    """获取conv的diff值
    只计算非bad_capture的diff
    """
    try:
        v = [ex.diff for ex in conv.exons if not ex.contrast.flag_bad_capture]
        if not v:
            raise
        res = round(np.mean(v), 4)
    except:
        res = round(np.mean([ex.diff for ex in conv.exons]), 4)
    return res




