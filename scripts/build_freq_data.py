#-*-coding:utf-8-*-
"""构建频率库"""

import os
import logging
import fire
from itertools import *
from exon_pipeline.core.contrast import *
from exon_pipeline.core.controlled import *
from exon_pipeline.core.conv import *
from exon_pipeline.utils.database_store import DatabaseStore
from mod_tools import *
from mod_tools.mod_data_IO import *


def func_filter(conv, key=None):
    if key is None:
        key = lambda x: x
    conv._data = [i for i in conv if key(i)]
    conv.make_index()

filter_key = lambda x: (x.tag!='NA' and x.flag_diff.classify != False)


def build(ana_paths, output, 
          file_name='flag_info.pic', key=itemgetter(2), 
          filter_key=filter_key):
    data = load_file(ana_paths, '\t', '#')
    
    ds = DatabaseStore()
    p = ProgressBar(len(data))
    for i in data:
        pth = os.path.join(key(i), file_name)
        conv = load_pickle(pth)
        func_filter(conv, filter_key)
        if len(conv) < 10000:
            ds.add_sample(conv)
        p.read()
    ds.save(output)
