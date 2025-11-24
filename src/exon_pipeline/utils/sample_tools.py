#-*-coding:utf-8-*-
"""样本工具"""

import configparser
from mod_tools.mod_data_IO import load_file

def get_info(path, key, type_cls=float):
    """获取统计文件的某个值"""
    data = load_file(path, '\t', '#')
    data_dict = {line[0]: line[1] for line in data if len(data)>1}
    tag = data_dict[key]
    try:
        tag = type_cls(tag)
    except:
        pass
    return tag

def get_info_reverse(path, key, type_cls=float):
    """获取统计文件的某个值 倒转的文件"""
    data = load_file(path, '\t', '#')
    data = list(zip(*data))
    data_dict = {line[0]: line[1] for line in data if len(data)>1}
    tag = data_dict[key]
    try:
        tag = type_cls(tag)
    except:
        pass
    return tag

def read_config(config_file):
    """读取config文件"""
    config = configparser.ConfigParser()
    config.read(config_file)
    return config

def fmt_id_string(chr, start, end, tag1, tag2):
    map1 = {'loss1': 'D',
            'loss2': 'D',
            'NA': 'wild',
            'wild': 'wild',
            'gain1': 'I',
            'gain2': 'I',
            'loss': 'D',
            'gain': 'I', }
    map2 = {'loss1': 'loss1',
            'loss2': 'loss2',
            'NA': 'wild',
            'wild': 'wild',
            'gain1': 'gain1',
            'gain2': 'gain2',
            'loss': 'loss',
            'gain': 'gain', }
    return '%s_%s_%s_%s_%s_C'%(chr, start, end, map1.get(tag1, 'wild'), map2.get(tag2, 'wild'))
