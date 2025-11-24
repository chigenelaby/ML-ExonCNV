#-*-coding:utf-8-*-
"""exon模块格式输出工具"""

import json
from operator import attrgetter

def fmt_id_string(chr, start, end, tag1, tag2):
    """CNVseq标签"""
    map1 = {'loss1': 'D',
            'loss2': 'D',
            'loss': 'D',
            'NA': 'wild',
            'wild': 'wild',
            'gain1': 'I',
            'gain2': 'I', 
            'gain': 'I', } 
    map2 = {'loss1': 'loss1',
            'loss2': 'loss2',
            'loss': 'loss', 
            'NA': 'wild',
            'wild': 'wild',
            'gain1': 'gain1',
            'gain2': 'gain2', 
            'gain': 'gain', }
    return '%s_%s_%s_%s_%s_C'%(chr, start, end, map1[tag1], map2[tag2])

def get_fmt_tag(conv):
    """格式化tag值"""
    tag = conv.tag
    if conv.is_male_sex_chr:
        tag = tag[:4]
    return tag

def get_fmt_tag2(tag, is_male_sex_chr=False):
    """格式化tag值"""
    if is_male_sex_chr:
        tag = tag[:4]
    return tag

def fmt_infos_json(conv, attrs=None):
    """格式化信息，生成json"""
    res = {}
    if attrs is None:
        pass
    else:
        for attr in attrs:
            try:
                info = getattr(conv, attr)
                res[attr] = info
            except:
                pass
    return json.dumps(res, ensure_ascii=False)

def update_infos_json(json_info, conv, attrs=None):
    """更新conv的json字段"""
    if attrs is None:
        return json_info
    try:
        res = json.loads(json_info, encoding='utf-8')
    except:
        res = json.loads(json_info)
    # res = json.loads(json_info, encoding='utf-8')
    for attr in attrs:
        info = getattr(conv, attr)
        res[attr] = info
    return json.dumps(res, ensure_ascii=False)

def update_infos_json_dict(json_info, **attrs_dict):
    """更新conv的json字段"""
    if attrs_dict is None:
        return json_info
    try:
        res = json.loads(json_info, encoding='utf-8')
    except:
        res = json.loads(json_info)
    for attr, info in attrs_dict.items():
        res[attr] = info
    return json.dumps(res, ensure_ascii=False)
