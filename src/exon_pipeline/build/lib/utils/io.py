#-*-coding:utf-8-*-
"""读写工具"""
import re
from collections import abc
import pandas as pd

def fmt_write_file(data, filename, key=None, headers=None, infos=None):
    """
    写数据(二维list)
    添加表头和信息头
    """
    end = '\n'
    sep = '\t'
    if key is None:
        key = lambda x: x
    with open(filename, 'w', encoding='utf-8') as f:
        # infos
        if infos is None:
            pass
        elif isinstance(infos, dict):
            for k, v in infos.items():
                print('##', k, v, sep=' ', file=f)
        elif isinstance(infos, abc.Iterable):
            for line in infos:
                print('##', *line, sep=' ', file=f)
        # headers
        if headers:
            print('#', end='', file=f)
            print(*headers, sep=sep, file=f)
        for line in data:
            line = key(line)
            print(*line, sep=sep, end=end, file=f)


def show(data, attrs=None, errors='strict', fill='NA'):
    """将数据按属性列表格式化
    data list(obj) obj数据
    attrs list(str) 属性列表
    errors 当属性不存在时的操作
        strict  引发 ValueError
        replace 替换为fill值
    fill 当errors=replace 时的默认值
    return : list(tuple())
    """
    if not attrs:
        return data
    res = []
    for obj in data:
        line = []
        for attr in attrs:
            try:
                if '.' in attr:
                    ele = obj
                    for attr_stage in attr.split('.'):
                        ele = getattr(ele, attr_stage)
                else:
                    ele = getattr(obj, attr)
            except:
                if errors == 'strict':
                    if attr == "reliability_value":
                        # print("ValueError", obj, attr)
                        ele = "未知"
                    else:
                        print("ValueError::", obj, attr)
                        raise ValueError
                elif errors == 'replace':
                    ele = fill
            line.append(ele)
        res.append(line)
    return res

def round_table(data, ndigits=None):
    """将二维列表每一个元素table进行round格式化"""
    def fmt(v):
        if not ndigits:
            return v
        if (isinstance(v, int) or isinstance(v, float)) and not isinstance(v, bool):
            v = round(v, ndigits)
        return v
    
    res = []
    for i in data:
        res.append((fmt(j) for j in i))
    return res

def read_table_to_df(path):
    """读取标准table文件，生成df
    跳过双#注释
    单#为表头
    """
    with open(path, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('##'):
                continue
            if line.startswith('#'):
                headers = line[1:].split()
                break
            else:
                raise Exception('格式错误无表头或者无效数据')
        else:
            raise Exception('格式错误无表头或者无效数据')
        df = pd.read_table(f, names=headers)
    return df

reg = re.compile('([1-9xyXY]\d*)(?=(\.\d+)*$)')
_chromosome_map = {'23': 'X', '24': 'Y', 'x': 'X', 'y': 'Y'}

def chr_conv(chr):
    """MC格式转化chr"""
    num = reg.search(chr).group(1)
    fmt_num = _chromosome_map.get(num, num)
    return 'chr%s'%fmt_num

