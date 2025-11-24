#-*-coding:utf-8-*-
"""
NQB数据
只有父母没有先证者
合并父母相互矫正的gene文件
"""
import re
from mod_tools.mod_data_IO import IO_format
from intervaltools.core import IntervalBase, BucketIndexIntervalList
from exon_pipeline.utils.conv_tools import conv_filter


# reg = re.compile('[0-9A-Z][A-Z]{2}\d{3}')
# 用以兼容外来数据文库号，如“000GGG”
reg = re.compile('[A-Z]{3}[0-9]{3}')
reg1 = re.compile('[0-9]{3}[A-Z]{3}')


def load_tsv(file):
    """读取tsv文件返回headers和数据res"""
    headers = []
    with open(file, 'r', encoding='utf-8') as f:
        for i in f:
            if i.startswith('##'):
                continue
            if i.startswith('#'):
                headers = i[1:].strip().split()
                break
        res = IO_format.stdin_format_2DArray(f)
    return headers, res

def get_wkcode_index_from_headers(headers):
    """从gene文件表头获取wkcode索引"""
    r = {}
    for index, h in enumerate(headers):
        # print(">>", headers)
        m = reg.search(h)
        m1 = reg1.search(h)
        if m:
            r.setdefault(m.group(), []).append(index)
        elif m1:
            r.setdefault(m1.group(), []).append(index)
    return r

def check_swap(headers1, headers2):
    """检测headers是否可以交换如果可以则返回交换列index
    检测 
    1. 有且只有两个文库
    2. 两文件文库互补
    3. 两文件文库index互补
    检测不通过 raise AssertionError
    通过返回互补的indexs items
    """
    r_dict1 = get_wkcode_index_from_headers(headers1)
    r_dict2 = get_wkcode_index_from_headers(headers2)
    assert len(r_dict1) == 2
    assert len(r_dict2) == 2
    assert set(r_dict1) == set(r_dict2)
    key_list = list(r_dict1)
    assert r_dict1[key_list[0]] == r_dict2[key_list[1]]
    assert r_dict1[key_list[1]] == r_dict2[key_list[0]]
    return list(r_dict1.values())

def swap_line(indexs_items):
    """返回line的交换函数"""
    def func(line):
        line = list(line)
        for x,y in zip(*indexs_items):
            line[x], line[y] = line[y], line[x]
        return line
    return func


class IntervalLine(IntervalBase):
    def __init__(self, line):
        self.line = line
        super().__init__(*line[:3])

def merge_conv(conv1, conv2):
    """合并区间 直接在conv1上修改"""
    conv1._data.extend(conv2._data)
    conv1._data = sorted(conv1._data)
    conv1.make_index()



def merge_gene_file(gene_file1, gene_file2, no_filter=False):
    """
    NQB数据
    只有父母没有先证者
    合并父母相互矫正的gene文件
    以gene_file1为视角。
    gene_file2补充与gene_file1无交集的部分，并且交换文件内文库内容
    返回items: 表头 合并后数据列表
    """
    # 读取文件
    headers1, res1 = load_tsv(gene_file1)
    headers2, res2 = load_tsv(gene_file2)
    
    # 文件检测
    try:
        indexs_items = check_swap(headers1, headers2)
    except:
        raise Exception("文件检测未通过")
    swap = swap_line(indexs_items)
    
    # 预处理
    conv1 = BucketIndexIntervalList(res1, key=IntervalLine, default_base=IntervalLine)
    conv2 = BucketIndexIntervalList(res2, key=IntervalLine, default_base=IntervalLine)
    conv1.make_index()
    conv2.make_index()
    
    # gene_file2过滤以及交换位置
    if not no_filter:
        conv_filter(conv2, key=lambda x: not conv1.find_one(x, action='intersect'))
    conv2.add_attr('line', key=lambda x:swap(x.line))
    
    #合并
    merge_conv(conv1, conv2)
    res = [i.line for i in conv1]
    
    return headers1, res



    
