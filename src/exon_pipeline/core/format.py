#-*-coding:utf-8-*-
"""格式化"""
import os
import re
import logging
from operator import *

from intervaltools.core import IntervalBase, BucketIndexIntervalList
from exon_pipeline.utils.database_store import get_index

logger = logging.getLogger(__name__)

# 在bed文件中，提取exon的编号改为提取reg的编号
reg_cds1 = re.compile(r'CDS:(\d+)')
reg_cds2 = re.compile(r'REG:(\d+)')
reg_exon1 = re.compile(r'EXON:(\d+)')
reg_exon2 = re.compile(r'REG:(\d+)')


class ExonCDS(IntervalBase):
    """
    get_index('/share/chg2master/train/housy/CV_combin/NT01T_target.bed_mod', default_base=ExonCDS, key=itemgetter(0,1,2,3,4))
    """
    def __init__(self, chr, start, end, cds, gene):
        super().__init__(chr, start, end)
        self.reg_cds1 = reg_cds1
        self.reg_exon1 = reg_exon1
        self.reg_cds2 = reg_cds2
        self.reg_exon2 = reg_exon2
        self.cds_num = self.get_info(self.reg_cds1, self.reg_cds2, cds)
        self.exon_num = self.get_info(self.reg_exon1, self.reg_exon2,  cds)
        self.cds = cds
        self.gene = gene

    def get_info(self, reg1, reg2, string):
        """获取格式化数据"""
        try:
            res = reg1.search(string).group(1)
            res = int(res)
        except:
            try:
                res = reg2.search(string).group(1)
                res = int(res)
            except:
                res = None
        return res


class GeneInterval(IntervalBase):
    def __init__(self, exon_cds):
        chrs = {i.chr for i in exon_cds}
        assert len(chrs) == 1
        chr = exon_cds[0].chr
        genes = {i.gene for i in exon_cds}
        assert len(genes) == 1
        gene = exon_cds[0].gene
        exon_cds = sorted(exon_cds)
        start = min([i.start for i in exon_cds])
        end = max([i.end for i in exon_cds])
        super().__init__(chr, start, end)
        self.exons = exon_cds
        self.gene = gene
        self.exon_info = self.get_str_info('exon_num', 'EXON')
    
    def get_str_info(self, attr='exon_num', value='EXON'):
        """设置连接信息"""
        no_res = '%s:NA'%(value)
        counter = [getattr(exon, attr) for exon in self.exons]
        counter = [c for c in counter if c]
        if not counter:
            return no_res
        min_num = min(counter)
        max_num = max(counter)
        if min_num == max_num:
            return '%s:%s'%(value, min_num)
        return '%s:%s-%s'%(value, min_num, max_num)


def make_gene_intervals(gene_intervals_path, default_base=ExonCDS, key=itemgetter(0,1,2,3,4)):
    """构建gene区间文件"""
    exon_cds_list = get_index(gene_intervals_path, default_base=default_base, key=key)
    
    gene_dict = {}
    for i in exon_cds_list:
        gene_dict.setdefault((i.gene, i.chr), []).append(i)
    
    genes = [GeneInterval(exs) for info, exs in gene_dict.items()]
    gene_intervals = BucketIndexIntervalList(genes)
    gene_intervals.make_index()
    return gene_intervals



class GeneFormat:
    """
    >>> gene_intervals = make_gene_intervals(gene_intervals_path)
    >>> gf = GeneFormat(convs, gene_intervals)
    >>> gf.forward()
    >>> convs[3].gene_divides[0].exon_info
    'EXON:10-11'
    """
    def __init__(self, convs, gene_intervals):
        self.convs = convs
        self.gene_intervals = gene_intervals
    
    def intersect_exons(self, exons1, exons2):
        """返回exons1相同的外显子集合"""
        s1 = {i.interval() for i in exons1}
        s2 = {i.interval() for i in exons2}
        s = s1 & s2
        return [exon for exon in exons1 if exon.interval() in s]
    
    def set_one_conv(self, conv):
        """设置一个conv的gene划分"""
        gene_intervals = self.gene_intervals
        divides = []
        genes = gene_intervals.find_all(conv, action='intersect')
        for gene in genes:
            ist_exons = self.intersect_exons(gene.exons, conv.exons) ## by housy 20230105
            if not ist_exons:
                continue
            gene_divide = GeneInterval(ist_exons)
            setattr(gene_divide, 'is_all_gene', (gene_divide.exon_info==gene.exon_info))
            # print(gene_divide.source_exons)
            divides.append(gene_divide)
        return divides
    
    def forward(self):
        convs = self.convs
        convs.add_attr('gene_divides', key=self.set_one_conv)
        convs.check_attr('gene_divides', errors='strict')
        
