#!/user/bin/python
# -*- coding: utf-8 -*-

import sys
import os
import argparse
import collections
import glob
from time import strftime, gmtime, localtime
import gzip
import math, statistics
from io import StringIO

# =========================================
def print_log(msg):
    c_time = strftime("%Y-%m-%d %H:%M:%S", localtime())
    print(f">>> {c_time}\t{msg}")


def print_listlog(s, l):
    print(s)
    for a in l:
        print(a)


def write_tofile(f, buffer_o):
    with open(f, 'w') as outF:
        for o in buffer_o:
            outF.write(f"{o}\n")

# =========================================
def fun_avg1(l):
    if len(l) == 0:
        return 0
    t = sum([x.diff for x in l])
    a = t / len(l)
    return a


def fun_avg2(l, mut_r1=1/2.3, mut_r2=1.3/2.3):
    if len(l) == 0 or None in [mut_r1, mut_r2]:
        return ""
    
    q1, q2 = expd_q(mut_r1, [0.05, 0.45], [0.15, 0.075])
    q3, q4 = expd_q(mut_r2, [0.55, 0.95], [0.075, 0.15])
    # print("mut-lmt-q", len(l), q1, q2, q3, q4)
    
    l1 = [x.get_mutratio() for x in l]
    l1a, l1b, l1c, l1d, l1e = [], [], [], [], []
    for x in l1:
        if x >= 0 and x <= q1:
            l1a.append(x)
        elif x > q1 and x <= q2:
            l1b.append(x)
        elif x > q2 and x <= q3:
            l1c.append(x)
        elif x > q3 and x <= q4:
            l1d.append(x)
        else:
            l1e.append(x)

    t1, t1a, t1b, t1c, t1d, t1e = len(l1), len(l1a), len(l1b), len(l1c), len(l1d), len(l1e)
    r1a, r1b, r1c, r1d, r1e = t1a/t1, t1b/t1, t1c/t1, t1d/t1, t1e/t1
    if t1b + t1d == 0:
        r_q_half = None
    elif t1c == 0:
        r_q_half = 1000
    else:
        r_q_half = (t1b + t1d) / t1c
    
    s = "{}|{}-{}-{}-{}-{}|{:.2f}-{:.2f}-{:.2f}-{:.2f}-{:.2f}|{}".format(t1, t1a, t1b, t1c, t1d, t1e, r1a, r1b, r1c, r1d, r1e, r_q_half)
    return s


def stat_chr(d, prop, fun_avg):
    for i in range(1, 25):
        l = d[i]
        l1 = list(map(lambda x: getattr(x, prop), l))
        if len(l1) == 0:
            print(f"chr{i} empty {prop} list")
            continue
        a1 = min(l1)
        a2 = max(l1)
        reg = 20_000_000
        step = 5_000_000
        
        a = a1
        while a <= a2 + step:
            b = a + reg
            l2 = list(filter(lambda x: getattr(x, prop) >= a and getattr(x, prop) < b, l))
            print(f"chr{i}", a, b, len(l2), fun_avg(l2))
            a = a + step

# =========================================
def chr_tonum(tmpchr):
    if isinstance(tmpchr, int):
        return tmpchr

    tmpchr = tmpchr.strip()

    if tmpchr.startswith("NC"):
        return int(tmpchr[7:9])
    elif tmpchr.startswith("chr"):
        if tmpchr not in {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}:
            return 0
        if tmpchr == "chrX":
            return 23
        elif tmpchr == "chrY":
            return 24
        else:
            return int(tmpchr[3:])
    elif tmpchr == "X":
        return 23
    elif tmpchr == "Y":
        return 24
    else:
        return int(tmpchr)


def chrnum_tochr(i, fmt='long'):
    if i == 23:
        i = 'X'
    elif i == 24:
        i = 'Y'
    else:
        pass
    
    if fmt == 'long':
        i = f"chr{i}"

    return i


def get_chr_startend(data_type='all'):
    info_str1 = """
NC_000001.10    1  249250621
NC_000002.11    1  243199373
NC_000003.11    1  198022430
NC_000004.11    1  191154276
NC_000005.9     1  180915260
NC_000006.11    1  171115067
NC_000007.13    1  159138663
NC_000008.10    1  146364022
NC_000009.11    1  141213431
NC_000010.10    1  135534747
NC_000011.9     1  135006516
NC_000012.11    1  133851895
NC_000013.10    1  115169878
NC_000014.8     1  107349540
NC_000015.9     1  102531392
NC_000016.9     1   90354753
NC_000017.10    1   81195210
NC_000018.9     1   78077248
NC_000019.9     1   59128983
NC_000020.10    1   63025520
NC_000021.8     1   48129895
NC_000022.10    1   51304566
NC_000023.10    1  155270560
NC_000024.9     1   59373566
"""
    info_str2 = """
NC_000001.10    721369  249228211
NC_000002.11    10001   243059133
NC_000003.11    60001   197897252
NC_000004.11    10001   190982280
NC_000005.9     30001   180718728
NC_000006.11    140001  170919973
NC_000007.13    50001   159128663
NC_000008.10    150001  146302589
NC_000009.11    210001  141046998
NC_000010.10    100001  135297528
NC_000011.9     180001  134946516
NC_000012.11    145740  133836993
NC_000013.10    19300001        115109878
NC_000014.8     19160001        107289540
NC_000015.9     20160001        102284474
NC_000016.9     140001  90189384
NC_000017.10    1       81029050
NC_000018.9     50001   78001821
NC_000019.9     240001  59101783
NC_000020.10    160001  62913370
NC_000021.8     10777897        48098051
NC_000022.10    16090001        51144778
NC_000023.10    321385  155157100
NC_000024.9     271385  59257657
"""
    info_str3 = """
NC_000001.10    721369  249228211
NC_000002.11    10001   243059133
NC_000003.11    60001   197897252
NC_000004.11    10001   190982280
NC_000005.9     30001   180718728
NC_000006.11    140001  170919973
NC_000007.13    50001   159128663
NC_000008.10    150001  146302589
NC_000009.11    210001  141046998
NC_000010.10    100001  135297528
NC_000011.9     180001  134946516
NC_000012.11    145740  133836993
NC_000013.10    19300001        115109878
NC_000014.8     19160001        107289540
NC_000015.9     20160001        102284474
NC_000016.9     140001  90189384
NC_000017.10    1       81029050
NC_000018.9     50001   78001821
NC_000019.9     240001  59101783
NC_000020.10    160001  62913370
NC_000021.8     10777897        48098051
NC_000022.10    16090001        51144778
NC_000023.10    2699520  154931044
NC_000024.9     2699520  28819361
"""
# 去除伪常区、Y染色体末尾N区
# NC_000023.10    321385  155157100
# NC_000024.9     271385  59257657

    list_chr_start = []
    list_chr_end = []

    tmp = 0

    if data_type == 'all':
        info_str = info_str1
    elif data_type == 'nn':
        info_str = info_str2
    elif data_type == 'cut':
        info_str = info_str3
    else:
        info_str = info_str1

    for tmp_line in StringIO(info_str):
        if tmp_line.strip() != "":
            list_lineparts = tmp_line.strip().split()
            # print(len(list_lineparts))
            # tmp += 1

            list_chr_start.append(int(list_lineparts[1]))
            list_chr_end.append(int(list_lineparts[2]))

    # print(tmp)
    # print(list_chr_start)
    # print(list_chr_end)

    return list_chr_start, list_chr_end

# =========================================
def which_mostype(a, chrnum, sex):
    t = None
    if sex == 'XY' and chrnum in [23, 24]:
        if a > 0.1 and a < 0.95:
            t = 'lossmos'
        if a > 1.1 and a < 1.9:
            t = 'gainmos'
    else:
        if a > 0.55 and a < 0.95:
            t = 'lossmos'
        if a > 1.05 and a < 1.4:
            t = 'gainmos'
    return t


def whether_refu_mostype(mostype, list_otherreg):
    def fun_refuse(t, a, total_exon):
        # if total_exon > 150:
        #     thld_refu_gain = 1
        # else:
        #     thld_refu_gain = 0.98
        
        if t == 'lossmos':
            b = a > 1
        elif t == 'gainmos':
            b = a < 1
        else:
            b = None
        return b
    l_refuse = [len(x.list_exon) > 20 and fun_refuse(mostype, x.diff, len(x.list_exon)) for x in list_otherreg]
    
    return any(l_refuse)


def get_qs_diff(list_exon):
    l = [x.diff for x in list_exon]
    l.sort()
    if len(l) >= 3:
        i1, i2, i3 = int(0.25*len(l)), int(0.5*len(l)), int(0.75*len(l))
        l1 = [l[i1], l[i2], l[i3]]
    elif len(l) > 0:
        avg = sum(l)/len(l)
        l1 = [avg, avg, avg]
    else:
        l1 = [None, None, None]
    return l1

# =========================================
class PART():
    chrnum: int
    start: int
    end: int
    
    def __init__(self, chrnum, start, end):
        self.chrnum = chrnum
        self.start = start
        self.end = end
    
    def get_len(self):
        return self.end - self.start + 1

    def whether_ovlap(self, a):
        b1 = self.chrnum == a.chrnum
        b2 = self.start < a.end
        b3 = a.start < self.end
        return b1 and b2 and b3

    def get_ovlap_len(self, a):
        if self.whether_ovlap(a):
            o1 = max(self.start, a.start)
            o2 = min(self.end, a.end)
            l = o2 - o1 + 1
        else:
            l = 0
        return l
    
    def get_unovlap_part(self, a):
        if not self.whether_ovlap(a):
            l = [self, a]
        else:
            start = min(self.start, a.start)
            o1 = max(self.start, a.start)
            o2 = min(self.end, a.end)
            end = max(self.end, a.end)
            if start == o1:
                p1 = None
            else:
                p1 = PART(self.chrnum, start, o1-1)
            if o2 == end:
                p2 = None
            else:
                p2 = PART(self.chrnum, o2+1, end)
            
            l = [p1, p2]
        return l

    def get_unovlap_len(self, a):
        if self.chrnum != a.chrnum:
            tmplen = None
        else:
            tmplen = 0
            l = self.get_unovlap_part(a)
            for p in l:
                if p is not None:
                    tmplen += p.get_len()
        
        return tmplen


class EXON(PART):
    diff: float
    depth: float
    
    def __repr__(self):
        return "{}-{}-{}-{:.2f}-{:.1f}".format(self.chrnum, self.start, self.end, self.diff, self.depth)
    
    def __init__(self, chrnum, start, end, diff, depth):
        super().__init__(chrnum, start, end)
        self.diff = diff
        self.depth = depth


class REG(EXON):
    cnvtype: str
    prt_cell: float
    # 没有exon区域所占比例
    prt_uncover: float
    # diff特别高的exon比例，用以表征数据质量差的区段
    prt_high: float
    total_high: int
    list_exon: list
    qs_diff: list
    evi_exon: bool
    evi_mut: bool
    info: dict

    def __init__(self, chrnum, start, end, diff, depth, cnvtype, list_exon):
        super().__init__(chrnum, start, end, diff, depth)
        self.cnvtype = cnvtype
        self.prt_cell = 0
        self.prt_uncover = 1
        self.prt_high = 1
        self.total_high = None
        self.list_exon = list_exon
        self.qs_diff = get_qs_diff(self.list_exon)
        self.evi_exon = None
        self.evi_mut = None
        self.info = {}
    
    def __repr__(self):
        s = "|{}-prtcell:{:.2f}-exontotal:{}-prtuncover:{:.2f}-prthigh:{:.4f}-qsdiff:{}-{}-{}-{}".format(self.cnvtype, self.prt_cell, self.get_exontotal(), self.prt_uncover, self.prt_high, self.qs_diff, self.evi_exon, self.evi_mut, self.info)
        return super().__repr__() + s
    
    def short_rep(self):
        return super().__repr__()

    def merge(self, tmpreg):
        self.start = min(self.start, tmpreg.start)
        self.end = max(self.end, tmpreg.end)
        
        self.list_exon.extend(tmpreg.list_exon)
        self.list_exon.sort(key=lambda x: x.start)
        # 合并的阈值比较宽，所以不用平均数
        list_diff = [x.diff for x in self.list_exon]
        list_diff.sort()
        if self.chrnum == 24:
            factor = 0.5
        elif self.diff > 1:
            if self.get_len() > 50_000_000 or max(self.diff, tmpreg.diff) > 1.25:
                factor = 0.5
            elif self.get_len() > 25_000_000 or max(self.diff, tmpreg.diff) > 1.15:
                factor = 0.6
            else:
                factor = 0.75
        else:
            factor = 0.5
        if len(list_diff) > 10:
            self.diff = list_diff[int(factor*len(self.list_exon))]
        else:
            self.diff = fun_avg1(self.list_exon)
            
        self.qs_diff = get_qs_diff(self.list_exon)

    def get_exontotal(self):
        return len(self.list_exon)

    def whether_ovlap(self, a):
        return super().whether_ovlap(a) and self.cnvtype == a.cnvtype

    def get_partlen(self, a, list_otherreg):
        l1 = [self, a]
        l1.sort(key=lambda x: x.start)
        mid_len = l1[1].start - l1[0].end - 1
        reg_t_len = l1[1].end - l1[0].start + 1
        
        list_otherreg.sort(key=lambda x: x.start)
        
        a = None
        list_gap = []
        list_gap1 = []
        gap_len = 0
        for x in list_otherreg:
            # print("info-get_partlen1", x.short_rep(), len(x.list_exon))
            # 不能只用起止位点，gap区也可能被化成bin
            if a is not None and x.start - a.end > 500_000:
                list_gap.append([a.end+1, x.start-1])
            if len(x.list_exon) < 3:
                list_gap.append([x.start, x.end])
            a = x
                            
        for a in list_gap:
            bool_add = False
            if len(list_gap1) > 0:
                if a[0] <= list_gap1[-1][1] and list_gap1[-1][0] <= a[1]:
                    list_gap1[-1] = [min(a[0], list_gap1[-1][0]), max(a[1], list_gap1[-1][1])]
                    bool_add = True
            if not bool_add:
                list_gap1.append(a)
        # print("info-get_partlen", self.short_rep(), self.short_rep(), len(list_otherreg), list_gap1)
        
        for a in list_gap1:
            gap_len += a[1] - a[0] + 1
            
        return gap_len, mid_len, reg_t_len


    def whether_supp(self, a, list_otherreg):
        b = False
        if self.chrnum != a.chrnum:
            return b
        
        if self.whether_ovlap(a):
            b = True
        else:            
            gap_len, mid_len, reg_t_len = self.get_partlen(a, list_otherreg)
            # if self.chrnum == 9:
            #     print("info-whether_supp", self.short_rep(), a.short_rep(), gap_len, mid_len, reg_t_len, (mid_len - gap_len)/reg_t_len, whether_refu_mostype(self.cnvtype, list_otherreg), len(list_otherreg))
            
            if (mid_len - gap_len)/reg_t_len < 0.05 and not whether_refu_mostype(self.cnvtype, list_otherreg):
                b = True
        return b


class MUT():
    chrnum: int
    pos: int
    dp_ref: int
    dp_mut: int
    
    def __repr__(self):
        return f"{self.chrnum}-{self.pos}-{self.dp_ref}-{self.dp_mut}"
    
    def __init__(self, chrnum, pos, dp_ref, dp_mut):
        self.chrnum = chrnum
        self.pos = pos
        self.dp_ref = dp_ref
        self.dp_mut = dp_mut
        
    def get_dp(self):
        return self.dp_ref + self.dp_mut

    def get_mutratio(self):
        t = self.get_dp()
        if t == 0:
            r = 0
        else:
            r = self.dp_mut / t
        return r

# =========================================
def read_exondetail(f_exondetail):
    d = collections.defaultdict(list)
    with open(f_exondetail) as inF:
        for line in inF:
            if line.startswith("#"):
                continue
    
            lparts = line.rstrip("\n").split("\t")
            chrnum = chr_tonum(lparts[0])
            if chrnum == 0:
                continue
            start = int(lparts[1])
            end = int(lparts[2])
            diff = float(lparts[5])
            depth = float(lparts[6])
            d[chrnum].append(EXON(chrnum, start, end, diff, depth))
    # print("dict exondetail", d)
    return d


def whichsex(dictlib):
    # list_reg = list(filter(lambda x: x.start>2649520 and x.end<59034050, dictlib[24]))
    # 去除伪常区和AZF区
    list_reg = list(filter(lambda x: x.start>2_649_520 and x.end<17_000_000, dictlib[24]))
    # 如果没有Y染色体，则相邻位点的值会大幅度振荡
    a = 0
    list_abnorm = []
    for reg in list_reg:
        dp = reg.depth
        b = reg.diff
        bool_depth_abnorm = dp is None or math.isinf(dp) or math.isnan(dp) or dp<20
        bool_diff_abnorm = b is None or math.isinf(b) or math.isnan(b) or b<0.1
        bool_fluc = abs(b-a)>0.25 and min(a, b) < 0.2
        if bool_depth_abnorm or bool_diff_abnorm or bool_fluc:
            list_abnorm.append(reg)
        a = b
        
    total = len(list_reg)
    total_ab = len(list_abnorm)
    # print(list_copydiff)
    # print(list_copydiff_abnorm)
    if total == 0:
        ratio = 1
    else:
        ratio = total_ab / total
    if ratio < 0.4:
        sex = 'XY'
    else:
        sex = 'XX'
    print(f"||| sex {sex}::bintotalchrY {total}::bintotalchrY abnormal {total_ab}")
    return sex

# =========================================
def get_dpinfo(fmt, info):
    dp_ref = 0
    dp_mut = 0
    for fmt1, info1 in zip(fmt.split(":"), info.split(":")):
        if fmt1 == 'AD':
            l = info1.split(',')
            # print(fmt, info, l)
            if l[0] == '.':
                return dp_ref, dp_mut
            elif len(l) > 1:
                dp_ref = int(l[0])
                dp_mut = int(l[1])
            else:
                dp_ref = int(l[0])
    return dp_ref, dp_mut            


def read_gvcfmut(gvcf_mut):
    d = collections.defaultdict(list)
    with gzip.open(gvcf_mut, 'rb') as inF:
        for line in inF:
            line = line.decode(encoding='utf8')
                
            if line.startswith("#"):
                continue
    
            lparts = line.rstrip("\n").split("\t")
            chrnum = chr_tonum(lparts[0])
            if chrnum == 0:
                continue
            pos = int(lparts[1])
            dp_ref, dp_mut = get_dpinfo(lparts[8], lparts[9])
            d[chrnum].append(MUT(chrnum, pos, dp_ref, dp_mut))
    # print("dict gvcfmut", d)
    return d

# =========================================
def get_fun_cnvratio(chrnum, cnvtype, diff, sex, qs_diff=None, factor1=1, factor2=2):
    if sex == "XY" and chrnum in [23, 24]:
        # if cnvtype == 'gainmos':
        #     l = [1.1, 1.9]
        # elif cnvtype == 'lossmos':
        #     l = [0.1, 0.95]
        # else:
        #     l = None
        rng = 0.2
        if diff > 1.45:
            rng = 0.4
    else:
        # if cnvtype == 'gainmos':
        #     l = [1.05, 1.4]
        # elif cnvtype == 'lossmos':
        #     l = [0.6, 0.95]
        # else:
        #     l = None
        rng = 0.1
    
    if "gain" in cnvtype:
        if qs_diff is None or chrnum == 24:
            l = [max(diff-rng*factor1, 1.05), diff+rng*factor2]
        else:
            l = [max(qs_diff[0], 1.01), diff+rng*factor2]
    elif "loss" in cnvtype:
        l = [diff-rng*factor2, min(diff+rng*factor1, 0.95)]
    else:
        l = None
    
    if l is None:
        f = None
    else:
        f = lambda x: x>l[0] and x<l[1]
        
    return f


def whether_merge(a, b, list_otherreg, sex):   
    bool_merge = False
    bool_break = False
    # print("info-whether-in", reg.chrnum, reg.cnvtype, reg.diff)
    whether_in = get_fun_cnvratio(a.chrnum, a.cnvtype, a.diff, sex, factor1=2, factor2=2)
    
    if a.cnvtype != b.cnvtype:
        bool_break = True        
    elif whether_in is not None and whether_in(b.diff):
        list_otherreg1 = list(filter(lambda x: x.start>=min(a.start, b.start) and x.end<=max(a.end, b.end), list_otherreg))
        # print("info-whether_merge", len(list_otherreg1))
        if b.whether_supp(a, list_otherreg1):
            bool_merge = True
        else:
            bool_break = True
    else:
        pass

    return bool_merge, bool_break


def grow_mos(l_reg, reg1, list_otherreg, sex):
    l = []
    a = reg1
    while len(l_reg) > 0:
        b = l_reg.pop(-1)
        bool_merge, bool_break = whether_merge(a, b, list_otherreg, sex)
        
        if bool_merge:
            a.merge(b)
        elif bool_break:
            l = l_reg + [b] + l
            l_reg = []
        else:
            l = [b] + l
    
    l.append(a)
    return l


def cut_ovlap(l_reg):
    l_reg.sort(key=lambda x: x.get_len())
    l = []
    bool_rm = False
    while len(l_reg) > 0:
        a = l_reg.pop(-1)
        for b in l:
            if a.whether_ovlap(b):
                ovlap_point = [max(a.start, b.start), min(a.end, b.end)]
                if a.start == ovlap_point[0] and a.end == ovlap_point[1]:
                    bool_rm = True
                    break
                elif a.start == ovlap_point[0]:
                    a.start = ovlap_point[1] + 1
                elif a.end == ovlap_point[1]:
                    a.end = ovlap_point[0] - 1
                else:
                    pass
        if not bool_rm:
            l.append(a)
    
    return l


def filter_merge(dict_chrnum_exonlist, sex, prop='start'):
    l = []
    for i in range(1, 25):
        l_reg = []
        list_otherreg = []
        whether_in = None
        
        list_exon = dict_chrnum_exonlist[i]
        list_pos = list(map(lambda x: getattr(x, prop), list_exon))
        if len(list_pos) == 0:
            print(f"||| chr{i} empty list")
            continue

        # chrY去除伪常区和AZFc区
        # reg_end 15M，不是划分区段的终点，还要加step及bin的冗余
        if i == 24:
            reg_start = 2_649_520
            reg_end = 15_000_000
        else:
            reg_start = min(list_pos)
            reg_end = max(list_pos)
        len_bin = 5_000_000
        step = 1_000_000
        
        a = reg_start
        while a <= reg_end + step:
            b = a + len_bin
            list_exon1 = list(filter(lambda x: getattr(x, prop) >= a and getattr(x, prop) < b, list_exon))
            total_exon1, avg_diff1 = len(list_exon1), fun_avg1(list_exon1)
            mostype = which_mostype(avg_diff1, i, sex)
            # if i == 11 and a>=0 and a<=400_000_000:
            #     print("step-run", i, a, b, total_exon1, avg_diff1, mostype)
            
            reg1 = REG(i, a, b, avg_diff1, float('nan'), mostype, list_exon1)
            if total_exon1 > 20:
                if mostype in ['lossmos', 'gainmos']:
                    l_reg = grow_mos(l_reg, reg1, list_otherreg, sex)
                else:
                    list_otherreg.append(reg1)
            else:
                list_otherreg.append(reg1)
            
            # if reg is not None:
            #     print("filter_merge", i, a, b, total_exon1, avg_diff1, mostype, reg)
        
            a = a + step
        
        l.extend(cut_ovlap(l_reg))
    return l

# =========================================
def stat_1st(list_mos):
    total_reg1st = len(list_mos)
    l_gain = list(filter(lambda x: x.cnvtype=="gainmos", list_mos))
    l_loss = list(filter(lambda x: x.cnvtype=="lossmos", list_mos))
    if len(l_loss) == 0 and len(l_gain) == 0:
        ratio_gainloss = 1
    elif len(l_loss) == 0:
        ratio_gainloss = 1000
    else:
        ratio_gainloss = len(l_gain) / len(l_loss)
    return total_reg1st, round(ratio_gainloss, 1)
# =========================================

def refresh_pos(l, reg_part):
    if reg_part == 'start':
        pos = l[0].start
    else:
        pos = l[-1].end
   
    return pos


def cal_prtcell(mos, sex):
    p = abs(mos.diff - 1)
    if sex == "XY" and mos.chrnum in [23, 24]:
        pass
    else:
        p = p * 2
    # print("prtcell", mos, p)
    return p


def cut_step(l, reg_part):
    if reg_part == "start":
        a = 0
    else:
        a = -1    
    l.pop(a)
    return l


def get_cutpos(list_exon):
    i = int(len(list_exon) / 2)
    tmplen_max = 0
    bool_gapcut = False
    for i1 in range(1, len(list_exon)):
        tmplen = list_exon[i1].start - list_exon[i1-1].end
        if tmplen > 1_000_000 and tmplen > tmplen_max:
            i = i1
            tmplen_max = tmplen
            bool_gapcut = True
    return i, bool_gapcut


def cut_move(list_exon, mos, sex, pos, reg_part):
    whether_in1 = get_fun_cnvratio(mos.chrnum, mos.cnvtype, mos.diff, sex, qs_diff=mos.qs_diff, factor1=1, factor2=2)
    whether_in2 = get_fun_cnvratio(mos.chrnum, mos.cnvtype, mos.diff, sex, qs_diff=mos.qs_diff, factor1=1.5, factor2=2)
    if whether_in1 is None:
        return pos
    
    while(len(list_exon) > 1):
        i, bool_gapcut = get_cutpos(list_exon)
        l1 = list_exon[0:i]
        l2 = list_exon[i:]
        whether_in = whether_in1 if bool_gapcut else whether_in2

        if reg_part == "start":
            l_t = l2
            l_other = l1
        else:
            l_t = l1
            l_other = l2

        if bool_gapcut and len(l_other) < 20:
            list_exon = l_t
            continue

        avgdiff_t = sum([x.diff for x in l_t]) / len(l_t)
        avgdiff_other = sum([x.diff for x in l_other]) / len(l_other)
        
        # if mos.chrnum == 11:
        #     cutpos = l1[-1].end
        #     print("cut-move", f"{reg_part}-{mos.chrnum}-diff:{mos.diff}-pos:{pos}-cutpos:{cutpos}-{bool_gapcut}-{avgdiff_t}-{avgdiff_other}", "reg", list_exon[0].start, list_exon[-1].end, len(list_exon), len(l_t), len(l_other))
        
        if whether_in(avgdiff_t):
            if whether_in(avgdiff_other):
                list_exon = l_other
                pos = refresh_pos(list_exon, reg_part)
            else:
                list_exon = cut_step(list_exon, reg_part)
                pos = refresh_pos(list_exon, reg_part) 
        else:
            list_exon = l_t
            pos = refresh_pos(list_exon, reg_part) 

    return pos


def get_bkp(mos, dict_chrnum_exonlist, sex, reg_part='start'):
    pos = getattr(mos, reg_part)
    if mos.get_len() < 10_000_000:
        # print("====get_bkp_len_filter", mos)
        return pos
    
    posmid = 0.5 * (mos.start + mos.end)

    list_exon = list(filter(lambda x: x.start>min(pos-5_000_000, posmid) and x.end<max(pos+5_000_000, posmid), dict_chrnum_exonlist[mos.chrnum]))
    list_exon = sorted(list_exon, key=lambda x: x.start)
    if len(list_exon) < 20:
        # print("====get_bkp_totalexon_filter", mos, len(list_exon))
        return pos
    
    # print("get_bkp_info", pos-5_000_000, pos+5_000_000, reg_part)
    # print(mos, list_exon[0], list_exon[1], list_exon[-2], list_exon[-1])
    
    pos = cut_move(list_exon, mos, sex, pos, reg_part)
    return pos


def get_chr_ratiothld(prt_cell):
    if prt_cell > 0.4:
        thld = 0.9
    elif prt_cell > 0.3:
        thld = 0.85
    else:
        thld = 0.8
    return thld


def whether_whole_chr(mos, endpoint_thld=500_000):
    tmpchr_num, start, end, prt_cell = mos.chrnum, mos.start, mos.end, mos.prt_cell
    # 去除伪常区和AZF区
    # if tmpchr_num == 24 and start < 5_000_000 and end > 26_000_000:
    if tmpchr_num == 24 and start < 5_500_000 and end > 16_000_000:
        return True

    list_chr_start, list_chr_end = get_chr_startend(data_type='cut')
    chrstart, chrend = list_chr_start[tmpchr_num - 1], list_chr_end[tmpchr_num - 1]

    chr_startpos_thld = chrstart + endpoint_thld
    chr_endpos_thld = chrend - endpoint_thld

    tmplen = end - start + 1
    tmpchrlen = chrend - chrstart + 1
    ratio_cnv = tmplen / tmpchrlen
    # print("?wholechr", tmplen, tmpchrlen, ratio_cnv)

    ratio_thld = get_chr_ratiothld(prt_cell)

    tmpbool = False
    if (start < chr_startpos_thld and end > chr_endpos_thld) or ratio_cnv > ratio_thld:
        tmpbool = True
    else:
        pass

    return tmpbool


def calib_bkp(list_mos, dict_chrnum_exonlist, sex):
    for mos in list_mos:
        mos.start = get_bkp(mos, dict_chrnum_exonlist, sex, reg_part='start')
        mos.end = get_bkp(mos, dict_chrnum_exonlist, sex, reg_part='end')

        mos.list_exon = list(filter(lambda x: x.start>=mos.start and x.end<=mos.end, dict_chrnum_exonlist[mos.chrnum]))
        mos.diff = fun_avg1(mos.list_exon)
        mos.prt_cell = cal_prtcell(mos, sex)

        list_chr_start, list_chr_end = get_chr_startend(data_type='all')
        chrstart, chrend = list_chr_start[mos.chrnum - 1], list_chr_end[mos.chrnum - 1]
        # print(">>>calib_bkb1", mos)
        
        if whether_whole_chr(mos, endpoint_thld=5_000_000):
            mos.start = chrstart
            mos.end = chrend
        else:
            if mos.start < 500_000:
                mos.start = 1
            if mos.end > chrend - 500_000:
                mos.end = chrend
        
        # print(">>>calib_bkb2", mos)
    return list_mos

# =========================================
# 统计同一样本中，和该mos有相似分布的mos总数
def stat_samedist(list_mos):
    for mos in list_mos:
        # if mos.get_len() < 50_000_000:
        #     factor_low, factor_mid, factor_high = 0.3, 0.75, 0.95
        # else:
        #     factor_low, factor_mid, factor_high = 0.1, 0.5, 0.7
        factor_low, factor_mid, factor_high = 0.1, 0.5, 0.7
            
        if mos.qs_diff[1] < 1:
            thld_a = 0.5*(mos.qs_diff[1]+mos.qs_diff[2])
            tmp_fun1 = lambda x: x.qs_diff[1]<thld_a and x.get_len()>factor_low*mos.get_len()
            tmp_fun2 = lambda x: x.qs_diff[1]<thld_a and x.get_len()>factor_mid*mos.get_len()
            tmp_fun3 = lambda x: x.qs_diff[1]<thld_a and x.get_len()>factor_high*mos.get_len()
        else:
            thld_a = 0.5*(mos.qs_diff[0]+mos.qs_diff[1])
            tmp_fun1 = lambda x: x.qs_diff[1]>thld_a and x.get_len()>factor_low*mos.get_len()
            tmp_fun2 = lambda x: x.qs_diff[1]>thld_a and x.get_len()>factor_mid*mos.get_len()
            tmp_fun3 = lambda x: x.qs_diff[1]>thld_a and x.get_len()>factor_high*mos.get_len()

        tmp_l1 = list(filter(tmp_fun1, list_mos))
        tmp_l2 = list(filter(tmp_fun2, list_mos))
        tmp_l3 = list(filter(tmp_fun3, list_mos))
        
        # len_samedist = sum([x.get_len() for x in tmp_l1]) - mos.get_len()
        # ratio_samedist_len = len_samedist / mos.get_len()
        
        info = {'total_samedist_short':len(tmp_l1)-1, 'total_samedist_mid':len(tmp_l2)-1, 'total_samedist_long':len(tmp_l3)-1}
        mos.info.update(info)
    return list_mos

#============================
def check_exondist(list_mos, dict_chrnum_exonlist, sex):
    for mos in list_mos:        
        whether_in = get_fun_cnvratio(mos.chrnum, mos.cnvtype, mos.diff, sex)
        if whether_in is None:
            continue
        
        list_exon_f = list(filter(lambda x: whether_in(x.diff), mos.list_exon))
        total_exon_f, total_exon = len(list_exon_f), len(mos.list_exon)
        # print("mos-info", mos)
        if total_exon == 0:
            prt_exon_f = 0
        else:
            prt_exon_f = total_exon_f / total_exon
        info = {'total_exon_f':total_exon_f, 'total_exon':total_exon, 'prt_exon_f':prt_exon_f}
        
        if mos.chrnum == 24 and total_exon_f > 30:
            mos.evi_exon = True
        # elif (total_exon_f > 10 and prt_exon_f > 0.5) or (total_exon_f > 1000 and prt_exon_f > 0.4):
        elif prt_exon_f >= 0.5:
            mos.evi_exon = True
        elif total_exon > 20:
            mos.evi_exon = False
        else:
            pass
        mos.info.update(info)
    return list_mos

#============================
# 嵌合体理论上，点突变的突变率
# 只考虑二倍体情况
def model_mutratio(mos):
    r1, r2 = None, None
    if mos.prt_cell > 1:
        pass
    elif mos.cnvtype == "lossmos":
        r1 = (1-mos.prt_cell)/(2-mos.prt_cell)
        r2 = 1/(2-mos.prt_cell)
    elif mos.cnvtype == "gainmos":
        r1 = 1/(2+mos.prt_cell)
        r2 = (1+mos.prt_cell)/(2+mos.prt_cell)
    else:
        pass
    return r1, r2


def expd_q(mut_r, list_lmt, list_step):
    qa = max([mut_r-list_step[0], list_lmt[0]])
    qb = min([mut_r+list_step[1], list_lmt[1]])
    return qa, qb


def stat_mut(l, mut_r1, mut_r2):
    if len(l) == 0 or None in [mut_r1, mut_r2]:
        return 0, 0, 0, 0
    
    q1, q2 = expd_q(mut_r1, [0.05, 0.45], [0.25, 0.075])
    q3, q4 = expd_q(mut_r2, [0.55, 0.95], [0.075, 0.25])
    # print("mut-lmt-q", len(l), q1, q2, q3, q4)
    
    l1 = [x.get_mutratio() for x in l]
    l1a, l1b, l1c, l1d, l1e = [], [], [], [], []
    for x in l1:
        if x >= 0 and x <= q1:
            l1a.append(x)
        elif x > q1 and x <= q2:
            l1b.append(x)
        elif x > q2 and x <= q3:
            l1c.append(x)
        elif x > q3 and x <= q4:
            l1d.append(x)
        else:
            l1e.append(x)
    t1, t1a, t1b, t1c, t1d, t1e = len(l1), len(l1a), len(l1b), len(l1c), len(l1d), len(l1e)
    r1a, r1b, r1c, r1d, r1e = t1a/t1, t1b/t1, t1c/t1, t1d/t1, t1e/t1
    if t1b + t1d == 0:
        r_q_half = None
    elif t1c == 0:
        r_q_half = 1000
    else:
        r_q_half = (t1b + t1d) / t1c
    
    return t1b, t1c, t1d, r_q_half


def get_mut_avg(total_mut, tmplen):
    return total_mut / (tmplen / 1_000_000)


def check_mutdist(list_mos, dict_chrnum_mutlist, sex):
    for mos in list_mos:
        list_mut = list(filter(lambda x: x.pos>mos.start and x.pos<mos.end, dict_chrnum_mutlist[mos.chrnum]))
        
        # if sex == 'XY' and mos.chrnum in [23, 24] and mos.cnvtype == 'lossmos':
        # gainmos可能是完全相同的chr，则mut比例不能成为支持或否定证据
        if sex == 'XY' and mos.chrnum in [23, 24]:
            evi_mut = None
            total_mut_q2, total_mut_q3, total_mut_q4, ratio_q_half = None, None, None, None
        else:
            mut_r1, mut_r2 = model_mutratio(mos)
            total_mut_q2, total_mut_q3, total_mut_q4, ratio_q_half = stat_mut(list_mut, mut_r1, mut_r2)
            if total_mut_q2 + total_mut_q4 < 4:
                evi_mut = None
            elif ratio_q_half > 1 or total_mut_q2 + total_mut_q4 > 300:
                evi_mut = True
            elif get_mut_avg(total_mut_q3, mos.get_len()) < 2:
                evi_mut = None
            else:
                evi_mut = False
            
        info = {'total_mut_q2':total_mut_q2, 'total_mut_q3':total_mut_q3, 'total_mut_q4':total_mut_q4, 'ratio_q_half':ratio_q_half}
        mos.evi_mut = evi_mut
        mos.info.update(info)
    return list_mos

#============================
def get_prop_uncover(mos):
    a = mos.start
    step = 500_000
    t, t_c, t_nc = 0, 0, 0
    while a < mos.end + step:
        l = list(filter(lambda x: x.start>=a and x.start<a+step, mos.list_exon))

        t += 1
        if len(l) < 10:
            t_nc += 1
        else:
            t_c += 1
        
        a = a + step
    
    r = 0 if t == 0 else t_nc / t
    # print("prop_uncover", f"chr{mos.chrnum}", t, t_c, t_nc, r)
    return r


def get_prop_high(mos):
    t1 = len(list(filter(lambda x: x.diff>2.1+(mos.diff-1), mos.list_exon)))
    t = mos.get_exontotal()
    r = 0 if t==0 else t1/t
    return t1, r


def stat_reg(list_mos):
    for mos in list_mos:
        mos.prt_uncover = get_prop_uncover(mos)
        mos.total_high, mos.prt_high = get_prop_high(mos)
    return list_mos

# =========================================
# 该统计包含长度和细胞比例均仅在阈值附近的杂mos区段
def get_dataqual(total_reg1st):
    if total_reg1st <= 15:
        g = 1
    elif total_reg1st <= 50:
        g = 2
    # elif total_reg1st <= 100:
    elif total_reg1st <= 130:
        g = 3
    else:
        g = 4
    return g


def prod_funfilter(data_qual):
    # d = {"2-gainmos": 0.2, "4-gainmos": 0.2, "5-gainmos": 0.2, "6-gainmos": 0.2, "8-gainmos": 0.2, "13-gainmos": 0.2, "19-gainmos": 0.3, 
    #      "2-lossmos": 0.2, "4-lossmos":0.2, "5-lossmos":0.2, "6-lossmos":0.2, "8-lossmos":0.2, "13-lossmos": 0.2, "18-lossmos": 0.2, "19-lossmos": 0.3}
    d1 = {"19-lossmos":0.3, "23-lossmos":0.25}
    d2 = {"19-lossmos":0.4, "23-lossmos":0.25}
    
    def f1(mos):
        b = mos.get_len()>20_000_000 or (mos.get_len()>10_000_000 and mos.prt_cell>0.24)
        return b
    
    def f2(mos):
        b = (mos.get_len()>40_000_000 and mos.prt_cell>0.1) or (mos.get_len()>10_000_000 and mos.prt_cell>0.25) or (mos.get_len()>20_000_000 and mos.prt_cell>0.12)
        return b
    
    def f3(mos):
        b1 = (mos.get_len()>60_000_000 and mos.prt_cell>0.12) or (mos.get_len()>10_000_000 and mos.prt_cell>0.25)
        
        tmptag = f"{mos.chrnum}-{mos.cnvtype}"
        if tmptag in d1.keys():
            b2 = mos.prt_cell >= d1[tmptag]
        else:
            b2 = True
        return b1 and b2
    
    def f4(mos):
        b1 = (mos.get_len()>60_000_000 and mos.prt_cell>0.15) or (mos.get_len()>20_000_000 and mos.prt_cell>0.25) or (mos.get_len()>25_000_000 and mos.prt_cell>0.2)
        
        tmptag = f"{mos.chrnum}-{mos.cnvtype}"
        if tmptag in d2.keys():
            b2 = mos.prt_cell >= d2[tmptag]
        else:
            b2 = True
        return b1 and b2
    
    def f5(mos):
        b = whether_whole_chr(mos) and mos.prt_cell>0.49
        return b
    
    if data_qual == 1:
        f = f1
    elif data_qual == 2:
        f = f2
    elif data_qual == 3:
        f = f3
    elif data_qual == 4:
        f = f4
    else:
        f = f5
    return f


def prod_funfilter_chrY(data_qual):
    if data_qual == 1:
        f = lambda mos: mos.prt_cell<0.95 and ((mos.get_len()>=30_000_000 and mos.prt_cell>0.25) or (mos.get_len()<30_000_000 and mos.get_len()>12_000_000 and mos.prt_cell>0.3))
    elif data_qual == 2:
        f = lambda mos: mos.prt_cell<0.95 and (mos.get_len()>40_000_000 and mos.prt_cell>0.15)
    else:
        f = lambda mos: mos.prt_cell<0.95 and (mos.get_len()>40_000_000 and mos.prt_cell>0.25)
        
    return f


# 部分区域外显子稀少，易产生假阳性波动
def whether_sparsereg(mos):
    d = [PART(19, 19_500_000, 30_500_000)]
    b = False
    for p in d:
        unlap_len = p.get_unovlap_len(mos)
        if unlap_len is not None and unlap_len < 5_000_000:
            # print("whether-sparse", mos)
            # print("unlap-len", unlap_len)
            b = True
    return b


def whether_pass(l):
    l1 = list(filter(lambda x: x.chrnum not in [23, 24], l))
    # 去除chrY染色体的AZF区干扰
    l2 = list(filter(lambda x: x.chrnum == 23 or (x.chrnum == 24 and x.start>17_000_000), l))
    if len(l1) > 1 or (len(l1) > 0 and len(l2) > 0):
        b = False
    else:
        b = True
    # print("info-whether_pass", len(l1), len(l2), b)
    return b


def filter_mos(list_mos, sex, bool_mut, total_reg1st):
    l = list_mos.copy()
    bool_1stfilter = True
    data_quality = get_dataqual(total_reg1st)
    
    while bool_1stfilter or not whether_pass(l):
        if not bool_1stfilter:
            data_quality += 1
        bool_1stfilter = False
        if data_quality >= 6:
            if not no_dataqual_filter:
                l = list(filter(lambda x: x.prt_cell > 0.45, l))
            break
        
        fun_filter1 = prod_funfilter(data_quality)
        fun_filter1_chrY = prod_funfilter_chrY(data_quality)
        
        list_mos_t = l.copy()
        l = []
        for mos in list_mos_t:
            b1 = sex != "XX" or mos.chrnum != 24
            b2 = (mos.chrnum != 24 and fun_filter1(mos)) or (mos.chrnum == 24 and fun_filter1_chrY(mos))
            b3 = not whether_sparsereg(mos)
            b4 = mos.chrnum == 24 or (mos.prt_uncover is not None and mos.prt_uncover <= 0.7)
            b5 = mos.prt_high < 0.002 or mos.total_high<5 or (mos.get_len()>100_000_000 and mos.prt_high<0.0035) or (mos.prt_cell>0.25 and mos.prt_high<0.0035)
            b6 = mos.info['total_samedist_mid']<=1 and mos.info['total_samedist_short']<=10
            b7 = mos.evi_exon is None or mos.evi_exon
            b8 = not bool_mut or mos.evi_mut is None or mos.evi_mut or (sex == 'XY' and mos.chrnum in [23, 24])
            #if str(mos.chrnum) == "5":
            #    print(mos.chrnum, mos.start, mos.end, b1, b2, b3, b4,b5,b6,b7,b8)
            if b1 and b2 and b3 and b4 and b5 and b6 and b7 and b8:
                l.append(mos)
    
    return l, data_quality


# 避免不同区间间的重合
def fix_end(list_mos):
    l = []
    for mos in list_mos:
        for mos1 in l:
            if mos.chrnum == mos1.chrnum:
                if mos.start < mos1.end and mos.start > mos1.start and mos.end > mos1.end:
                    if mos1.end - mos.start < 500_000:
                        mos.start = mos1.end + 1
                if mos.end < mos1.end and mos.end > mos1.start and mos.start < mos1.start:
                    if mos.end - mos1.start < 500_000:
                        mos.end = mos1.start - 1
        l.append(mos)
    return l
# =========================================
def tr_key_none(d, k):
    a = d.get(k, 0)
    a = 0 if a is None else a
    return a


def output_mos(list_mos, sex, total_reg1st, data_quality):
    info1 = f"##sex:{sex}"
    # info2 = f"##STEP1reg-total:{total_reg1st}##RATIO_gainloss:{ratio_gainloss}##DATAqual:{data_quality}"
    info2 = f"##STEP1reg-total:{total_reg1st}##DATAqual:{data_quality}"
    header = "#chr\tstart\tend\tcnvtype\tdiff\tprt_cell\tprt_exon_inmos\ttotal_exon\tratio_mutinmos_half\ttotal_mut10-40\ttotal_mut40-60\ttotal_mut60-90\ttotal_samedist_short\ttotal_samedist_mid\ttotal_samedist_long"
    buffer_o = [info1, info2, header]

    for mos in list_mos:
        s1 = "{}\t{}\t{}\t{}\t{:.2f}\t{:.2f}".format(chrnum_tochr(mos.chrnum), mos.start, mos.end, mos.cnvtype, mos.diff, mos.prt_cell)
        s2 = "{:.2f}\t{}".format(tr_key_none(mos.info, 'prt_exon_f'), tr_key_none(mos.info, 'total_exon'))
        s3 = "{:.2f}\t{}\t{}\t{}".format(tr_key_none(mos.info, 'ratio_q_half'), tr_key_none(mos.info, 'total_mut_q2'), tr_key_none(mos.info, 'total_mut_q3'), tr_key_none(mos.info, 'total_mut_q4'))
        s4 = "{}\t{}\t{}".format(tr_key_none(mos.info, 'total_samedist_short'), tr_key_none(mos.info, 'total_samedist_mid'), tr_key_none(mos.info, 'total_samedist_long'))

        s = f"{s1}\t{s2}\t{s3}\t{s4}"
        
        buffer_o.append(s)

    return buffer_o

# =========================================
# =========================================
#  python3 0-script/find-mos.py /share/chg1fs1b/prod/project/DDN22022000-22022999/DDN22022177/AHJ833_NTb01F_NT01T_7321/ExonResult/exons_details.txt r1/AHJ833.txt -v /share/chg1fs1b/prod/project/DDN22022000-22022999/DDN22022177/AHJ833_NTb01F_NT01T_7321/Variants/variants.vcf.gz

# python3 /share/chg1fs1b/train/liuyb/65-wesmos/0-script/find-mos.py /share/chg1fs1b/prod/project/DDN23005000-23005999/DDN23005258/AIV981_NT01F_NT01T_7708/ExonResult/exons_details.txt /share/chg1fs1b/train/liuyb/65-wesmos/r1/DDN23005258-AIV981-mos.txt -v /share/chg1fs1b/prod/project/DDN23005000-23005999/DDN23005258/AIV981_NT01F_NT01T_7708/Variants/variants.vcf.gz

intro_str = ""

parser = argparse.ArgumentParser(description=intro_str)

parser.add_argument('exon_detail', type=str, help='')
parser.add_argument('f_o', type=str, help='')
parser.add_argument('-v', '--mut_gvcf', type=str, help='')
parser.add_argument('--verbose', action='store_true', help='')
parser.add_argument('--nofilter', action='store_true', help='not filter based on data qual')

args = parser.parse_args()

f_exondetail = args.exon_detail
f_o = args.f_o
gvcf_mut = args.mut_gvcf
verbose = args.verbose
no_dataqual_filter = args.nofilter

bool_mut = gvcf_mut is not None
#===================

print_log("mosaic find::START")

print_log("mosaic find::read data")
dict_chrnum_exonlist = read_exondetail(f_exondetail)
sex = whichsex(dict_chrnum_exonlist)
if bool_mut:
    dict_chrnum_mutlist = read_gvcfmut(gvcf_mut)

# stat_chr(dict_chrnum_exonlist, prop='start', fun_avg=fun_avg1)
# stat_chr(dict_chrnum_mutlist, prop='pos', fun_avg=fun_avg2)

print_log("mosaic find::filter bin and merge")
list_mos = filter_merge(dict_chrnum_exonlist, sex)
total_reg1st, ratio_gainloss = stat_1st(list_mos)
if verbose:
    print_listlog("|| list_mos1", list_mos)

print_log("mosaic find::calibrate breadpoint")
list_mos = calib_bkp(list_mos, dict_chrnum_exonlist, sex)

print_log("mosaic find::stat region")
list_mos = stat_samedist(list_mos)
list_mos = check_exondist(list_mos, dict_chrnum_exonlist, sex)

if bool_mut:
    list_mos = check_mutdist(list_mos, dict_chrnum_mutlist, sex)

list_mos = stat_reg(list_mos)

if verbose:
    print_listlog("|| list_mos2", list_mos)

print_log("mosaic find::filter region and fix end")
#print(len(list_mos))
#for i in list_mos:
#    print(i.chrnum, i.start, i.end)
list_mos, data_quality = filter_mos(list_mos, sex, bool_mut, total_reg1st)
list_mos = fix_end(list_mos)
if verbose:
    print_listlog("|| list_mos3", list_mos)
print(f"||| STEP1reg-total::{total_reg1st}::ratio_gainloss::{ratio_gainloss}::DATA quality::{data_quality}")

print_log("mosaic find::write results")
buffer_o = output_mos(list_mos, sex, total_reg1st, data_quality)
write_tofile(f_o, buffer_o)

print_log("mosaic find::END")

