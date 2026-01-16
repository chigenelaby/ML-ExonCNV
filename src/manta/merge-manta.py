#!/user/bin/python
# -*- coding: utf-8 -*-
import os, sys, argparse
import collections, math, statistics
import glob, gzip, subprocess, json
from time import strftime, gmtime, localtime
from io import StringIO
import sqlite3

# =========================================
def print_log(msg):
    c_time = strftime("%Y-%m-%d %H:%M:%S", localtime())
    print(f">>> {c_time}\t{msg}")


def write_tofile(f, buffer_o):
    with open(f, 'w') as outF:
        for o in buffer_o:
            outF.write(f"{o}\n")


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


def chrnum_tochr(i):
    if i == 23:
        s = "X"
    elif i == 24:
        s = "Y"
    else:
        s = str(i)
    return f"chr{s}"

# =========================================
class REG():
    chrnum: int
    start: int
    end: int

    def __init__(self, chrnum, start, end):
        self.chrnum = chrnum
        self.start = start
        self.end = end

    def __repr__(self):
        return "{}-{}-{}".format(self.chrnum, self.start, self.end)
    
    def whether_ovlap(self, a):
        return self.chrnum == a.chrnum and self.start <= a.end and a.start <= self.end
    
    def whether_cover(self, a):
        return self.chrnum == a.chrnum and self.start <= a.start and a.end <= self.end
    
    def get_len(self):
        return self.end - self.start + 1

    def cal_ovlap_len(self, a):
        tmplen_o = 0
        if self.whether_ovlap(a):
            start_o = max(self.start, a.start)
            end_o = min(self.end, a.end)
            tmplen_o = end_o - start_o + 1
        return tmplen_o

    def cal_ovlap_prop(self, a):
        p = self.cal_ovlap_len(a) / self.get_len()
        return p


class REG1(REG):
    diff: float
    def __init__(self, chrnum, start, end):
        super().__init__(chrnum, start, end)
        self.diff = float('nan')

    def __repr__(self):
        return "{}-{}".format(super().__repr__(), self.diff)


class EXON(REG):
    exon_name: str

    def __init__(self, chrnum, start, end, exon_name):
        super().__init__(chrnum, start, end)
        self.exon_name = exon_name

    def __repr__(self):
        return "{}-{}".format(super().__repr__(), self.exon_name)


class GENE(REG):
    gene_name: str
    var_name: str
    direction: str
    list_exon: list

    def __init__(self, chrnum, start, end, gene_name, var_name, direction):
        super().__init__(chrnum, start, end)
        self.gene_name = gene_name
        self.var_name = var_name
        self.direction = direction
        self.list_exon = []

    def __repr__(self):
        return "{}-{}-{}={}".format(super().__repr__(), self.gene_name, self.var_name, ",".join([str(x) for x in self.list_exon]))


class CNV(REG):
    cnvtype: str
    freq: float
    list_geneovlap: list
    list_regcover: list
    list_regtwoside: list

    # 若恰好是一个内含子，疑似RNA插入
    bool_intron: bool
    diff_avg: float
    bool_dpevi: bool

    def __init__(self, chrnum, start, end, cnvtype):
        super().__init__(chrnum, start, end)
        self.cnvtype = cnvtype
        self.freq = 0
        self.list_geneovlap = []
        self.list_regcover = []
        self.list_regtwoside = []
        self.bool_intron = False
        self.diff = float('nan')
        self.diff_avg = float('nan')    #yzh
        self.bool_dpevi = False

    def __repr__(self):
        # return "{}-{}||{}-{}||{}".format(super().__repr__(), self.cnvtype, "|".join([str(x) for x in self.list_regcover]), "|".join([str(x) for x in self.list_regtwoside]), "|".join([str(x) for x in self.list_geneovlap]))
        return "{}-{}-{}|intron|{}-{}|dp|{}-{}-{:.2f}-{}".format(super().__repr__(), self.cnvtype, self.freq, len(self.list_geneovlap), self.bool_intron, len(self.list_regcover), len(self.list_regtwoside), self.diff_avg, self.bool_dpevi)
    
    def cal_ovlap_prop(self, a):
        p = 0
        if self.cnvtype == a.cnvtype:
            p = super().cal_ovlap_prop(a)
        return p

    def toid(self):
        return "{}={}={}={}".format(self.chrnum, self.start, self.end, self.cnvtype)


def id_tocnv(tmpid):
    l = tmpid.split("=")
    return CNV(int(l[0]), int(l[1]), int(l[2]), l[3])


# =========================================
def read_sex(f_sex):
    s = "XY"
    with open(f_sex) as inF:
        line = inF.readline()
        s1 = line.rstrip("\n").split("\t")[1]
        s = "XY" if "Y" in s1 else "XX"
    return s

# =========================================
def read_db(db_exon):
    d_chrnum_genes = collections.defaultdict(list)
    d_gene = {}
    
    with open(db_exon) as inF:
        for line in inF:
            if line.startswith("#"):
                continue
    
            lparts = line.rstrip("\n").split("\t")
            chrnum = chr_tonum(lparts[0])
            start = int(lparts[1])
            end = int(lparts[2])
            gene_name = lparts[3]
            direction = lparts[4]
            var_name = lparts[5]
            gene_start = int(lparts[6])
            gene_end = int(lparts[7])
            exon_name = lparts[8]

            exon = EXON(chrnum, start, end, exon_name)
            if chrnum == 0:
                continue
            if gene_name in d_gene.keys():
                d_gene[gene_name].list_exon.append(exon)
            else:
                gene = GENE(chrnum, gene_start, gene_end, gene_name, var_name, direction)
                gene.list_exon.append(exon)
                d_chrnum_genes[chrnum].append(gene)
                d_gene[gene_name] = gene

    # # 查看内含子的大小
    # for gene in d_gene.values():
    #     list_exon = gene.list_exon
    #     set_chrnum = set([x.chrnum for x in list_exon])
    #     for chrnum in set_chrnum:
    #         l1 = list(filter(lambda x: x.chrnum == chrnum, list_exon))
    #         l1.sort(key=lambda x: x.start)

    #         for i in range(1, len(l1)):
    #             d = l1[i].start - l1[i-1].end
    #             print(chrnum, l1[i].start, l1[i].end, d)

    return d_chrnum_genes, d_gene

# =========================================
def tr_cnvtype(s):
    if s.startswith("MantaDUP"):
        s1 = "gain"
    elif s.startswith("MantaDEL"):
        s1 = "loss"
    else:
        s1 = None
    return s1


def get_endinfo(s):
    s1 = None
    for a in s.split(";"):
        if a.startswith("END="):
            s1 = int(a[4:])
    return s1


def read_mantavcf(f_manta):
    d = collections.defaultdict(list)
    with open(f_manta) as inF:
        for line in inF:
            if line.startswith("#"):
                continue

            if not line.startswith("chr"):
                continue
 
            lparts = line.rstrip("\n").split("\t")

            chrnum = chr_tonum(lparts[0])
            start = int(lparts[1])
            cnvtype = tr_cnvtype(lparts[2])

            if cnvtype is None or chrnum == 0:
                continue

            end = get_endinfo(lparts[7])
            
            if end is None:
                print("WARNING::NO END info::{}".format(line.rstrip()))
                continue

            cnv = CNV(chrnum, start, end, cnvtype)
            d[chrnum].append(cnv)

    return d


# =========================================
def read_detail(f):
    d = collections.defaultdict(list)
    with open(f) as inF:
        for line in inF:
            if line.startswith("#"):
                continue
    
            lparts = line.rstrip("\n").split("\t")
            chrnum = chr_tonum(lparts[0])
            start = int(lparts[1])
            end = int(lparts[2])
            diff = float(lparts[5])
            if chrnum == 0:
                continue
            reg = REG1(chrnum, start, end)
            reg.diff = diff
            d[chrnum].append(reg)
    return d


# =========================================
def get_ovlap(x, l_regs):
    l = []
    for reg in l_regs:
        if x.whether_ovlap(reg):
            l.append(reg)

    # if len(l) > 0:
    #     print(l)
    return l

def get_cover(x, l_regs):
    l = []
    for reg in l_regs:
        if x.whether_cover(reg):
            l.append(reg)

    # if len(l) > 0:
    #     print(l)
    return l

# =========================================
def cut_intron(l_e):
    l_e.sort(key=lambda x: x.start)
    
    l = []
    if len(l_e) < 2:
        return l
    
    for i in range(1, len(l_e)):
        l.append((l_e[i-1].end, l_e[i].start))
    return l


def cmp_intropos(s, e, list_intronpair):
    tmpbool = False
    for a, b in list_intronpair:
        if abs(a-s) <= 3 and abs(b-e) <= 3:
            tmpbool = True
            break
    return tmpbool


def whether_intron(mantacnv, list_gene):
    # 可能是重叠基因，所以是list
    l_bool = []
    for g in list_gene:
        # 同一个基因可能在不同染色体上分布
        list_exon1 = list(filter(lambda x: x.chrnum == mantacnv.chrnum, g.list_exon))
        list_intronpair = cut_intron(list_exon1)
        # print(g, list_intronpair)
        l_bool.append(cmp_intropos(mantacnv.start, mantacnv.end, list_intronpair))

    return any(l_bool)


def determine_intron(dict_chrnum_mantacnv, dict_chrnum_dbgenes):
    for i in dict_chrnum_mantacnv.keys():
        for mantacnv in dict_chrnum_mantacnv[i]:
            mantacnv.list_geneovlap = get_ovlap(mantacnv, dict_chrnum_dbgenes[i])
            mantacnv.bool_intron = whether_intron(mantacnv, mantacnv.list_geneovlap)
    return dict_chrnum_mantacnv

# =========================================
def cal_diffavg(l_reg):
    l_diff = [x.diff for x in l_reg]
    l_diff = list(filter(lambda x: x not in [float('nan'), float('inf'), float('-inf')], l_diff))

    if len(l_diff) == 0:
        return float('nan')
    
    avg_diff = sum(l_diff) / len(l_diff)
    return avg_diff


def whether_dpevi(cnvtype, a):
    b = False
    if a in [float('nan'), float('inf'), float('-inf')]:
        return b

    if cnvtype == "gain" and a >= 1.2:
        b = True
    elif cnvtype == "loss" and a <= 0.8:
        b = True
    else:
        pass

    return b


def determine_dp(dict_chrnum_mantacnv, dict_chrnum_regdetail):
    for i in dict_chrnum_mantacnv.keys():
        for mantacnv in dict_chrnum_mantacnv[i]:
            l_cover = get_cover(mantacnv, dict_chrnum_regdetail[i])
            l_ovlap = get_ovlap(mantacnv, dict_chrnum_regdetail[i])
            # print(l_cover)
            # print(l_ovlap)
            
            mantacnv.list_regcover = l_cover
            mantacnv.list_regtwoside = list(set(l_ovlap) - set(l_cover))
            # print(mantacnv.list_regtwoside)

            mantacnv.diff_avg = cal_diffavg(mantacnv.list_regcover)
            if math.isnan(mantacnv.diff_avg) and len(mantacnv.list_regtwoside) != 0:
                for a in mantacnv.list_regtwoside:
                    # print(mantacnv.get_len())
                    # print(a.get_len())
                    r1 = REG.cal_ovlap_prop(mantacnv, a)
                    r2 = REG.cal_ovlap_prop(a, mantacnv)
                    # print(r1, r2)

            mantacnv.bool_dpevi = whether_dpevi(mantacnv.cnvtype, mantacnv.diff_avg)
            # print(mantacnv)

    return dict_chrnum_mantacnv

# =========================================
def cal_total(cursor, sam, sex):
    sql_checkin = f"SELECT COUNT(*) FROM INFO WHERE sam = ?;"
    cursor.execute(sql_checkin, (sam,))
    num_in = cursor.fetchone()[0]
    bool_in = False if num_in == 0 else True

    sql_caltotal = f"SELECT COUNT(*) FROM INFO WHERE sex = ?;"
    cursor.execute(sql_caltotal, (sex,))
    num_rows = cursor.fetchone()[0]

    # print("In db_freq {} {}, total {}".format(sex, bool_in, num_rows))
    return bool_in, num_rows


def cal_freq(cursor, cnv, sex, bool_in, total_sam, delta=3):
#     sql_s_raw = """
# SELECT COUNT(*) FROM {table_name}
# WHERE chrnum = ? AND 
# start >= ? AND 
# start <= ? AND 
# end >= ? AND
# end <= ?
# """
    sql_s_raw = """
SELECT COUNT(DISTINCT sam) FROM {table_name}
WHERE chrnum = ? AND 
start >= ? AND 
start <= ? AND 
end >= ? AND
end <= ?
"""
    sql_s = sql_s_raw.format(table_name=sex)
    cursor.execute(sql_s, (cnv.chrnum, cnv.start-delta, cnv.start+delta, cnv.end-delta, cnv.end+delta))
    total_s = cursor.fetchone()[0]

    total_s = total_s - 1 if bool_in else total_s
    total_sam = total_sam - 1 if bool_in else total_sam
    if total_sam < 50:
        freq = 0
    else:
        freq = total_s / total_sam
    # print("Freq {}, total_s {}".format(freq, total_s))
    freq = round(freq, 3)
    return freq


def update_freq(dict_chrnum_mantacnv, db_freq, sam, sex):
    conn = sqlite3.connect(db_freq)
    cursor = conn.cursor()

    bool_in, total_sam = cal_total(cursor, sam, sex)

    for chrnum in dict_chrnum_mantacnv:
        for cnv in dict_chrnum_mantacnv[chrnum]:
            cnv.freq = cal_freq(cursor, cnv, sex, bool_in, total_sam)
            # print_log(cnv)

    conn.close

    return dict_chrnum_mantacnv

# =========================================
def filter_split_cnvgroup(dict_chrnum_mantacnv):
    list_highfreq = []
    list_rm = []
    list_add = []

    for i in dict_chrnum_mantacnv.keys():
        for cnv in dict_chrnum_mantacnv[i]:
            if str(cnv.chrnum) == "22":
                print(cnv, cnv.freq, cnv.bool_intron, cnv.bool_dpevi)
            if cnv.freq >= 0.05:
                # print(f"info::cnv-remove-freq>0.05::{cnv}")
                list_highfreq.append(cnv)
            elif cnv.bool_intron or not cnv.bool_dpevi:
                list_rm.append(cnv)
            else:
                list_add.append(cnv)
    return list_highfreq, list_rm, list_add

# =========================================
# l_rm中收录的皆为内含子区域的cnv，所以不需要考虑复杂情形
# l_rm另收录无外显子支持的manta结果（废弃）
# def whether_rm(cnv, l_rm):
#     b = False
#     for cnv1 in l_rm:
#         if cnv.cal_ovlap_prop(cnv1) > 0.98:
#             b = True
#             break
#     return b


def whether_fix(cnv, l_add, dict_fix, lparts):
    b = False
    for cnv1 in l_add:
        if cnv.cal_ovlap_prop(cnv1) > 0.4:
            b = True
            dict_fix[cnv1.toid()].append(lparts)
            break
        elif cnv.cal_ovlap_prop(cnv1) > 0:
            print("info:cnvovlap-manta-wesexon:", cnv1, cnv)
            # pass
        else:
            pass

    return b, dict_fix


def to_o(f_in_fmtgene, l_add):
    buffer_o = []
    # 收集输入文件中需要更改断点的行
    dict_fix = collections.defaultdict(list)
    total_col = 0
    bool_XY = False

    with open(f_in_fmtgene) as inF:
        for line in inF:
            if line.startswith("#"):
                buffer_o.append(line.rstrip("\n"))

                if line.startswith("#chr"):
                    total_col = len(line.rstrip("\n").split("\t"))

                continue

            if line.startswith("chrY"):
                bool_XY = True
    
            lparts = line.rstrip("\n").split("\t")
            if chr_tonum(lparts[0]) == 0:
                continue
# #chr    start   end     gene_name       gene_info_str   best_exon_str   freq    all_freq        infos   ALX488_exon     ALX488_tag      ALX488_result
# chr1    368582  368702  OR4F29  (REG:2) NA      0.103   0.2     {"reliability_show": ["满足外显子个数条件"], "report_tag": "NA"}        0/1=1.4966      gain1   不可靠
            cnv = CNV(chr_tonum(lparts[0]), int(lparts[1]), int(lparts[2]), lparts[10][0:4])
            # if not hasattr(cnv, 'diff_avg'):
            #     print(1)
            #     continue  # 新bed区间不覆盖，直接跳过不输出 yzh
            # bool_rm = whether_rm(cnv, l_rm)

            # if bool_rm:
            #     # print("info:cnv-remove:intron:", cnv)
            #     continue

            bool_fix, dict_fix = whether_fix(cnv, l_add, dict_fix, lparts)

            if bool_fix:
                # print("info:cnv-fix:breakpoint:", cnv)
                continue


            buffer_o.append("\t".join(lparts))

    #怎么写的loss1比例 为什么设定diffavg NAN 会变成NA yzh
    for cnv in l_add:
        
        cnvid = cnv.toid() 
        if cnvid in dict_fix.keys():
            lparts = dict_fix[cnvid][0]
            if lparts[11] == "不可靠":
                lparts[11] = "相对可靠"
        else:
            lparts = ["NA"] * total_col
            lparts[0] = chrnum_tochr(cnv.chrnum)
            lparts[8] = "{\"reliability_show\": [], \"report_tag\": \"NA\"}"
            if bool_XY and cnv.chrnum in [23, 24]:
                lparts[10] = f"{cnv.cnvtype}"
            else:
                lparts[10] = f"{cnv.cnvtype}1"
            lparts[11] = "不可靠"
        
        lparts[1] = str(cnv.start)
        lparts[2] = str(cnv.end)

        # manta结果，添加标注
        dict_info = json.loads(lparts[8], encoding='utf-8')
        if "reliability_show" not in dict_info.keys():
            dict_info["reliability_show"] = []

        if cnvid in dict_fix.keys() or cnv.bool_dpevi:
            dict_info["reliability_show"].append("断点外显子支持")
            # dict_info["reliability_show"].append("疑似RNA插入")
        else:
            dict_info["reliability_show"].append("断点外显子不支持")

        lparts[8] = json.dumps(dict_info, ensure_ascii=False)

        buffer_o.append("\t".join(lparts))
    return buffer_o

# =========================================
# =========================================
# python3 /share/chg1fs1b/train/liuyb/60-wes-exon/0-script1/merge-manta.py /share/chg1fs1b/prod/project/DDN24002000-24002999/DDN24002710/ALX488_NT01_NT01T_8500/Manta/diploidSV.vcf /share/cg01-08-fsd/public/liuyb/60-wes-exon/r-t/DDN24002710/ALX488/ExonResult/exons_details.txt /share/cg01-08-fsd/public/liuyb/60-wes-exon/r-t/DDN24002710/ALX488/ExonResult/qc_karyotype.txt /share/cg01-08-fsd/public/liuyb/60-wes-exon/r-t/DDN24002710/ALX488/ExonResult/format_gene_info_pre.txt format_gene_info-manta.txt --sam ALX488 --db_exon /share/chg1fs1b/train/liuyb/30-cnvseq-new/0-script/db/EXON.longestNM.bed --db_freq /share/chg1fs1b/train/liuyb/60-wes-exon/0-db/wesmanta.db

# 添加的manta结果，在结果文件中添加标注"断点外显子支持", "断点外显子不支持"
intro_str = "将manta结果合并到exon-gene结果中，1.若是内含子对应区段，则删除；2. 若没有外显子支持则删除"

parser = argparse.ArgumentParser(description=intro_str)

parser.add_argument('f_manta', type=str, help='')
parser.add_argument('f_detail', type=str, help='')
parser.add_argument('f_sex', type=str, help='')
parser.add_argument('f_in_fmtgene', type=str, help='')
parser.add_argument('f_o_fmtgene', type=str, help='')
parser.add_argument('--sam', type=str, help='')
parser.add_argument('--db_exon', type=str, help='')
# sqlite3 db
parser.add_argument('--db_freq', type=str, help='')

args = parser.parse_args()

f_manta = args.f_manta
f_detail = args.f_detail
f_sex = args.f_sex
f_in_fmtgene = args.f_in_fmtgene
f_o_fmtgene = args.f_o_fmtgene
sam = args.sam
db_exon = args.db_exon
db_freq = args.db_freq

#=========================================
sex = read_sex(f_sex)

dict_chrnum_dbgenes, dict_dbgene = read_db(db_exon)
# for a in dict_chrnum_dbgenes[1][0:10]:
#     print(a)
#     print(dict_dbgene[a.gene_name])

dict_chrnum_regdetail = read_detail(f_detail)
# print(len(dict_regdetail.keys()))
# print(dict_regdetail[1][0:10])

dict_chrnum_mantacnv = read_mantavcf(f_manta)

dict_chrnum_mantacnv = determine_intron(dict_chrnum_mantacnv, dict_chrnum_dbgenes)

dict_chrnum_mantacnv = determine_dp(dict_chrnum_mantacnv, dict_chrnum_regdetail)
# for i in dict_chrnum_mantacnv.keys():
#     print(dict_chrnum_mantacnv[i])
# print(dict_chrnum_mantacnv[1][0])

if os.path.exists(db_freq):
    dict_chrnum_mantacnv = update_freq(dict_chrnum_mantacnv, db_freq, sam, sex)

l_hfreq, l_rm, l_add = filter_split_cnvgroup(dict_chrnum_mantacnv)
# print("rm-num", len(l_rm), l_rm)
total_mantacnv = len(l_hfreq) + len(l_rm) + len(l_add)
print_log("merge-mantacnv: {}: total {}; hfreq {}; rm-num {}; add-num {}".format(sam, total_mantacnv, len(l_hfreq), len(l_rm), len(l_add)))
# for a in l_rm + l_add:
#     print(a)

buffer_o = to_o(f_in_fmtgene, l_add)

write_tofile(f_o_fmtgene, buffer_o)
