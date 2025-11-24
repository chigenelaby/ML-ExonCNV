#!/usr/bin/python
#-*- coding: utf-8 -*-

from __future__ import print_function
from __future__ import unicode_literals
from __future__ import absolute_import
from __future__ import division

import sys
from io import StringIO

#===============================
# 修改核型标注

#===============================

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

list_chr_start, list_chr_end = get_chr_startend('cut')


#===============================

def chr_num(tmpchr):
    tmpchr = tmpchr.strip()

    if tmpchr.startswith("NC"):
        return int(tmpchr[7:9])
    elif tmpchr.startswith("chr"):
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


def split_nultype(nultype):
    tmp_list = nultype.split(",")

    if len(tmp_list) > 2:
        chr_total = tmp_list[0]
        sexchr = tmp_list[1]
        list_normalchr_change = tmp_list[2:]
    else:
        chr_total = tmp_list[0]
        sexchr = tmp_list[1]
        list_normalchr_change = []
    return chr_total, sexchr, list_normalchr_change


def cat_nultype(chr_total, sexchr, list_normalchr_change):
    str_normalchr_change = ",".join(list_normalchr_change)
    nultype = ""

    if str_normalchr_change != "":
        nultype = ",".join([chr_total, sexchr, str_normalchr_change])
    else:
        nultype = ",".join([chr_total, sexchr])

    return nultype


# ===============================
# 嵌合体对染色体端粒附近判断不好，放宽阈值
def whether_whole_chr(chrnum, start, end, start_thld=3_000_000, end_thld=3_000_000):
    chr_startpos_thld = list_chr_start[chrnum - 1] + start_thld
    chr_endpos_thld = list_chr_end[chrnum - 1] - end_thld

    tmpbool = False

    if start < chr_startpos_thld and end > chr_endpos_thld:
        tmpbool = True
    # print(chrnum, start, end)
    # print(chr_startpos_thld, chr_endpos_thld)
    return tmpbool


def get_chrchange(mosfile):
    total_chrchange = 0
    # ["+2", "-4"]
    list_normalchr_change = []
    # ["+23", "-24"]
    list_sexchr_change = []

    with open(mosfile) as inF:
        for line in inF:
            if line.startswith("NC") or line.startswith("chr"):
                list_lineparts = line.strip().split("\t")
                tmpchr_num = chr_num(list_lineparts[0])
                tmpstart = int(list_lineparts[1])
                tmpend = int(list_lineparts[2])
                tmpmostype = list_lineparts[3]

                # print(tmpchr_num)
                # print(tmpstart, list_chr_start[tmpchr_num - 1])
                # print(tmpend, list_chr_end[tmpchr_num - 1])
                # print(tmpmostype)

                if whether_whole_chr(tmpchr_num, tmpstart, tmpend):
                    if "gainmos" in tmpmostype:
                        # print("check_code")
                        total_chrchange += 1

                        if tmpchr_num < 23:
                            list_normalchr_change.append("+" + str(tmpchr_num))
                        else:
                            list_sexchr_change.append("+" + str(tmpchr_num))

                    elif "lossmos" in tmpmostype:
                        total_chrchange -= 1

                        if tmpchr_num < 23:
                            list_normalchr_change.append("-" + str(tmpchr_num))
                        else:
                            list_sexchr_change.append("-" + str(tmpchr_num))

                    else:
                        pass

    return total_chrchange, list_sexchr_change, list_normalchr_change


def update_mos_nultype(nultype_ori, sex, total_chrchange, list_sexchr_change, list_normalchr_change):
    # print(list_sexchr_change)
    # print(list_normalchr_change)
    if "Y" in sex:
        nultype_normal = "46,XY"
    else:
        nultype_normal = "46,XX"


    if list_normalchr_change == [] and list_sexchr_change == []:
        return nultype_ori
    elif list_normalchr_change == [] and set(list_sexchr_change) == set(["+23", "-24"]):
        l = nultype_ori.split(",")
        l[1] = "XX"
        a1 = ",".join(l)
        l[1] = "XY"
        a2 = ",".join(l)
        return f"mos {a1}/{a2}"
    else:
        chr_total_ori, sexchr_ori, list_normalchr_change_ori = split_nultype(nultype_normal)

        chr_total_alt = str(int(chr_total_ori) + int(total_chrchange))
        
        sexchr_alt = sexchr_ori
        if list_sexchr_change != []:
            for tmpsexchr_change in list_sexchr_change:
                if "-" == tmpsexchr_change[0]:
                    if "23" == tmpsexchr_change[1:] and "X" == sexchr_alt[0]:
                        sexchr_alt = sexchr_alt[1:]
                    if "24" == tmpsexchr_change[1:] and "Y" == sexchr_alt[-1]:
                        sexchr_alt = sexchr_alt[:-1]
                if "+" == tmpsexchr_change[0]:
                    if "23" == tmpsexchr_change[1:]:
                        sexchr_alt = "X" + sexchr_alt
                    if "24" == tmpsexchr_change[1:]:
                        sexchr_alt = sexchr_alt + "Y"

        list_normalchr_change_alt = list_normalchr_change
        list_normalchr_change_alt.extend(list_normalchr_change_ori)

        list_normalchr_change_alt.sort(key=lambda x: int(x[1:]))

        nultype_alt = cat_nultype(chr_total_alt, sexchr_alt, list_normalchr_change_alt)

        if nultype_ori == nultype_alt:
            nultype_mos = f"mos {nultype_normal}/{nultype_alt}"
        else:
            nultype_mos = f"mos {nultype_ori}/{nultype_alt}"

        return nultype_mos



# ===============================

if __name__ == "__main__":
    str_intro = """
    Add the ISCN notes of mos, such as mos 46,XX/47,XX,+5

    python fix-sex.py [input_sexfile] [input_mosfile] [output_sexfile]

    [input_sexfile] ExonResult/qc_karyotype.txt
    [input_mosfile] ExonResult/mos.txt
    [output_sexfile] ExonResult/qc_karyotype_mos.txt
"""
    if len(sys.argv) != 4:
        print(str_intro)
        sys.exit()


    input_sexfile = sys.argv[1]
    input_mosfile = sys.argv[2]
    output_sexfile = sys.argv[3]


    total_chrchange, list_sexchr_change, list_normalchr_change = get_chrchange(input_mosfile)

    with open(input_sexfile) as inF, open(output_sexfile, 'w') as outF:
        for i, line in enumerate(inF):
            if i == 0:
                nultype, sex = line.strip("\n").split("\t")[0:2]

        if nultype not in ["46,XX", "46,XY"]:
            nultype2 = nultype
        else:
            nultype2 = update_mos_nultype(nultype, sex, total_chrchange, list_sexchr_change, list_normalchr_change)

        tmps = f"{nultype2}\t{sex}"
        outF.write(tmps)