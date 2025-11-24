import sys
import os
import re
import math
from operator import add
#合并的区间拆分为单个外显子以其特征如balance,gc,z_score等
def compare(a1,a2,b1,b2): #比较区间是否有交集
  Max  = max([a1,b1])
  Min  = min([a2,b2])
  if(Min - Max > 0):
    return 1
  else:
    return 0


def compareinclued(a1,a2,b1,b2):
    if (a1 <= b1 and a2 >= b2):
        return 1
    else:
        return 0

def compareregion(a1,a2,b1,b2):
    Max  = max([a1,b1])
    Min  = min([a2,b2])
    temp = 0
    if(Min - Max > 0):
       temp = Min-Max
       return temp
    else:
       return temp

def revisregion(Chr,start,end,G):
    transregion = {}
    matchregion = ""
    for key,values in G.items():
        cchr,cst,sen = key.split("_")
        cst = int(cst)
        sen = int(sen)
        length = compareregion(start,end,cst,sen)
        if(Chr == cchr and length):
             transregion[key] = length
    #print(Chr,start,end,transregion)
    for key,value in sorted(transregion.items(), key=lambda x:x[1]):
       matchregion = key
    return matchregion
#gms   gcmappabilityscore.txt
#F     format_gene_info.txt

def diff_dist_para(diff,dist,posdistpara = 0.001,negdistpara = 0.05,const = 0.5):
    natrual = math.e
    diff = float(diff)
    dist = float(dist)
    if diff == 0 and dist == 0:
        return 0
    if dist>0:
        newdist = math.log(posdistpara*dist)
    else:

        newdist = negdistpara*dist
    para = diff / (const + natrual**(newdist))

    # print(diff,dist,newdist,para)
    return para


# def calculate_differences_distance(intersect_region):
#     n = len(intersect_region)
#     pre_diff_diff, pre_distance, next_diff_diff, next_distance = 0, 0, 0, 0
#     result = []
#     for inx, content in enumerate(intersect_region):
#         pre_diff_diff = int(content[inx][5]) - int(content[inx-1][5]) if inx > 0 else 0  # 首元素无前驱
#         next_diff_diff = int(content[inx+1][5]) - int(content[inx][5]) if inx < n-1 else 0  # 末元素无后继
#         pre_distance = int(content[inx][1]) - int(content[inx-1][2]) if inx > 0 else 0
#         next_distance = int(content[inx+1][1]) - int(content[inx][2]) if inx < n-1 else 0
        
#         pre_para = diff_dist_para(pre_diff_diff,pre_dist,posdistpara = 0.001,negdistpara = 0.03,const = 0.5)
#         next_para = diff_dist_para(next_diff_diff,next_dist,posdistpara = 0.001,negdistpara = 0.03,const = 0.5)
#         neighbourset = {pre_para,next_para}

#         result.append([n, pre_diff_diff, pre_distance, next_diff_diff, next_distance, min(neighbourset), max(neighbourset)])
#     return result


def demergeinfor(F,samplename,exonmode):
 #chrnameJoin
    allinfor = []
    alluqinfor = []
    if str(exonmode) == "2":
        gms = os.path.join(os.path.dirname(os.path.abspath(__file__)),"gcmappabilityscore_m2.txt")
    else:
        gms = os.path.join(os.path.dirname(os.path.abspath(__file__)),"gcmappabilityscore.txt")
    #print(gms)
    fp5 = open(gms,"r")  #  gc mappabilityscore  size 信息添加到每个外显子后面
    G = {}
    for ggss in fp5:
        ggs = ggss.strip().split("\t")
        key = "_".join([ggs[0],ggs[1],ggs[2]])  # chr_start_end
        G[key] = ggs[5:]

    header = ["libId","chrnameJoin","gene_name","freq","all_freq","nucleotideChangeName","_tag","_result","balance","FOLD80","#chr","start","end","cds","gene","diff","depth","sample_value","flag_outlier","z_score","contrast.mean","contrast.cv","contrast.std","contrast.dps_mean","contrast.flag_bad_capture","tag","GC","mappabilityscore","size", "exon_num", "pre_diff_diff", "pre_dist", "next_diff_diff", "next_dist", "min_diff_dist", "max_diff_dist"]

    #print("\t".join(header))
    f = F
    qc = os.path.join(os.path.dirname(f),"qc_info.txt") #样品均衡性文件
    no = os.path.join(os.path.dirname(f),"exons_details.txt")  #单外显子文件
    fold80file = ""
    if(os.path.exists(os.path.join(os.path.dirname(f)+'/../Stat/',"fold80.txt"))):
        fold80file = os.path.join(os.path.dirname(f)+'/../Stat/',"fold80.txt")
    else:
        fold80file = os.path.join(os.path.dirname(f)+'/',"fold80.txt")
    #fold80file = os.path.join(os.path.dirname(f)+'/',"fold80.txt")
    outfile = os.path.join(os.path.dirname(f),"demergetable.txt")
    with open(fold80file,"r") as f80:
        for fn,f80ls in enumerate(f80):
            f80l = f80ls.strip().split("\t")
            if(fn==0):
                FOLD80i = f80l.index("FOLD80")
                continue
            FOLD80 = f80l[FOLD80i]
                  
    fw = open(outfile,"w")
    fw.write("\t".join(header)+"\n")
    fp2 = open(qc,"r") #记录样品均衡性文件
    #name = f.split("/")[-2].split("_")[0]
    name = samplename
    for ggs in fp2:
        if(re.findall("^#",ggs.strip())):
            continue
        gg = ggs.strip().split("\t")
        bal = re.findall(r'\((.*)\)',gg[-1])[0]  #取出均衡值
    #print(bal,f,no,name)
    
    singlechr = {} #单个外显子文件
    with open(no,"r") as out1:
                for j,lines in enumerate(out1):
                    line = lines.strip().split("\t")
                    if(j==0 or j ==1):
                        continue
                    # if(line[-1]=="NA"):
                    #     continue
                    singlechr.setdefault(line[0],[]).append(line) #键染色体，键值每个单外显子的详细信息
                    # if [*line[:3], line[5]] in singlechr_uq[line[0]]:  # 键染色体, 键值每个单外显子的起始终止和diff值，去重
                    #     continue
                    # singlechr.setdefault(line[0],[]).append([*line[:3], line[5]])
    #format manta
    with open(f,"r") as out: #拆分合并区间成单个exon
        for i,tts in enumerate(out):
            if(i == 0 or i == 1):
                continue
            tt = tts.strip().split("\t")
            if(i == 2):
                chri = tt.index("#chr")
                si = tt.index("start")
                ei = tt.index("end")
                gi = tt.index("gene_name")
                fri = tt.index("freq")
                afri = tt.index("all_freq")
                exoni = tt.index("gene_info_str")
                tagi = tt.index(name +"_tag")
                rei = tt.index(name +"_result")
                continue
            Chr = tt[chri]
            a1 = int(tt[si]) #start
            a2 = int(tt[ei]) #end
            gene = tt[gi]
            fre = tt[fri]
            freall = tt[afri]
            exon = tt[exoni]
            tag = tt[tagi]
            result = tt[rei]
            chrjoin = Chr + ":" + str(a1)  + "-" + str(a2)
            tempcomm = [name,chrjoin,gene,fre,freall,exon,tag,result,bal,FOLD80]  #合并区间文件对每个exon共有的部分
            if Chr not in singlechr:        # 整条染色体都不在exondetails，相当于下面的compare（计算区间是否有交集）也不可能通过，直接跳过就行，此情况初见于AKX186，WESCNVex流程中exondetails额外过滤掉了chrY变异，导致manta这里的chrY变异键值报错 2025.4.11 yzh
                continue
            
            # 取出相交的区间和+-2的区间、计算相交的个数
            intersect_region, exon_num = [], 0
            for idx, line in  enumerate(singlechr[Chr]):
                b1, b2, consisgene, diff = int(line[1]), int(line[2]), line[4], line[5]  #start #end
                size = str(b2-b1+1)
                key1 = "_".join([line[0],line[1],line[2]])
                if(compareinclued(a1,a2,b1,b2) and (consisgene == gene)): #合并区间和单个外显子相交的区间
                    exon_num += 1
                    if idx - 2 >= 0 and singlechr[Chr][idx-2] not in intersect_region:
                        intersect_region.append(singlechr[Chr][idx-2])
                    if idx - 1 >= 0 and singlechr[Chr][idx-1] not in intersect_region:
                        intersect_region.append(singlechr[Chr][idx-1])
                    if line not in intersect_region:
                        intersect_region.append(line)
                    if idx + 1 < len(singlechr[Chr]) and singlechr[Chr][idx+1] not in intersect_region:
                        intersect_region.append(singlechr[Chr][idx+1])
                    if idx + 2 < len(singlechr[Chr]) and singlechr[Chr][idx+2] not in intersect_region:
                        intersect_region.append(singlechr[Chr][idx+2])
            
            # print(intersect_region)
            # 根据相交区间，计算每个区间的前后距离、前后diff值
            for idx, line in  enumerate(intersect_region):
                b1, b2, consisgene, diff = int(line[1]), int(line[2]), line[4], line[5]  #start #end
                size = str(b2-b1+1)
                key1 = "_".join([line[0],line[1],line[2]])

                ## 前一个探针的距离、后一个探针的距离、前一个探针的diff值、后一个探针的diff值
                pre_start, pre_end, pre_diff = b1, b2, diff
                if idx >= 1:
                    pre_start, pre_end, pre_diff = int(intersect_region[idx-1][1]), int(intersect_region[idx-1][2]), intersect_region[idx-1][5]
                    if (pre_start, pre_end) == (b1, b2):
                        if idx > 1:
                            pre_start, pre_end, pre_diff = int(intersect_region[idx-2][1]), int(intersect_region[idx-2][2]), intersect_region[idx-2][5]
                pre_diff_diff = abs(float(diff) - float(pre_diff))
                pre_dist = b1 - pre_end

                next_start, next_end, next_diff = b1, b2, diff
                if idx < len(intersect_region)-1:
                    next_start, next_end, next_diff = int(intersect_region[idx+1][1]), int(intersect_region[idx+1][2]), intersect_region[idx+1][5]
                    if (next_start, next_end) == (b1, b2):
                        if idx < len(intersect_region)-2:
                            next_start, next_end, next_diff  = int(intersect_region[idx+2][1]), int(intersect_region[idx+2][2]), intersect_region[idx+2][5]
                next_diff_diff = abs(float(diff) - float(next_diff))
                next_dist = next_start - b2

                pre_para = diff_dist_para(pre_diff_diff,pre_dist,posdistpara = 0.001,negdistpara = 0.03,const = 0.5)
                next_para = diff_dist_para(next_diff_diff,next_dist,posdistpara = 0.001,negdistpara = 0.03,const = 0.5)
                neighbourset = {pre_para,next_para}
                min_diff_dist = min(neighbourset)
                max_diff_dist = max(neighbourset)

                if(compareinclued(a1,a2,b1,b2) and (consisgene == gene)): #合并区间和单个外显子相交的区间
                    #if revisregion(Chr,b1,b2,G) not in G: print( line , Chr,b1,b2,G, size)
                    if(key1 in G):
                        out = tempcomm + line + G[key1] + [size, exon_num, pre_diff_diff, pre_dist, next_diff_diff,next_dist, min_diff_dist, max_diff_dist]
                    else:
                        out = tempcomm + line + G[revisregion(Chr,b1,b2,G)] + [size, exon_num, pre_diff_diff, pre_dist, next_diff_diff,next_dist, min_diff_dist, max_diff_dist]
                    out = list(map(lambda x: str(x), out))
                    allinfor.append(out)


    for value in allinfor:
        if value not in alluqinfor:
            alluqinfor.append(value)
    for uvalue in alluqinfor:
        fw.write("\t".join(uvalue)+"\n")
    return float(bal)
#demergeinfor(sys.argv[1],sys.argv[2])
