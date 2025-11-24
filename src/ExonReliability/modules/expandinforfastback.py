import sys
import os
import re
#合并的区间拆分为单个外显子以其特征如balance,gc,z_score等
def compare(a1,a2,b1,b2): #比较区间是否有交集
  Max  = max([a1,b1])
  Min  = min([a2,b2])
  if(Min - Max > 0):
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
def demergeinfor(F,samplename):
 #chrnameJoin
    allinfor = []
    alluqinfor = []
    gms = os.path.join(os.path.dirname(os.path.abspath(__file__)),"gcmappabilityscore.txt")
    #print(gms)
    fp5 = open(gms,"r")  #  gc mappabilityscore  size 信息添加到每个外显子后面
    G = {}
    for ggss in fp5:
        ggs = ggss.strip().split("\t")
        key = "_".join([ggs[0],ggs[1],ggs[2]])  # chr_start_end
        G[key] = ggs[5:]

    header = ["libId","chrnameJoin","gene_name","freq","all_freq","nucleotideChangeName","_tag","_result","balance","#chr","start","end","cds","gene","diff","depth","sample_value","flag_outlier","z_score","contrast.mean","contrast.cv","contrast.std","contrast.dps_mean","contrast.flag_bad_capture","tag","GC","mappabilityscore","size"]

    #print("\t".join(header))
    f = F
    qc = os.path.join(os.path.dirname(f),"qc_info.txt") #样品均衡性文件
    no = os.path.join(os.path.dirname(f),"exons_details.txt")  #单外显子文件
    outfile = os.path.join(os.path.dirname(f),"demergetable.txt")
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
                        #flag_outlieri = line.index("flag_outlier")
                        #contrast_flag_bad_capturei = line.index("contrast.flag_bad_capture")
                        continue
                    #line = lines.strip().split("\t")
                    if(line[-1]=="NA"):
                        continue
                    #Chr1 = line[0]
                    #b1 = int(line[1]) #start
                    #b2 = int(line[2]) #end
                    #consisgene = line[4]
                    #size = str(b2-b1+1)
                    #key1 = "_".join([line[0],line[3],line[4]])
                    singlechr.setdefault(line[0],[]).append(line) #键染色体，键值每个单外显子的详细信息

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
            tempcomm = [name,chrjoin,gene,fre,freall,exon,tag,result,bal]  #合并区间文件对每个exon共有的部分
            for line in singlechr[Chr]:
                #Chr1 = line[0]
                b1 = int(line[1]) #start
                b2 = int(line[2]) #end
                consisgene = line[4]
                size = str(b2-b1+1)
                key1 = "_".join([line[0],line[1],line[2]])
                if(compare(a1,a2,b1,b2) and (consisgene == gene)): #合并区间和单个外显子相交的区间
                    if(key1 in G):
                        allinfor.append(tempcomm + line + G[key1] + [size])
                    else:
                        allinfor.append(tempcomm + line + G[revisregion(Chr,b1,b2,G)] + [size])

    for value in allinfor:
        if value not in alluqinfor:
            alluqinfor.append(value)
    for uvalue in alluqinfor:
        fw.write("\t".join(uvalue)+"\n")
    return float(bal)
demergeinfor(sys.argv[1],sys.argv[2])
