import sys
from collections import Counter
from collections import defaultdict
import os
#c = defaultdict(dict)

#合并区间时发现，相对可靠单外显子错误的比较多，然后进一步划分为不可靠
def unre(freq,balance,contrast_mean,contrast_std,diff,all_freq,GC,sample_value,z_score):
    '''
    unre1 = freq > 0.0055 and balance > 0.00735 and contrast_mean <= 4.7773 # 141 neg  13 pos
    unre2 = freq <= 0.0055 and contrast_std <= 0.19735 and diff > 0.64165 and diff <= 1.378  #21 neg,0 pos
    unre3 = freq <= 0.0055 and contrast_std <= 0.19735 and diff <= 0.64165 and all_freq > 0.0045 and GC > 0.409213  # 16 neg ,1 pos
    unre4 = freq > 0.0055 and balance <=  0.00735 and all_freq > 0.0165 and sample_value > 0.33355 # 12 neg  0 pos
    unre5 = freq <= 0.0055 and contrast_std > 0.19735 and GC > 0.582045 and balance > 0.01445 # 5 neg 0 pos
    unre6 = freq <= 0.0055 and contrast_std <= 0.19735 and diff > 0.64165 and diff > 1.378 and balance > 0.0176 # 7 neg 0 pos
    unre7 = freq <= 0.0055 and contrast_std <= 0.19735 and diff <= 0.64165 and all_freq <= 0.0045 and diff > 0.5709 and GC <= 0.540771 # 5 neg  0 pos
    '''
    unre1 = z_score > -3.5044 and freq > 0.0005 and diff > 0.6195 and sample_value <= 5.7547
    if(any([unre1])):
        return 1
    else:
        return 0
def finetuneregion(F,orif,samplename):
    ve1 = "特别可靠"
    re1 = "相对可靠"
    #un1 = "可能不可靠"
    un2 = "不可靠"
    #unk = "未知"
    balancedict = {}
    a = {}  #存每个区间有多少单个外显子
    g = {}  #记录一下样品，位置，判断的条件，基因
    outfile = os.path.join(os.path.dirname(F),"regionexonfinetune.txt")
    fw = open(outfile,"w")
    with open(F) as out:
        for i,lines in enumerate(out):
            line = lines.strip().split("\t")
            if(i == 0):
                libIdi = line.index("libId")
                chrnameJoini = line.index("chrnameJoin")
                resulti = line.index("_result")
                tagi = line.index("_tag")
                #flagi = line.index("flag")
                freqi = line.index("freq")
                gi = line.index("gene_name")
                #sangerCarryi = line.index("sangerCarry")
                #ngsCarryResulti = line.index("ngsCarryResult")
                #taskCodei = line.index("taskCode")
                balancei = line.index("balance")
                contrast_meani = line.index("contrast.mean")
                contrast_stdi = line.index("contrast.std")
                GCi = line.index("GC")
                all_freqi = line.index("all_freq")
                diffi = line.index("diff")
                yuani = line.index("原始可靠性标注")
                continue
            libId = line[libIdi]
            chrnameJoin = line[chrnameJoini]
            result = line[resulti]
            #flag = line[flagi]
            yuan = line[yuani]
            freq = float(line[freqi])
            all_freq = float(line[all_freqi])
            balance = float(line[balancei])
            diff = float(line[diffi])
            contrast_mean = float(line[contrast_meani])
            contrast_std = float(line[contrast_stdi])
            GC = float(line[GCi])
            gene = line[gi]
            key = libId +"/" + chrnameJoin + "/" + gene
            #a.setdefault(key,[]).append([result,flag])
            #print(key)
            balancedict[key] = balance
            a.setdefault(key,[]).append(result)
            #g.setdefault(key,set()).add(yuan)
            g[key] = [line[libIdi],line[chrnameJoini],line[tagi],gene]
    uniq = []
    for key,values in a.items():
        typecount = dict(Counter(values))
        if(re1 in typecount and typecount[re1] ==1):
            uniq.append(key)
    mdict = defaultdict(dict)
    #单外显子判断条件存储
    with open(F) as out:
        for i,lines in enumerate(out):
            line = lines.strip().split("\t")
            if(i == 0):
                libIdi = line.index("libId")
                chrnameJoini = line.index("chrnameJoin")
                resulti = line.index("_result")
                #flagi = line.index("flag")
                freqi = line.index("freq")
                gi = line.index("gene_name")
                balancei = line.index("balance")
                contrast_meani = line.index("contrast.mean")
                contrast_stdi = line.index("contrast.std")
                GCi = line.index("GC")
                all_freqi = line.index("all_freq")
                diffi = line.index("diff")
                z_scorei = line.index("z_score")
                sample_valuei = line.index("sample_value")
                yuani = line.index("原始可靠性标注")
                continue
            libId = line[libIdi]
            chrnameJoin = line[chrnameJoini]
            result = line[resulti]
            #flag = line[flagi]
            yuan = line[yuani]
            gene = line[gi]
            balance = float(line[balancei])
            key = libId +"/" + chrnameJoin + "/" + gene
            if(key in uniq):
                mdict[key]["freq"] = float(line[freqi])
                mdict[key]["all_freq"] = float(line[all_freqi])
                mdict[key]["balance"] = float(line[balancei])
                mdict[key]["diff"] = float(line[diffi])
                mdict[key]["contrast_mean"] = float(line[contrast_meani])
                mdict[key]["contrast_std"] = float(line[contrast_stdi])
                mdict[key]["GC"] = float(line[GCi])
                mdict[key]["sample_value"] = float(line[sample_valuei])
                mdict[key]["z_score"] = float(line[z_scorei])
            #a.setdefault(key,[]).append([result,flag]
    judgere = []
    storecon = []
    storeconre = []
    #print("libId","chrnameJoin","ngsCarryResult","whetherreliabilty","basedconditon",sep="\t")
    #根据目标调整合并区间规则
    for key,values in a.items():
        #print(key,Counter(values))
        typecount = dict(Counter(values))
        #id = "\t".join(key.split("/"))
        #print(key,values)
        #print(key,typecount)
        posgene = ":".join(key.split("/")[1:])
        #特别可靠
        if(ve1 in typecount and len(typecount)==1):
            #print(id,ve1,typecount[ve1],sep="\t")
            judgere.append([posgene,ve1]) #特别可靠
            storecon.append(g[key] + [ve1] + ['1'])
            #pass
        elif((ve1 in typecount and re1 in typecount) and len(typecount)==2):
            if(typecount[ve1] >= 2):
                judgere.append([posgene,ve1]) #特别可靠
                storecon.append(g[key] + [ve1] + ['2'])
                #pass
            else:
                judgere.append([posgene,re1])  #相对可靠
                storecon.append(g[key] + [re1] + ['8'])
                storeconre.append(g[key] + [re1] + ['8']+[ve1,typecount[ve1],re1,typecount[re1],un2,"0"])
                #judgere.append([posgene,ve1])
                #print(id,ve1,typecount[ve1],re1,typecount[re1],sep="\t")
        elif((ve1 in typecount and un2 in typecount) and len(typecount)==2):
            if(typecount[ve1]/(typecount[un2]+typecount[ve1]) >= 0.9):
                judgere.append([posgene,ve1]) #特别可靠
                storecon.append(g[key] + [ve1] + ['3'])
            elif(typecount[ve1]/(typecount[un2]+typecount[ve1]) < 0.9 and typecount[ve1]/(typecount[un2]+typecount[ve1]) >= 0.5):
                judgere.append([posgene,re1]) #相对可靠
                storecon.append(g[key] + [re1] + ['14'])
                storeconre.append(g[key] + [re1] + ['14']+[ve1,typecount[ve1],re1,"0",un2,typecount[un2]])
            else:
                judgere.append([posgene,un2]) #不可靠
                storecon.append(g[key] + [un2] + ['5'])
        #else:
        #    pass
        #不可靠
        elif(un2 in typecount and len(typecount)==1):
            #print(id,un2,typecount[un2],sep="\t")
            judgere.append([posgene,un2])
            storecon.append(g[key] + [un2] + ['4'])
            #pass
        elif((un2 in typecount and re1 in typecount) and len(typecount) ==2):
            #if(typecount[un2] >=2 and int(typecount[un2]) > int(typecount[re1])):
            #print(key,typecount)
            if(typecount[un2]/(typecount[un2]+typecount[re1]) > 0.3 and balancedict[key] > 0.01):
                judgere.append([posgene,un2])
                storecon.append(g[key] + [un2] + ['6'])
                #print(key,typecount)
            elif(typecount[un2] >=2 and int(typecount[un2]) > int(typecount[re1])):
            #elif(typecount[un2]/(typecount[un2]+typecount[re1]) > 0.3 and balance > 0.01):
                judgere.append([posgene,un2])
                #print(id,un2,typecount[un2],re1,typecount[re1])
                storecon.append(g[key] + [un2] + ['6'])
                #print(key,typecount)
                #pass
            else:
                 judgere.append([posgene,re1]) #相对可靠
                 storecon.append(g[key] + [re1] + ['9'])
                 storeconre.append(g[key] + [re1] + ['9']+[ve1,"0",re1,typecount[re1],un2,typecount[un2]])
         #       pass
        #单外显子注释掉
        elif((re1 in typecount and typecount[re1] == 1)):
               judgere.append([posgene,re1])  #相对可靠
               storecon.append(g[key] + [re1] + ['10'])
               storeconre.append(g[key] + [re1] + ['10']+[re1,typecount[re1],ve1,"0",un2,"0"])
        #相对可靠
        elif((re1 in typecount and typecount[re1] > 1 and len(typecount) == 1)):
            judgere.append([posgene,re1])
            storecon.append(g[key] + [re1] + ['11'])
            storeconre.append(g[key] + [re1] + ['11']+[ve1,"0",re1,typecount[re1],un2,"0"])
        elif((ve1 in typecount and re1 in typecount and un2 in typecount) and len(typecount)==3):
            if((typecount[un2]/(typecount[ve1] + typecount[re1] +typecount[un2]) <= 0.50)):
                if(typecount[ve1] > typecount[re1]):
                    judgere.append([posgene,ve1])
                    storecon.append(g[key] + [ve1] + ['12'])
                else:
                    judgere.append([posgene,re1])
                    storecon.append(g[key] + [re1] + ['12'])
                    storeconre.append(g[key] + [re1] + ['12']+[ve1,typecount[ve1],re1,typecount[re1],un2,typecount[un2]])
            else:
                judgere.append([posgene,un2])
                storecon.append(g[key] + [un2] + ['13'])

            #print(key,values)
    #print(len(judgere))
    '''
    ldict = {}
    for rf in judgere:
        #print(rf)
        ys,key = rf
        ldict.setdefault(key,[]).append(ys)
    j = 0
    for key,values in ldict.items():
        #print(key,Counter(values))
        for value in dict(Counter(values)).values():
            j += value
    print(j)
    '''
    corlabel = {}
    for k in storecon:
        #pass
        samplename,pos,tagw,gene,wreun,condi = k
        key = pos+":"+gene
        value = wreun
        corlabel[key] = value
        print(samplename,pos,tagw,gene,wreun,sep="\t")
    #for k in storeconre:
    #    print("\t".join([str(i) for i in k]))
    #原始文件format_gene_info.txt
    '''
    name = samplename
    with open(orif,"r") as fgi:
        for i,oris in enumerate(fgi):
            if(i == 0 or i == 1):
                fw.write(oris)
                continue
            ori = oris.strip().split("\t")
            if(i == 2): ##chr    start   end     gene_name       gene_info_str   best_exon_str   freq    all_freq        infos
                chri = ori.index("#chr")
                si = ori.index("start")
                ei = ori.index("end")
                gi = ori.index("gene_name")
                gisi = ori.index("gene_info_str")
                besi  = ori.index("best_exon_str")
                fri = ori.index("freq")
                afri = ori.index("all_freq")
                infosi = ori.index("infos")
                exoni = ori.index(name +"_exon")
                tagi = ori.index(name +"_tag")
                rei = ori.index(name +"_result")
                header = "\t".join(["#chr","start","end","gene_name","gene_info_str","best_exon_str","freq","all_freq","infos",name +"_exon",name +"_tag",name +"_result"])
                fw.write(header +"\n")
                continue
            keyori = ori[chri] + ":" + ori[si] + "-" + ori[ei] + ":" + ori[gi]
            if(keyori in corlabel):
                ori[rei] = corlabel[keyori]
                fw.write("\t".join(ori[:12]) + "\n")
            else:
                fw.write("\t".join(ori[:12]) + "\n")



     '''
finetuneregion(sys.argv[1],sys.argv[2],sys.argv[3])
