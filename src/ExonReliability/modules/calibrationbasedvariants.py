import os
import sys
import re
import vcfpy
import json
import bisect
#
def coverageregion():
        goodcovregion = {}
        cov = os.path.join(os.path.dirname(os.path.abspath(__file__)),"coverageregion1000sample.txt")
        with open(cov,"r") as out:
                for i,lines in enumerate(out):
                        line = lines.strip().split("\t")
                        if(i==0):
                                Posi = line.index("Pos")
                                depthXi = line.index("30X")
                                nsamplei = line.index("total")
                                continue
                        depthX = int(line[depthXi])
                        Chr,start,end = line[Posi].split("_")
                        key = "_".join([Chr,start,end])
                        nsample = int(line[nsamplei])
                        if(depthX/nsample > 0.99):
                            goodcovregion.setdefault(Chr,[]).append([int(start),int(end),round(depthX/nsample,3)])
                            #print(depthX,nsample,lines.strip())
        return goodcovregion

def regionoverlap(goodregion,Chr,a1,a2):
                total = 0
                for poslist in goodregion[Chr]:
                        b1,b2,score = poslist
                        #print(Chr,a1,a2,b1,b2)
                        if((min(a2,b2) - max(a1,b1)) > 0):
                                total += 1
                                break
                if(total == 0):
                    return 0
                if(total == 1):
                    return 1

def ssrfilter():
        region = os.path.join(os.path.dirname(os.path.abspath(__file__)),"ssr.bed")
        ssrregion = {}
        ssrregionsort = {}
        ssrregionfirstelemnt = {}
        ssrregionsecondelemnt = {}
        with open(region,"r") as out:
                for lines in out:
                        line = lines.strip().split("\t")
                        Chr = line[0]
                        Start = int(line[1])
                        End = int(line[2])
                        ssrregion.setdefault(Chr,[]).append([Start,End])
        for key,values in ssrregion.items():
                ssrregionsort[key] = sorted(values, key=lambda x: x[0])
        for key,sten in ssrregionsort.items():
                ssrregionfirstelemnt[key] = [r[0] for r in sten]
                ssrregionsecondelemnt[key] = [r[1] for r in sten]
        return ssrregionfirstelemnt,ssrregionsecondelemnt



def vcftranstable(vcffile):
    # filter some parameters,accoding to requires
    # use pyvcf module to do it
    sample_info = {}
    infile = vcffile
    #for infile in glob.glob(path+"/*.vcf"):
    #        print(infile)
    vcf_reader = vcfpy.Reader.from_path(infile)
    #batch = infile.split("/")[-1].split("_")[0]
    fre_int = []
    fre_str = []
    for record in vcf_reader:
            Chr = record.CHROM
            Pos = str(record.POS)
            Ref = record.REF
            Alt = "|".join(str(alt.value) for alt in record.ALT)
            Qual = record.QUAL
            if("MQRankSum" not in record.INFO.keys()):
                MQRankSum = "-"
            else:
                MQRankSum = record.INFO['MQRankSum']
            if("ReadPosRankSum" not in record.INFO.keys()):
                ReadPosRankSum = "-"
            else:
                ReadPosRankSum = record.INFO["ReadPosRankSum"]
            if("ClippingRankSum" not in record.INFO.keys()):
                ClippingRankSum = "-"
            else:
                ClippingRankSum = record.INFO["ClippingRankSum"]
            if("QD" not in record.INFO.keys()):
                QD = "-"
            else:
                QD = record.INFO["QD"]
            if("FS" not in record.INFO.keys()):
                FS = "-"
            else:
                FS = record.INFO["FS"]
            if("FS" not in record.INFO.keys()):
                SOR = "-"
            else:
                SOR = record.INFO["SOR"]

            if "MQ" in record.INFO.keys() and record.INFO['MQ']<55:
                continue
            
            for samp in record.calls:
                
                if(samp.data.get('AD') is None):
                    pass
                else:
                #    Sum = sum(samp['AD'])
                    Sum = samp.data.get('DP')
                    alt = samp.data.get('AD')[1:]
                    depth_list = [str(i) for i in alt]
                    variant_depth = "|".join(depth_list)
                    if(Sum == 0):
                        fre = 0
                    else:
                        fre_str = [str(round(i/Sum,3)) for i in alt ]
                        fre = "|".join(fre_str)
                    key = "|".join([samp.sample,Chr,Pos])

                    sample_info[key] = [Chr,Pos,Ref,Alt,variant_depth,
                    Sum,fre,Qual,FS,SOR,
                    record.INFO['MQ'],ClippingRankSum,samp.data.get('GQ'),
                    ReadPosRankSum,MQRankSum,QD,samp.data.get("GT")]
    #for key,value in sample_info.items():
    #    print(key,value)
    #return sample_info
    #fo = open(output,"w")
    #header =["sample","Chr","Pos","参考碱基","突变碱基","突变深度","总深度","突变率","QUAL","FS","SOR","MQ","ClippingRankSum","GQ","ReadPosRankSum","MQRankSum","QD","GT"]
    #header = "\t".join(header)
    #fo.write(header + "\n")
    ssrregion = {}
    ssrregionfirstelemnt,ssrregionsecondelemnt  = ssrfilter()
    variantsdict = {}
    for key,value in sample_info.items():
        #print(key,value)
        sample = key.split("|")[0]
        #ssrnumber = 0
        #for msted in ssrregion[value[0]]:
        #    ssrst,ssred = msted
        #    if(int(value[1])>=ssrst and int(value[1])<=ssred):
        #        ssrnumber += 1
        #if(ssrnumber !=0):
        #    continue
        value = [str(i) for i in value]
	#for msted in ssrregion[value[0]]:
        genotype = value[-1]
        hm1= [float(value[5]) > 20 and float(value[10]) < 58 and float(frevar) > 0.35 and float(frevar) < 0.65 for base,devar,frevar in zip(value[3].split("|"),value[4].split("|"),value[6].split("|"))]
        hm2= [float(value[5]) > 20 and float(value[10]) >= 58 and float(frevar) > 0.2 and float(frevar) < 0.8 for base,devar,frevar in zip(value[3].split("|"),value[4].split("|"),value[6].split("|"))]
        if(genotype == "1/2"): #["Chr","Pos","参考碱基","突变碱基","突变深度","总深度","突变率","QUAL","FS","SOR","MQ","ClippingRankSum","GQ","ReadPosRankSum","MQRankSum","QD","GT"]
            for base,devar,frevar in zip(value[3].split("|"),value[4].split("|"),value[6].split("|")):
                value[3],value[4],value[6] = base,devar,frevar
                #ssrnumber1 = 0 
                if any(hm1 + hm2): #总深度，MQ，突变率
                    #if(value[0] == "chr17"):
                    #    for msted in ssrregion[value[0]]:
                    #        ssrst,ssred = msted
                    #        if(int(value[1])>=ssrst and int(value[1])<=ssred):
                    #            ssrnumber1 += 1   
                    #if(ssrnumber1 == 0):
                    #ssrkey = value[0]+"_"+value[1]
                    Chr = value[0]
                    pos = int(value[1])
                    insert_index = bisect.bisect_left(ssrregionfirstelemnt[Chr],pos)
                    if(insert_index > 0 and ssrregionfirstelemnt[Chr][insert_index - 1] <= pos <= ssrregionsecondelemnt[Chr][insert_index - 1]):
                        pass
                    else:
                    #if(ssrkey not in ssrregion):
                        variantsdict.setdefault(value[0],[]).append(int(value[1]))
                    #print(value)
                    #pass
        elif(genotype =="0/1"):
            #ssrnumber2 = 0
            if any(hm1 + hm2):
                #if(value[0] == "chr17"):
                #    for msted in ssrregion[value[0]]:
                #        ssrst,ssred = msted
                #        if(int(value[1])>=ssrst and int(value[1])<=ssred):
                #            ssrnumber2 += 1
                #print(value)
                #if(ssrnumber2 == 0):
                #ssrkey = value[0]+"_"+value[1]
                Chr = value[0]
                pos = int(value[1])
                insert_index = bisect.bisect_left(ssrregionfirstelemnt[Chr],pos)
                if(insert_index > 0 and ssrregionfirstelemnt[Chr][insert_index - 1] <= pos <= ssrregionsecondelemnt[Chr][insert_index - 1]):
                    pass
                #if(ssrkey not in ssrregion):
                else:
                    variantsdict.setdefault(value[0],[]).append(int(value[1]))
    #t = 0
    #for c,ccs in variantsdict.items():
    #    t += len(ccs)
    #print(t)
    #for k,vs in variantsdict.items():
    #    print(k,vs)



    return  variantsdict
'''
def compare(a1,a2,b1,b2): #比较区间是否有交集
  Max  = max([a1,b1])
  Min  = min([a2,b2])
  if(Min - Max > 0):
    return 1
  else:
    return 0
'''
def calivariant(regionfile,samplename,variantfile):
    goodregion = coverageregion()
    variants = {}
    variants = vcftranstable(variantfile)
    outfile = os.path.join(os.path.dirname(regionfile),"regionexonfinetuneaddvariant.txt")
    #["gain","gain1","gain2","loss1","loss2","loss"]
    #2024/10/25
    homoloss = ["loss2","loss"]  #打标签，测序稳定区域，1000个样本中>99%样本平均深度>30X
    #discardtype = ["gain","gain1"]  #不做处理
    specialdealgain2 = ["gain2","gain"] # 特别可靠保留 ；只要有纯合(双倍)、点突变支持、断点 其中之一的1个或<300bp  保留为特别可靠，其他1个或<300bp的最多为相对可靠
    fw = open(outfile,"w")
    with open(regionfile,"r") as refile:
        for relines in refile:
            reline = relines.strip().split("\t")
            if(re.findall("^##",relines.strip())):
                fw.write(relines)
                continue
            if(re.findall("^#",relines.strip())):
                fw.write(relines)
                chri = reline.index("#chr")
                si = reline.index("start")
                ei = reline.index("end")
                freqi = reline.index("freq")
                infosi = reline.index("infos")
                rei = reline.index(samplename+"_result")
                flagi = reline.index(samplename+"_tag")
                gisi = reline.index("gene_info_str")
                continue
            Chr = reline[chri]
            a1 = int(reline[si])
            a2 = int(reline[ei])
            wereli = reline[rei]
            flag = reline[flagi]
            infos = reline[infosi]
            try:
                freq = float(reline[freqi])
            except:
                print(reline)
                freq = 0
            gis = reline[gisi]
            regionlegth = a2 - a1 + 1
            if(wereli=="特别可靠" and freq > 0.05):
                reline[rei] = "相对可靠"
            if(re.findall("存在至少一个家系成员支持",infos) and re.findall("wgs支持",infos) and wereli=="不可靠"):
                reline[rei] = "相对可靠"
            #if(flag in discardtype):
            #    fw.write("\t".join(reline) + "\n")
            #    continue
            if(wereli == "不可靠"):
                fw.write("\t".join(reline) + "\n")
                continue
            k = 0
            temp = []
            if(flag == "loss1"):
                for vars in variants[Chr]:
                    if(vars > a1 and vars < a2):
                        temp.append([Chr,vars])
                        k += 1
                if(k != 0):
                    if(wereli=="特别可靠"):
                        reline[rei] = "相对可靠"
                    if(wereli=="相对可靠"):
                        reline[rei] = "不可靠"
            if(flag in homoloss):
                #print(reline)
                if(regionoverlap(goodregion,Chr,a1,a2)):
                    seqstabregion = json.loads(reline[infosi])
                    seqstabregion['reliability_show'].append("测序稳定区域")
                    reline[infosi] = json.dumps(seqstabregion,ensure_ascii=False)
            if((not re.findall("-",reline[gisi]) and not re.findall("all",reline[gisi])) or regionlegth < 300) and wereli=="特别可靠" and flag not in specialdealgain2:
                #print(reline[infosi],reline)
                #print("##",re.findall("测序稳定区域",reline[infosi]),re.findall("断点外显子支持",reline[infosi]))
                if(re.findall("测序稳定区域",reline[infosi])): 
                    pass
                elif(re.findall("断点外显子支持",reline[infosi])):
                    pass
                    #print(reline[infosi],reline)
                else:
                    reline[rei] = "相对可靠"
            fw.write("\t".join(reline) + "\n")
            #if(len(temp) !=0):
            #    print(reline,wereli,temp)
            #print(temp)

#calivariant(sys.argv[1],sys.argv[2],sys.argv[3])



