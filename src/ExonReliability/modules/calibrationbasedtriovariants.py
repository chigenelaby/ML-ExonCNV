import os
import sys
import re
import vcfpy
import json

def vcftranstable(mvcffile):
    # filter some parameters,accoding to requires
    # use pyvcf module to do it
    sample_info = {}
    infile = mvcffile
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
            if("MQ" not in record.INFO.keys()):
                continue
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

            for samp in record.calls:
                if(samp['AD'] is None):
                    pass
                else:
                    #Sum = sum(samp['AD'])
                    Sum = sum(samp.data.get('DP'))
                    alt = samp.data.get('AD')[1:]
                    depth_list = [str(i) for i in alt]
                    variant_depth = "|".join(depth_list)
                    if(Sum == 0):
                        fre = 0
                    else:
                        fre_str = [str(round(i/Sum,3)) for i in alt]
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
    variantsdict = {}
    for key,value in sample_info.items():
        #print(key,value)
        sample = key.split("|")[0]
        #value = [str(i) for i in value]
        #print(sample)
        if(re.findall(r"\|",value[3])):
            continue
        if(float(value[5]) > 20 and float(value[10]) > 55 and (value[-1] == "0/0" or value[-1] == "1/1")): #总深度,MQ,genotype
        #    print("done")
            if(value[-1] == "0/0"):
                variantsdict.setdefault(sample,[]).append(value[:2] + ["wild"])
            else:
                variantsdict.setdefault(sample,[]).append(value[:2] + ["homo"])
        if(float(value[5]) > 20 and float(value[10]) > 55 and (value[-1] == "0/1") and (float(value[6]) > 0.2 and float(value[6]) < 0.8)):
            variantsdict.setdefault(sample,[]).append(value[:2] + ["heter"])

    return variantsdict
def overlap(provar,Chr,a1,a2):
    temp = set()
    for pvars in provar[Chr]:
                ppos = int(pvars[0])
                pflag = pvars[1]
                if(ppos > a1 and ppos < a2):
                    temp.add(pflag)
    return temp

def calivarianttrio(a,samplenamelist,variantfile):
    regionfile = os.path.join(os.path.dirname(a[0]),"format_gene_inforevis.txt")
    variants = {}
    variants = vcftranstable(variantfile)
    condition = [["wild","homo","wild"],["wild","wild","homo"],["homo","homo","wild"],["homo","wild","homo"]]
    conditionm = [["wild","homo","heter"],["wild","heter","homo"],["homo","heter","wild"],["homo","wild","heter"]] #modification 2024/02/23
    outfile = os.path.join(os.path.dirname(regionfile),"format_gene_inforevistrio.txt")
    fw = open(outfile,"w")
    provar = {} #先证者
    for pv in variants[samplenamelist[0]]:
        provar.setdefault(pv[0],[]).append(pv[1:])
    fmvar = {} #父或母
    for fv in variants[samplenamelist[1]]:
        fmvar.setdefault(fv[0],[]).append(fv[1:])
    mfvar = {} #父或母
    for mv in variants[samplenamelist[2]]:
        mfvar.setdefault(mv[0],[]).append(mv[1:])

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
                infosi = reline.index("infos")
                rei = reline.index(samplenamelist[0]+"_result")
                flagi = reline.index(samplenamelist[0]+"_tag")
                continue
            Chr = reline[chri]
            a1 = int(reline[si])
            a2 = int(reline[ei])
            wereli = reline[rei]
            #info = reline[infosi]
            flag = reline[flagi]
            if(flag != "loss1"):
                fw.write("\t".join(reline) + "\n")
                continue
            if(wereli == "特别可靠"):
                fw.write("\t".join(reline) + "\n")
                continue
            #k = 0
            ptemp = overlap(provar,Chr,a1,a2)
            ftemp = overlap(fmvar,Chr,a1,a2)
            mtemp = overlap(mfvar,Chr,a1,a2)
            if(len(ptemp) !=0 and len(ftemp) !=0 and len(mtemp) !=0):
                num = len(ptemp) + len(ftemp) + len(mtemp)
                if(len(ptemp)==len(ftemp)==len(mtemp)==1): #
                    combin = list(ptemp) + list(ftemp) + list(mtemp)
                    if((combin in condition or combin in conditionm) and wereli=="相对可靠"):
                        vafsupport = json.loads(reline[infosi])
                        vafsupport['reliability_show'].append("家系突变支持")
                        reline[infosi] = json.dumps(vafsupport,ensure_ascii=False)
                        reline[rei] = "特别可靠"
                    if((combin in condition or combin in conditionm) and wereli=="不可靠"):
                        vafsupport = json.loads(reline[infosi])
                        vafsupport['reliability_show'].append("家系突变支持")
                        reline[infosi] = json.dumps(vafsupport,ensure_ascii=False)
                        reline[rei] = "相对可靠"
                if(("heter" not in ptemp and "heter" not in ftemp and "heter" not in mtemp) and num == 4 and ptemp != ftemp and ptemp != mtemp and ftemp != mtemp): #{'homo'} {'wild', 'homo'} {'wild'}针对这种情况
                    if(wereli=="相对可靠"):
                        vafsupport = json.loads(reline[infosi])
                        vafsupport['reliability_show'].append("家系突变支持")
                        reline[infosi] = json.dumps(vafsupport,ensure_ascii=False)
                        reline[rei] = "特别可靠"
                    if(wereli=="不可靠"):
                        vafsupport = json.loads(reline[infosi])
                        vafsupport['reliability_show'].append("家系突变支持")
                        reline[infosi] = json.dumps(vafsupport,ensure_ascii=False)
                        reline[rei] = "相对可靠"
            fw.write("\t".join(reline) + "\n")



#vcftranstable(sys.argv[1])
#a= "/share/gwasfs3a/shcai/file/exonlossgc/scriptpack/finetuneexonlossgain/datatrio4members/DDN22001823/AFD922_exon/format_gene_inforevis.txt"
#b = ["AFD922","AFD924","AFD925","AFD923"]
#c = "/share/chg1fs1b/prod/project/DDN22001000-22001999/DDN22001823/WES_VCF/fam.vcf.gz"
#c = "/share/gwasfs3a/shcai/file/exonlossgc/scriptpack/finetuneexonlossgain/modules/fam.vcf.gz"
#calivarianttrio(a,b,c)

