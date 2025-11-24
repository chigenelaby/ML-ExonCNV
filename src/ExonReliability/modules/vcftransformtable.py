#!/share/ofs2a/train/liuyb/0-software/anaconda/bin/python
import sys
import os
import argparse
import vcf
from collections import defaultdict
from collections import Counter
import os
import glob
import xlwt
import xlsxwriter
import re
parser = argparse.ArgumentParser(description="count frequency of variants for samples")
parser.add_argument('-INPUT',help='vcf file',required=True)
parser.add_argument('-OUT',help='output file',required=True)
#parser.add_argument('-MQ',type=float,help='RMS mapping quality {default=50}',default=50)
#parser.add_argument('-DP',type=float,help='Combined depth across samples {default=10}',default=10)
#parser.add_argument('-QUAL',type=float,help='quality {default=10}',default=10)
#parser.add_argument('-zone',help='specific zone for vcf like chr1:0-249250621',required=True)
#parser.add_argument('-output',help='outfile',required=True)
args  = parser.parse_args()


def vcftranstable(path,output):
    # filter some parameters,accoding to requires
    # use pyvcf module to do it
    sample_info = {}
    infile = path
    #for infile in glob.glob(path+"/*.vcf"):
    #        print(infile)
    vcf_reader = vcf.Reader(filename=infile)
    #batch = infile.split("/")[-1].split("_")[0]
    fre_int = []
    fre_str = []
    for record in vcf_reader:
            Chr = record.CHROM
            Pos = str(record.POS)
            Ref = record.REF
            Alt = "|".join(str(i) for i in record.ALT)
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

            for samp in record.samples:
                if(samp['AD'] is None):
                    pass
                else:
                    Sum = sum(samp['AD'])
                    alt = samp['AD'][1:]
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
                    record.INFO['MQ'],ClippingRankSum,samp['GQ'],
                    ReadPosRankSum,MQRankSum,QD,samp["GT"]]
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
        value = [str(i) for i in value]
        genotype = value[-1]
        if(genotype == "1/2"): #["Chr","Pos","参考碱基","突变碱基","突变深度","总深度","突变率","QUAL","FS","SOR","MQ","ClippingRankSum","GQ","ReadPosRankSum","MQRankSum","QD","GT"]
            for base,devar,frevar in zip(value[3].split("|"),value[4].split("|"),value[6].split("|")):
                value[3],value[4],value[6] = base,devar,frevar
                if(float(value[5]) > 20 and float(value[10]) > 55 and (float(value[6]) > 0.2 and float(value[6]) < 0.8)): #总深度，MQ，突变率
                    variantsdict.setdefault(value[0],[]).append(value[1])
                    #print(value)
                    #pass
        elif(genotype =="0/1"):
            if(float(value[5]) > 20 and float(value[10]) > 55 and (float(value[6]) > 0.2 and float(value[6]) < 0.8)):
                #print(value)
                variantsdict.setdefault(value[0],[]).append(value[1])
    #t = 0
    #for c,ccs in variantsdict.items():
    #    t += len(ccs)
    #print(t)
    #for k,vs in variantsdict.items():
    #    print(k,vs)



    return  variantsdict
filter_vcf(args.INPUT,args.OUT)

