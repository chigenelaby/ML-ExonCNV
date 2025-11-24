import sys
import os
import re

def compare(a1,a2,b1,b2): #比较区间是否有交集
  Max  = max([a1,b1])
  Min  = min([a2,b2])
  if(Min - Max + 1 > 0):
    return (Min - Max + 1)
  else:
    return 0

def trio(a,b):
    k = 0
    pro = {} #先证者
    proinfo = {} #先证者，修改info信息
    other = {} #其他成员
    #print(a,b)
    outfile = os.path.join(os.path.dirname(a[0]),"format_gene_inforevis.txt")
    fw = open(outfile,"w")
    for i,j in zip(a,b):#返回文件路径，样品名
        #print(i,j)
        path = os.path.dirname(i)
        k += 1
        if(k == 1): #先证者finetune之后的单样品文件
            with open(os.path.join(path,"regionexonfinetuneaddvariant.txt"),"r") as prof:
                for prolines in prof:
                    if(re.findall("^##",prolines.strip())):
                        continue
                    proline = prolines.strip().split("\t")
                    if(re.findall("^#",prolines.strip())): ##chr    start   end     gene_name
                        pchri = proline.index("#chr")
                        psi = proline.index("start")
                        pei = proline.index("end")
                        pgenei = proline.index("gene_name")
                        prei = proline.index(j+"_result")
                        pinfosi = proline.index("infos")
                        continue
                    pchr = proline[pchri]
                    ps = proline[psi]
                    pe = proline[pei]
                    pgene = proline[pgenei]
                    pre = proline[prei]
                    pkey = "|".join([pchr,ps,pe,pgene])
                    pro[pkey] = pre
                    proinfo[pkey] = proline[pinfosi]
        else: #其他家庭成员finetune之后的单样品文件
            with open(os.path.join(path,"regionexonfinetuneaddvariant.txt"),"r") as otherf:
                for otherlines in otherf:
                    if(re.findall("^##",otherlines.strip())):
                        continue
                    otherline = otherlines.strip().split("\t")
                    if(re.findall("^#",otherlines.strip())): ##chr    start   end     gene_name
                        ochri = otherline.index("#chr")
                        osi = otherline.index("start")
                        oei = otherline.index("end")
                        ogenei = otherline.index("gene_name")
                        orei = otherline.index(j+"_result")
                        continue
                    otchr = otherline[ochri]
                    ots = otherline[osi]
                    ote = otherline[oei]
                    otgene = otherline[ogenei]
                    otre = otherline[orei]
                    otvalue = [otchr,ots,ote,otgene,otre]
                    other.setdefault(j,[]).append(otvalue)
    orfamilyfile = a[0] #原始的先证者家系文件
    orff = open(orfamilyfile,"r")
    recoderindex = {} #记录家系成员在先证者家系文件那一列
    for orlines in orff:
        orline = orlines.strip().split("\t")
        if(re.findall("^##",orlines.strip())):
            #print(orlines.strip())
            fw.write(orlines)
            continue
        #orline = orlines.strip().split("\t")
        if(re.findall("^#",orlines.strip())):
            #print(orlines.strip())
            fw.write(orlines)
            mchri = orline.index("#chr")
            msi = orline.index("start")
            mei = orline.index("end")
            mgenei = orline.index("gene_name")
            minfoi = orline.index("infos")
            for keyindex in other.keys():
                recoderindex[keyindex] = orline.index(keyindex+"_result")
            proindex = orline.index(b[0]+"_result")
            continue
        mchr = orline[mchri]
        ms = orline[msi]
        me = orline[mei]
        mgene = orline[mgenei]
        mkey = "|".join([mchr,ms,me,mgene])
        orline[minfoi] = proinfo[mkey] #修改先证者info信息
        orline[proindex] = pro[mkey] #修改先征者可靠性
        #print(mkey)
        for samplei,values in other.items():#修改家系成员可靠性
            if(orline[recoderindex[samplei]-1] == "wild"):
                continue
            #recordflag = [] #相交区间可以是多个，记录这些可靠性
            overlapsize = [] #记录overlap区间的大小
            overlappro = []  #记录overlap区间的可靠性
            
            for comvalue in values:
                comchr,coms,come,comgene,comreun = comvalue
                if(compare(int(ms),int(me),int(coms),int(come)) and (mchr==comchr) and (mgene == comgene)):
                    #recordflag.append(comreun+"|".join(comvalue))
                    overlapsize.append(compare(int(ms),int(me),int(coms),int(come)))
                    overlappro.append(comreun)
                    #print(mkey,samplei)
            if(len(overlapsize) !=0):
                maxindex = overlapsize.index(max(overlapsize))
            #print(samplei)
                recordflag = overlappro[maxindex]
            #if(len(recordflag) !=0):
                #recordflag.append("red")
            #    orline[recoderindex[samplei]] = "|".join(recordflag)
                orline[recoderindex[samplei]] = recordflag
        #print(samplei)
        #print("\t".join(orline))
        fw.write("\t".join(orline)+"\n")
