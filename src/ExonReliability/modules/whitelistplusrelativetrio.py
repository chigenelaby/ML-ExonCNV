import sys
import re
import os

#01/11/2023
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
				goodcovregion[key] = 0
	return goodcovregion
				
	
#01/31/2024  相对可靠修改，家系支持，满足外显子个数以及平衡系数，白名单
def whitelist():
	white = os.path.join(os.path.dirname(os.path.abspath(__file__)),"whitelistrelativereliability.txt")
	region =[]
	with open(white,"r") as out:
		for lines in out:
			line = lines.strip().split("\t")
			region.append([line[0],int(line[1]),int(line[2])]) #Chr,start,end
	return region
			


def modificationrelative(mfile,palvalue,outname,sexoncv):
	#palvalue = {'ACD573':0.0019,'ACD574':0.0138,'ACD575':0.0059}
	rei = {}
	specialdealgain2 = ["gain2","gain"] # 特别可靠保留
	goodcovregion = coverageregion()
	outfile = os.path.join(os.path.dirname(mfile),outname)
	fw = open(outfile,"w")
	with open(mfile,"r") as fgi:
		for i,oris in enumerate(fgi):
			if(i == 0 or i == 1):
				fw.write(oris)
				continue
			ori = oris.strip().split("\t")
			if(i == 2): ##chr    start   end     gene_name       gene_info_str   best_exon_str   freq    all_freq        infos
				chri = ori.index("#chr")
				si = ori.index("start")
				ei = ori.index("end")
				gisi = ori.index("gene_info_str")
				infosi = ori.index("infos")
				fw.write(oris)
				for name in palvalue.keys():
					#tagi = ori.index(name +"_tag")
					rei[name] = ori.index(name +"_result")
				continue
			info = ori[infosi]
			Chr = ori[chri]
			st = int(ori[si])
			en = int(ori[ei])
			regionlegth = en - st + 1
			J = 0
			for tvalues in whitelist():
				cChr,cst,cen = tvalues
				if(Chr==cChr and cst <= st and cen >= en):
					J += 1
			#noL = 0
			#vrL = 0
			pse = "_".join([ori[chri],ori[si],ori[ei]])
			if(not(re.findall("-",ori[gisi]))):
				noL = 0
				vrL = 0
				judecv = 99999
				changerno = 88888 #记录不可靠的下标
				temp = []
				for key1,value1 in rei.items():
					for cvlist in sexoncv[key1]: #取单显子的cv值
						Chr1,b1,b2,cv = cvlist
						if(Chr1==Chr and st==b1 and en==b2):
							judecv = cv
							#temp.append(cv)			
					if(re.findall("不可靠",ori[value1]) and palvalue[key1] < 0.08 and re.findall("loss",ori[value1 -1]) and judecv < 0.15 and pse in goodcovregion):	
						noL += 1
						changerno = value1
					if((re.findall("相对可靠",ori[value1]) or (re.findall("特别可靠",ori[value1]))) and palvalue[key1] < 0.08 and re.findall("loss",ori[value1 -1]) and judecv < 0.173 and pse in goodcovregion):
						vrL += 1
				#print(ori[chri],ori[si],ori[ei],temp,noL,vrL)
				if(noL == 1 and vrL >= 1):
					ori[changerno] = "相对可靠"
					#print(ori[chri],ori[si],ori[ei],temp,noL,vrL)
				#fw.write("\t".join(ori)+"\n")
				#continue	
			if(re.findall("-",ori[gisi])):
				exons,exone = re.findall(r'\(REG:(\d+-\d+)\)',ori[gisi])[0].split("-")
				exonnumber = int(exone) - int(exons) + 1
				if(exonnumber > 3 and re.findall("家系成员支持",info) and re.findall("满足外显子个数条件",info) and J > 0):
					index = []
					for key,value in rei.items():
						if(palvalue[key] < 0.02	and not(re.findall("wild",ori[value -1])) and not(re.findall("不可靠",ori[value]))):
							index.append(value)
					if(len(index) > 1):
						#print("\t".join(ori))
						for u in index:
							ori[u] = "特别可靠"	
					
					#fw.write("\t".join(ori)+"\n")
				#else:
					#fw.write(oris)
			#2024/10/25 只要有纯合(双倍)、点突变支持、断点 其中之一的1个或<300bp  保留为特别可靠，其他1个或<300bp的最多为相对可靠
			if((not(re.findall("-",ori[gisi])) or regionlegth < 300)):
				for keysample,valueindex in rei.items():
					if(ori[valueindex] == "特别可靠" and ori[valueindex-1] not in specialdealgain2):
						if(re.findall("测序稳定区域",ori[infosi])):
							pass
						elif(re.findall("断点外显子支持",ori[infosi])):
							pass
						elif(re.findall("家系突变支持",ori[infosi])):
							pass
						else:
							ori[valueindex] = "相对可靠"
			fw.write("\t".join(ori)+"\n")
			#print("\t".join(ori))
			
#modificationrelative(sys.argv[1])

