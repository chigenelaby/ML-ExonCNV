import sys
import re
import ast
import json
import os
import argparse
from  lossgainmeanvalue import *

def pgene():
	genelist = {}
	hgenelist = {}
	midgenelist = {}
	#print(file1)
	file1 = os.path.join(os.path.dirname(os.path.abspath(__file__)),"positivegenelist.txt")
	with open(file1,"r") as out:
		for i,lines in enumerate(out):
			line = lines.strip().split("\t")
			#print(lines)
			if(i == 0):
				genei = line.index("基因")
				numsamplei = line.index("阳性数")
				continue
			gene = line[genei]
			numsample = int(line[numsamplei])
			genelist[gene] = 0
			if(numsample >= 20):
				hgenelist[gene] = 0
			if(numsample >= 2):
				midgenelist[gene] = 0
			#genelist[gene] = 0
			
	return genelist,hgenelist,midgenelist



def whiteregion():
	#gene_attentions.txt_1：针对基因视角的白名单文件，只要是这个文件中的基因就保留
	#gene_attentions.txt_2:   针对基因视角的白名单文件，这个文件中的基因需要满足长度>=3个外显子，才保留
	focus1 = os.path.join(os.path.dirname(os.path.abspath(__file__)),"gene_attentions.txt_1")
	fp1 = open(focus1,"r")
	f1 = {}
	f2 = {}
	for i,f1lines in enumerate(fp1):
		if(i==0):
			continue
		f1line = f1lines.strip().split("\t")
		key = "_".join(f1line[:3])
		f1[key] = f1line[3]	
	focus2 = os.path.join(os.path.dirname(os.path.abspath(__file__)),"gene_attentions.txt_2")
	fp2 = open(focus2,"r")
	for j,f2lines in enumerate(fp2):
		if(j==0):
			continue
		f2line = f2lines.strip().split("\t")
		key = "_".join(f2line[:3])
		f2[key] = f2line[3]
	return f1,f2

#严重不均衡的样品只报白名单上的
#task = sys.argv[1].split("/")[-1].split(".")[0]
def judgeal(file1, freq_cutoff1, freq_cutoff2):
	#whitelistdict = {}
	setredict = {}
	setredictlow = {} #lowreport_loci
	#f1dict,f2dict = whiteregion()
	genelist, hgenelist, midgenelist = pgene()
	whitelistgene = ["DMD","CFHR1","CFHR3"]
	with open(file1,"r") as out:
		for i,lines in enumerate(out):	
			line = lines.strip().split("\t")
			if(i==0):
				Chri = line.index("#chr")
				starti = line.index("start")
				endi = line.index("end")
				genei = line.index("gene_name")
				freqi = line.index("freq")
				diffi = line.index("diff")
				z_scorei = line.index("z_score")
				sample_valuei = line.index("sample_value")
				ecounti = line.index("ecount")
				mappabilityscorei = line.index("mappabilityscore")
				GCi = line.index("GC")
				resulti = line.index("_result")
				infosi = line.index("infos")
				all_freqi = line.index("all_freq")
				contrast_stdi = line.index("contrast_std")
				inheritancei = line.index("inheritance")
				balancei = line.index("balance")
				contrast_cvi = line.index("contrast_cv") #
				contrast_dps_meani = line.index("contrast_dps_mean") #
				FOLD80i = line.index("FOLD80") #
				depthi = line.index("depth") #
				infosi = line.index("infos")
				gene_namei = line.index("gene_name")
				reliresulti = line.index("_result")
				continue
			infos = line[infosi]
			#if(re.findall("断点外显子支持",infos)):
			#	setredict[key] = "R1"
			#	line = line + ["R1"]
			#	continue
			key = "_".join(line[:5])
			contrast_cv  = float(line[contrast_cvi])
			contrast_dps_mean = float(line[contrast_dps_meani])
			FOLD80 = float(line[FOLD80i])
			depth = float(line[depthi])
			gene_name  = line[gene_namei]
			Chr = line[Chri]
			start = int(line[starti])
			end = int(line[endi])
			gene = line[genei]
			freq = float(line[freqi])
			reliresult = line[reliresulti]
			diff = float(line[diffi])
			z_score = float(line[z_scorei])
			sample_value = float(line[sample_valuei])
			ecount = float(line[ecounti])
			mappabilityscore = float(line[mappabilityscorei])
			GC = float(line[GCi])
			result = line[resulti]
			infos = line[infosi]
			all_freq = float(line[all_freqi])
			contrast_std = float(line[contrast_stdi])
			inheritance = float(line[inheritancei])
			balance = float(line[balancei])
			#key = "_".join(line[:5])
			#elen = end - start +1
			length = int(line[endi]) - int(line[starti])
			temp = json.loads(line[infosi])#
			temp['report_tag'] = 'NA'#
			line[infosi] = json.dumps(temp,ensure_ascii=False)#
			if(re.findall("nogene",gene_name)):
				continue
			if(re.findall("断点外显子支持",line[infosi])):
				#reline = addreportloic(infosi,line)
				#line = reline + ["R1"]
				setredict[key] = "R1"
				continue
			if(balance >= 0.08):
				if(re.findall("测序稳定区域",line[infosi])):
					#reline = addreportloic(infosi,line)
					#line = reline + ["R6"]
					setredict[key] = "R6"
				#if(re.findall("断点外显子支持",line[infosi])):
				#	reline = addreportloic(infosi,line)
				#	line = reline + ["R6"]
				if(re.findall("家系突变支持",line[infosi])):
					#reline = addreportloic(infosi,line)
					#line = reline + ["R6"]
					setredict[key] = "R6"
				#if(len(line) != rownumber):
				if(gene_name == "DMD"):
						#print("\t".join(line + ["R5"]))
					setredict[key] = "R7"
					#else:
					#	pass
						#print("\t".join(line + ["NA"]))
				#else:
					#pass
					#print("\t".join(line))
			else:
					
				if(reliresult == "特别可靠"):
					if(freq < freq_cutoff1):
						#reline = addreportloic(infosi,line)
						#line = reline + ["R1a"]
						setredict[key] = "R1a"
					if(freq >= freq_cutoff1 and freq < freq_cutoff2 and (ecount == 1 or length <300) and (inheritance == 1 or gene_name in genelist)):
						#reline = addreportloic(infosi,line)
						#line = reline + ["R1b"]
						setredict[key] = "R1b"
					if(freq >= freq_cutoff1 and freq < freq_cutoff2 and (ecount >= 2 and length >=300) and gene_name in genelist):
						#reline = addreportloic(infosi,line)
						#line = reline + ["R1c"]
						setredict[key] = "R1c"
					if(freq >= freq_cutoff1 and freq < freq_cutoff2 and (ecount >= 2 and length >=300) and gene_name not in genelist):
						setredictlow[key] = "R1c"

					if(freq < 0.05 and freq >= freq_cutoff2 and gene_name in genelist):
						#reline = addreportloic(infosi,line)
						#line = reline + ["R1d"]
						setredict[key] = "R1d"
				elif(reliresult == "相对可靠"):
					if(balance < 0.02):
						if(freq < freq_cutoff1 and (ecount == 1 or length <300) and ((inheritance == 1 and gene_name in midgenelist) or gene_name in midgenelist)):
							#if((ecount == 1 or length <300) and (inheritance == 1 or gene_name in genelist)):
							#reline = addreportloic(infosi,line)
							#line = reline + ["R2a"]
							setredict[key] = "R2a"
							#if((ecount >= 2 and length >=300)):
						elif(freq < freq_cutoff1 and (ecount >= 2 and length >=300) and gene_name in midgenelist):
							#reline = addreportloic(infosi,line)
							#line = reline + ["R2b"]
							setredict[key] = "R2b"
						elif(freq < freq_cutoff1 and (ecount >= 2 and length >=300) and gene_name not in midgenelist):
							setredictlow[key] = "R2b"

						elif(freq < 0.05 and freq >= freq_cutoff1 and (ecount == 1 or length <300) and ((inheritance == 1 and gene_name in midgenelist) or gene_name in hgenelist)):
							#if((ecount == 1 or length <300) and (inheritance == 1 or gene_name in hgenelist)):
							#reline = addreportloic(infosi,line)
							#line = reline + ["R2c"]
							setredict[key] = "R2c"
							#if((ecount == 2 and length >=300) and (inheritance == 1 or gene_name in genelist)):
						elif(freq < 0.05 and freq >= freq_cutoff1 and (ecount == 2 and length >=300) and (inheritance == 1 or gene_name in genelist)):
							#reline = addreportloic(infosi,line)
							#line = reline + ["R2d"]
							#setredict[key] = "R2d"
							if(gene_name in midgenelist):
								setredict[key] = "R2d"
							else:
								setredictlow[key] = "R2d"
						elif(freq >= freq_cutoff1 and freq < freq_cutoff2 and (ecount >= 3 and length >=300) and gene_name in midgenelist):
							#reline = addreportloic(infosi,line)
							#line = reline + ["R2e"]
							setredict[key] = "R2e"
						elif(freq >= freq_cutoff1 and freq < freq_cutoff2 and (ecount >= 3 and length >=300) and gene_name not in midgenelist):
							setredictlow[key] = "R2e"
						elif(freq < 0.05 and freq >= freq_cutoff2 and (ecount >= 3 and length >=300) and (gene_name in genelist)):
							#reline = addreportloic(infosi,line)
							#line = reline + ["R2f"]
							setredict[key] = "R2f"
						else:
							pass
					if(balance >= 0.02 and balance < 0.04):
						if(freq < freq_cutoff1 and (ecount < 3 or length <300) and ((inheritance == 1 and gene_name in midgenelist )or gene_name in midgenelist)):
							#if((ecount < 3 or length <300) and (inheritance == 1 or gene_name in genelist)):
							#reline = addreportloic(infosi,line)
							#line = reline + ["R3a"]
							setredict[key] = "R3a"
							#if(ecount >= 3):
						elif(freq < freq_cutoff1 and ecount >= 3 and gene_name in midgenelist):
							#reline = addreportloic(infosi,line)
							#line = reline + ["R3b"]
							setredict[key] = "R3b"
						elif(freq < freq_cutoff1 and ecount >= 3 and gene_name not in midgenelist):
							setredictlow[key] = "R3b"
						elif(freq >= freq_cutoff1 and freq < freq_cutoff2 and (ecount == 1 or length <300) and ((inheritance == 1 and gene_name in midgenelist ) or gene_name in hgenelist)):
							#if((ecount == 1 or length <300) and (inheritance == 1 or gene_name in hgenelist)):
							#reline = addreportloic(infosi,line)
							#line = reline + ["R3c"]
							setredict[key] = "R3c"
							#if((ecount >= 2 and length >=300) and ( gene_name in genelist or gene_name in genelist)):
						elif(freq >= freq_cutoff1 and freq < freq_cutoff2 and (ecount >= 2 and length >=300) and (inheritance == 1 or gene_name in genelist)):
							#reline = addreportloic(infosi,line)
							#line = reline + ["R3d"]
							if(gene_name in midgenelist):
								setredict[key] = "R3d"
							else:
								setredictlow[key] = "R3d"
						elif(freq < 0.05 and freq >= freq_cutoff2 and (ecount < 3 and length >=300) and gene_name in hgenelist):
							#reline = addreportloic(infosi,line)
							#line = reline + ["R3e"]
							setredict[key] = "R3e"
						elif(freq < 0.05 and freq >= freq_cutoff2 and (ecount >= 3 and length >=300) and gene_name in genelist):
							#reline = addreportloic(infosi,line)
							#line = reline + ["R3f"]
							setredict[key] = "R3f"
						else:
							pass
					if(balance >= 0.04 and balance < 0.08):
						if(freq < freq_cutoff1 and (ecount < 3 or length < 300) and gene_name in hgenelist):
							#if((ecount < 3 or length < 300) and gene_name in hgenelist):
							#reline = addreportloic(infosi,line)
							#line = reline + ["R4a"]
							setredict[key] = "R4a"
							#if((ecount >= 3 and length >= 300)and (inheritance == 1 or gene_name in genelist)):
						elif(freq < freq_cutoff1 and (ecount >= 3 and length >= 300)and (inheritance == 1 or gene_name in genelist)):
							#reline = addreportloic(infosi,line)
							#line = reline + ["R4b"]
							setredict[key] = "R4b"
						elif(freq >= freq_cutoff1 and freq < freq_cutoff2 and (inheritance == 1 or gene_name in hgenelist)):
							#reline = addreportloic(infosi,line)
							#line = reline + ["R4c"]
							setredict[key] = "R4c"
						elif(freq < 0.05 and freq >= freq_cutoff2 and gene_name in hgenelist):
							#reline = addreportloic(infosi,line)
							#line = reline + ["R4d"]
							setredict[key] = "R4d"
						else:
							pass
				else:
					pass		
					
						
				#if(len(line) != rownumber):
				if(gene_name in whitelistgene):
					setredict[key] = "R5"
						#print("\t".join(line + ["R5"]))
					#else:
						#pass
						#print("\t".join(line + ["NA"]))	
				#else:
					#pass
					#print("\t".join(line))		
			#print(rownumber)
			#print(len(line))

	 
			#if(balance < 0.02 and reportlocigoodbalance(freq,diff,z_score,sample_value,ecount,mappabilityscore,GC,all_freq,contrast_std,inheritance,contrast_cv,contrast_dps_mean,FOLD80,depth) and result == "相对可靠"):
				#print(line)
				#mre += 1
				#key = "_".join(line[:5])
				#setredict[key] = 0
			#if((balance >= 0.02 and balance < 0.08 ) and reportlociunbalance(freq,diff,z_score,sample_value,ecount,mappabilityscore,GC,all_freq,contrast_std,inheritance,contrast_cv,contrast_dps_mean,FOLD80,depth) and result == "相对可靠"):	
				#key = "_".join(line[:5])
				#setredict[key] = 0
						
	return setredict,setredictlow


parser = argparse.ArgumentParser()
parser.add_argument('input_file', help='format_gene_info文件')
parser.add_argument('handle_file', help='handle文件')
parser.add_argument('wkcode', help='文库号')
parser.add_argument('--freq_cutoff1', type=float, default=0.005, help='可选参数，默认0.001')
parser.add_argument('--freq_cutoff2', type=float, default=0.01, help='可选参数，默认0.005')
parser.add_argument('--filter_info', type=bool, default=False, help='可选参数，默认False')


args = parser.parse_args()

input_file = args.input_file
handle_file = args.handle_file
wkcode = args.wkcode
freq_cutoff1 = args.freq_cutoff1
freq_cutoff2 = args.freq_cutoff2
is_filter_info = args.filter_info

name = transformtable(handle_file,input_file,wkcode)

fp = open(os.path.dirname(input_file)+"/format_gene_info_back.txt","r")
setredict,setredictlow = judgeal(os.path.dirname(input_file)+"/format_gene_info_table.txt", freq_cutoff1, freq_cutoff2)
#print(whitelistdict)
fw = open(os.path.dirname(input_file)+"/format_gene_info.txt","w")
#name = transformtable(sys.argv[1],sys.argv[2])
for j,lines in enumerate(fp):
	if(j < 2):
		fw.write(lines)
		continue
	line = lines.strip().split("\t")
	#print(line)
	#name = line[-3]
	if(j == 2):
		fw.write(lines)
		#rei = line.index(name +"_result")
		infosi = line.index("infos")
		agi = line.index(name +"_tag")
		relii = line.index(name +"_result")
		continue
	#re = line[rei]
	ag = line[agi]
	key = "_".join(line[:5])
	infos = line[infosi]
	#if(len(setredict) != 0): 
	if(ag !='wild' and key in setredict):
		#print(line[infosi])
		#line[infosi]['report_tag'] = 'Report_loci'
		#temp = ast.literal_eval(line[infosi])
		temp = json.loads(line[infosi])
		temp['report_tag'] = 'Report_loci'
		temp['report_flag'] = setredict[key]
		#print(temp)
		line[infosi] = json.dumps(temp,ensure_ascii=False)
		#print(line[infosi])
		#print(line)
	#elif(re =="相对可靠" and key in setredict and ag !='wild'):
	#	temp = json.loads(line[infosi])
	#	temp['report_tag'] = 'Report_loci'
	#	temp['report_flag'] = setredict[key]
	#	line[infosi] = json.dumps(temp,ensure_ascii=False)
	#if(re =="相对可靠" and key not in setredict):
	#	temp = json.loads(line[infosi])
	#	temp['report_tag'] = 'NA'
	#	line[infosi] = json.dumps(temp,ensure_ascii=False)
	#if(re =="不可靠"):
	#	temp = json.loads(line[infosi])
	#	temp['report_tag'] = 'NA'
	#	line[infosi] = json.dumps(temp,ensure_ascii=False)
	elif(re.findall("断点外显子支持",line[infosi])):
		temp = json.loads(line[infosi])
		temp['report_tag'] = 'Report_loci'
		temp['report_flag'] = "断点"
		line[infosi] = json.dumps(temp,ensure_ascii=False)
	elif(ag !='wild' and key in setredictlow):
		temp = json.loads(line[infosi])
		temp['report_tag'] = 'lowReport_loci'
		temp['report_flag'] = setredictlow[key]
		line[infosi] = json.dumps(temp,ensure_ascii=False)
	else:
		temp = json.loads(line[infosi])
		temp['report_tag'] = 'NA'
		line[infosi] = json.dumps(temp,ensure_ascii=False)
	#if(key in whitelistdict):
	#	temp = json.loads(line[infosi])
	#	temp['report_tag'] = 'Report_loci'
	#	line[infosi] = json.dumps(temp,ensure_ascii=False)
	if is_filter_info:
		extracted_data = {}
		reli_dic = {'特别可靠':"highly_reliable", "相对可靠": "relatively_reliable", "不可靠":"unreliable"}
		loaded_dict = json.loads(line[infosi])
		keys_to_extract = ['report_tag']
		extracted_data = {key: loaded_dict[key] for key in keys_to_extract if key in loaded_dict}
		line[infosi] = json.dumps(extracted_data, ensure_ascii=False)
		if line[relii] in reli_dic:
			line[relii] = reli_dic[line[relii]]
	fw.write("\t".join(line)+"\n")
	#infos = line[infosi]
	
#for key in setredict.keys():
#	print(key)		
#print(task,mver,mre,0,mver+mre,pver,pre,pu,pver+pre+pu,sep="\t")
