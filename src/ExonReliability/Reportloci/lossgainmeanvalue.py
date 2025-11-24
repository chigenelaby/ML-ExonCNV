import sys
import re
import os
import numpy as np
from collections import defaultdict
import shutil
def triorelation(handle,proname):
	fp1 = open(handle,"r")
	relation = {}
	b = {}
	for tts in fp1:
		tt = tts.strip().split("\t")
		familyrole = tt[6] # 先证者，父，母
		taskcode = tt[3] #
		librarycode = tt[0]
		if(re.findall("父",familyrole) or re.findall("母",familyrole)): #or re.findall("先证者",familyrole)):
			b.setdefault(taskcode,[]).append([familyrole,librarycode])
		if(librarycode == proname):
			b.setdefault(taskcode,[]).append(["先证者",librarycode])
	#print(b)
	#print("task","先证者","母亲","父亲",sep="\t")
	for key,values in b.items():
		
		prohand = []
		mother = []
		father = []
		for gst in values:
			g,t = gst
			if(re.findall("先证者",g)):
				prohand.append([t])
			if(re.findall("母",g)):
				mother.append([t])
			if(re.findall("父",g)):
				father.append([t])
		prohandor = sorted(prohand, key=lambda x:x[0],reverse=True)
		motheror = sorted(mother, key=lambda x:x[0],reverse=True)
		fatheror = sorted(father, key=lambda x:x[0],reverse=True)
		#print(key,prohandor,motheror,fatheror)
		prohandp = ""
		motherp = ""
		fatherp = ""
		
		if(len(prohandor)!=0):
			prohandp = prohandor[0][0]
		if(len(motheror) !=0):
			motherp = motheror[0][0]
		if(len(fatheror) !=0):
			fatherp = fatheror[0][0]
		#print(key,prohandp,motherp,fatherp,sep="\t")
	relation[taskcode] = [prohandp,motherp,fatherp]
	return relation,prohandp 
		
#triorelation(sys.argv[1])	
def transformtable(file1,file2,proname):
	'''
	fp = open(sys.argv[1],"r")

	a = {}
	for i,lines in enumerate(fp):
		if(i == 0):
			continue
		line = lines.strip().split("\t")
		for j in line[2:]:
			if(j !=""):
				a.setdefault(line[1],[]).append(j)
	'''
	a = {}
	relation,prohandp = triorelation(file1,proname)
	for key,values in relation.items():
		for j in values:
			if(j !=""):
				a.setdefault(key,[]).append(j)	
	#for key,values in a.items():
	#	print(key,values)

	datan = list(a.keys())[0]  #get taskcode
	
	#print(a,datan)
	target = os.path.dirname(file2) + "/format_gene_info_back.txt"
	shutil.copyfile(file2, target)
	#print(target)	
	singleexon = os.path.dirname(file2) + '/singleexonfinetune.txt'
	#print(datan,singleexon)
	#with open(sys.argv[2],"r") as out:

	#b = {}
	b = defaultdict(lambda: defaultdict(list))
	with open(singleexon,"r") as out:
		for k,tts in enumerate(out):
			tt = tts.strip().split("\t")
			if(k == 0):
				chrnameJoini = tt.index("chrnameJoin")
				gene_namei = tt.index("gene_name")
				balancei = tt.index("balance")
				diffi = tt.index("diff")
				depthi = tt.index("depth")    
				sample_valuei = tt.index("sample_value")
				z_scorei =  tt.index("z_score")
				contrast_meani = tt.index("contrast.mean")  
				contrast_cvi = tt.index("contrast.cv")  
				contrast_stdi = tt.index("contrast.std")  
				contrast_dps_meani = tt.index("contrast.dps_mean")
				GCi = tt.index("GC")        
				mappabilityscorei = tt.index("mappabilityscore")  
				sizei = tt.index("size")
				fold80i = tt.index("FOLD80")
				#infosi = tt.index("infos")
				continue
			#info = tt[infosi]
			chrnameJoin = tt[chrnameJoini]	
			gene_name = tt[gene_namei]
			balance  = float(tt[balancei])
			diff = float(tt[diffi])
			depth = float(tt[depthi])
			sample_value = float(tt[sample_valuei])
			z_score = float(tt[z_scorei])
			contrast_mean = float(tt[contrast_meani])
			contrast_cv = float(tt[contrast_cvi])
			contrast_std = float(tt[contrast_stdi])
			contrast_dps_mean = float(tt[contrast_dps_meani])
			GC = float(tt[GCi])
			fold80 = float(tt[fold80i])
			mappabilityscore = float(tt[mappabilityscorei])
			size = float(tt[sizei])
			#if(re.findall("断点外显子支持",infos)):
			#	continue
			#sample[gene].setdefault(label,set()).add(name)
			b[chrnameJoin+"_"+gene_name]["balance"].append(balance)
			b[chrnameJoin+"_"+gene_name]["fold80"].append(fold80)
			b[chrnameJoin+"_"+gene_name]["diff"].append(diff)
			b[chrnameJoin+"_"+gene_name]["depth"].append(depth)
			b[chrnameJoin+"_"+gene_name]["sample_value"].append(sample_value)
			b[chrnameJoin+"_"+gene_name]["z_score"].append(z_score)
			b[chrnameJoin+"_"+gene_name]["contrast_mean"].append(contrast_mean)
			b[chrnameJoin+"_"+gene_name]["contrast_cv"].append(contrast_cv)
			b[chrnameJoin+"_"+gene_name]["contrast_std"].append(contrast_std)
			b[chrnameJoin+"_"+gene_name]["contrast_dps_mean"].append(contrast_dps_mean)
			b[chrnameJoin+"_"+gene_name]["GC"].append(GC)
			b[chrnameJoin+"_"+gene_name]["mappabilityscore"].append(mappabilityscore)
			b[chrnameJoin+"_"+gene_name]["size"].append(size)
	header = ["#chr","start","end","gene_name","gene_info_str","best_exon_str","freq","all_freq","infos","_tag","_result","lossgain","inheritance","ecount","balance","FOLD80","diff","depth","sample_value","z_score","contrast_mean","contrast_cv","contrast_std","contrast_dps_mean","GC","mappabilityscore","size"]
	fw = open(os.path.dirname(file2)+"/format_gene_info_table.txt","w")
	header = "\t".join(header)
	#print(header)
	fw.write(header+"\n")
	with open(file2,"r") as out1:
		for g,ggs in enumerate(out1):
			if(g==0 or g==1):
				continue
			gg = ggs.strip().split("\t")
			if(g==2):
				Q = []
				chri = gg.index("#chr")   
				starti = gg.index("start")   
				endi = gg.index("end")     
				gene_namemi = gg.index("gene_name") 
				freqi = gg.index("freq")    
				all_freqi = gg.index("all_freq")
				infosi = gg.index("infos")
				for L in a[datan]:
					Q.append(gg.index(L+"_tag"))
				continue
			infos = gg[infosi]
			if(re.findall("断点外显子支持",infos)):
				continue
			key = gg[chri]+":"+gg[starti]+"-"+gg[endi]+"_"+gg[gene_namemi]
			freq = gg[freqi]	
			all_freq = gg[all_freqi]	
			prohand = re.findall(r'gain|loss',gg[Q[0]])[0]
			count = 0
			label = 0
			gainlosslable = 0
			if(len(Q) > 1):
				for y in Q[1:]:
					if(re.findall(prohand,gg[y]) and gg[y] !="不可靠"):
						count += 1
			if(count > 0):
				label = 1
			else:
				label = 2
			if(prohand == "loss"):
				gainlosslable = 1
			if(prohand == "gain"):
				gainlosslable = 2
			Meanint = [gainlosslable,label,len(b[key]["balance"]),np.mean(b[key]["balance"]),np.mean(b[key]["fold80"]),np.mean(b[key]["diff"]),np.mean(b[key]["depth"]),np.mean(b[key]["sample_value"]),np.mean(b[key]["z_score"]),np.mean(b[key]["contrast_mean"]),np.mean(b[key]["contrast_cv"]),np.mean(b[key]["contrast_std"]),np.mean(b[key]["contrast_dps_mean"]),np.mean(b[key]["GC"]),np.mean(b[key]["mappabilityscore"]),np.mean(b[key]["size"])]
			Meanstr = [str(round(i,5)) for i in Meanint]
			temp = gg[:6] + gg[6:9] + [gg[Q[0]],gg[Q[0]+1]]+Meanstr
			fw.write("\t".join(temp)+"\n")
			#print("\t".join(temp))
	return prohandp
#transformtable(sys.argv[1],sys.argv[2])			
