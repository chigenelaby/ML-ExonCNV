import sys
import re
import json


def pgene(file1):
	genelist = {}
	hgenelist = {}
	#print(file1)
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
			#genelist[gene] = 0
			
	return genelist,hgenelist


def addreportloic(infosi,reline):

	temp = json.loads(reline[infosi])
	temp['report_tag'] = 'Report_loci'
	reline[infosi] = json.dumps(temp,ensure_ascii=False)
	return reline

def reportlocig(file1,file2):
	rownumber = 0
	genelist,hgenelist = pgene(file1)
	with open(file2,"r") as out:
		for j,lines in enumerate(out):
			line = lines.strip().split("\t")
			if(j==0):
				print("\t".join(line +['label']))
				freqi = line.index("freq")
				reliresulti = line.index("_result")
				inheritancei = line.index("inheritance")
				ecounti = line.index("ecount")
				balancei = line.index("balance")
				gene_namei = line.index("gene_name")
				infosi = line.index("infos")
				rownumber = len(line) + 1
				continue
			freq = float(line[freqi])
			reliresult = line[reliresulti]
			inheritance = int(line[inheritancei])
			ecount  = int(line[ecounti])
			balance  = float(line[balancei])
			gene_name  = line[gene_namei]
			temp = json.loads(line[infosi])#
			temp['report_tag'] = 'NA'#
			line[infosi] = json.dumps(temp,ensure_ascii=False)#
			if(balance > 0.08):
				continue
			if(reliresult == "特别可靠"):
				if(freq < 0.005):
					reline = addreportloic(infosi,line)
					line = reline + ["R1a"]
				if(freq >= 0.005 and gene_name in genelist):
					reline = addreportloic(infosi,line)
					line = reline + ["R1b"]
			elif(reliresult == "相对可靠"):
				if(balance <= 0.02):
					if(freq <= 0.001):
						if(ecount >= 2):
							reline = addreportloic(infosi,line)
							line = reline + ["R2a"]
						if(ecount==1 and gene_name in genelist):
							reline = addreportloic(infosi,line)
							line = reline + ["R2a"]
					elif(freq >= 0.001 and ecount < 3):
						if(ecount == 2 and (inheritance == 1 or gene_name in genelist)):
							reline = addreportloic(infosi,line)
							line = reline + ["R2b"]
						if(ecount == 1 and gene_name in genelist):
							reline = addreportloic(infosi,line)
							line = reline + ["R2b"]
							
					elif(freq >= 0.001 and freq < 0.005 and ecount >= 3):
						reline = addreportloic(infosi,line)
						line = reline + ["R2c"]
					elif(freq >= 0.005 and ecount >= 3 and (gene_name in genelist)):
						reline = addreportloic(infosi,line)
						line = reline + ["R2d"]
					else:
						pass
				if(balance > 0.02 and balance <= 0.04):
					if(freq < 0.001 and ecount < 3 and (inheritance == 1 or gene_name in genelist)):
						reline = addreportloic(infosi,line)
						line = reline + ["R3a"]
					elif(freq < 0.001 and ecount >= 3):
						reline = addreportloic(infosi,line)
						line = reline + ["R3b"]
					elif(freq >= 0.001 and freq < 0.005):
						if(ecount >= 2 and (inheritance == 1 or gene_name in genelist)):
							reline = addreportloic(infosi,line)
							line = reline + ["R3c"]
						if(ecount == 1 and gene_name in genelist):
							reline = addreportloic(infosi,line)
							line = reline + ["R3c"]
					elif(freq >= 0.005 and ecount >= 3 and gene_name in genelist):
						reline = addreportloic(infosi,line)
						line = reline + ["R3d"]
					elif(freq >= 0.005 and ecount < 3 and gene_name in hgenelist):
						reline = addreportloic(infosi,line)
						line = reline + ["R3e"]
					else:
						pass
				if(balance > 0.04 and balance <= 0.08):
					if(freq < 0.001 and ecount < 3 and gene_name in hgenelist):
						reline = addreportloic(infosi,line)
						line = reline + ["R4a"]
					elif(freq < 0.001 and ecount >= 3 and (inheritance == 1 or gene_name in genelist)):
						reline = addreportloic(infosi,line)
						line = reline + ["R4b"]
					elif(freq >= 0.001 and freq < 0.005 and (inheritance == 1 or gene_name in hgenelist)):
						reline = addreportloic(infosi,line)
						line = reline + ["R4c"]
					elif(freq >= 0.005 and gene_name in hgenelist):
						reline = addreportloic(infosi,line)
						line = reline + ["R4d"]
					else:
						pass
			else:
				pass		
					
						
			if(len(line) != rownumber):
				print("\t".join(line + ["NA"]))	
			else:
				print("\t".join(line))		
			#print(rownumber)
			#print(len(line))
				
				
reportlocig(sys.argv[1],sys.argv[2])
