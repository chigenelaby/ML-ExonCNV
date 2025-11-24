import argparse
from modules.expandinforfast  import *
# from modules.judgeconditionbasedMLsingleexon_xgboost import *
from modules.judgeconditionbasedMLsingleexon_xgboost_mutimodel import *
from modules.regionlossgainfinetune import *
from modules.relativesstriofinetune import *
from modules.calibrationbasedvariants import *
from modules.calibrationbasedtriovariants import *
from modules.whitelistplusrelativetrio import *
import sys
import os
import shutil

class CustomError(Exception):
    def __init__(self,ErrorInfo):
        super().__init__(self) #初始化父类
        self.errorinfo=ErrorInfo
    def __str__(self):
        return self.errorinfo


def finetune(args):
	a = args.input
	b = args.name
	c = args.svcf
	exonmode = args.exonmode
	balance = {} #均衡系数
	sexoncv = {} #外显子cv值
	for i,j,k in zip(a,b,c): #返回文件路径，样品名,vcf file
		path = os.path.dirname(i)
		#os.path.join(path,"demergetable.txt")
		balance[j] = demergeinfor(i,j, exonmode)
		sexoncv[j] = singleexonfinetune(os.path.join(path,"demergetable.txt"))
		
		# 新增xgboost可靠性校正
		fmt_xgboost_relative(os.path.join(path,"singleexonfinetune.txt"), os.path.join(path,"singleexonfinetune_fmt.txt"))

		finetuneregion(os.path.join(path,"singleexonfinetune_fmt.txt"),i,j)
		calivariant(os.path.join(path,"regionexonfinetune.txt"),j,k)#根据点突变，校正外显子缺失重复单样品
		#print(i,j)
	#print(balance)
	#以先证者角度合并外显子缺失重复
	if(len(a) > 1 and len(b) > 1 and len(a) == len(b)): #输入多个家系成员，需要对校正后的外显子缺失重复合并家系
		trio(a,b)
		#01/31/2024  相对可靠修改，家系支持，满足外显子个数以及平衡系数，白名单
		modificationrelative(os.path.join(os.path.dirname(a[0]),"format_gene_inforevis.txt"),balance,"format_gene_inforevis1.txt",sexoncv)
		untriosource = os.path.join(os.path.dirname(a[0]),"format_gene_inforevis1.txt")
		untriotarget = os.path.join(os.path.dirname(a[0]),"format_gene_info.txt")
		shutil.copy(untriosource,untriotarget)
	elif(len(a) == len(b) == 1):
		singlesource = os.path.join(os.path.dirname(a[0]),"regionexonfinetuneaddvariant.txt")
		singletarget = os.path.join(os.path.dirname(a[0]),"format_gene_info.txt")
		shutil.copy(singlesource,singletarget)
		#pass
	else: #输入的文件和输入的样品数目不一致
		raise CustomError("input error,maybe the number of input files are not equal the number files name or other reasons")
	#pass
#args = arg()
#finetune(args.input,args.name,args.svcf)
def finetunetrio(args):
	a = args.input
	b = args.name
	c = args.svcf
	d = args.mvcf
	exonmode = args.exonmode
	balance = {} #均衡系数
	sexoncv = {} #外显子cv值
	for i,j,k in zip(a,b,c): #返回文件路径，样品名,vcf file
		path = os.path.dirname(i)
		#os.path.join(path,"demergetable.txt")
		balance[j] = demergeinfor(i,j,exonmode)
		sexoncv[j] = singleexonfinetune(os.path.join(path,"demergetable.txt"))

		# 新增xgboost可靠性校正
		fmt_xgboost_relative(os.path.join(path,"singleexonfinetune.txt"), os.path.join(path,"singleexonfinetune_fmt.txt"))

		finetuneregion(os.path.join(path,"singleexonfinetune_fmt.txt"),i,j)
		calivariant(os.path.join(path,"regionexonfinetune.txt"),j,k) #根据点突变，校正外显子缺失重复单样品
		#print(i,j)
	#以先证者角度合并外显子缺失重复
	#print(balance)
	if(len(a) > 1 and len(b) > 1 and len(a) == len(b)): #输入多个家系成员，需要对校正后的外显子缺失重复合并家系
		trio(a,b)
	elif(len(a) == len(b) == 1):
		pass
	else: #输入的文件和输入的样品数目不一致
		raise CustomError("input error,maybe the number of input files are not equal the number files name or other reasons")
	calivarianttrio(a,b,d)
	#01/31/2024  相对可靠修改，家系支持，满足外显子个数以及平衡系数，白名单
	modificationrelative(os.path.join(os.path.dirname(a[0]),"format_gene_inforevistrio.txt"),balance,"format_gene_inforevistrio1.txt",sexoncv)
	triosource = os.path.join(os.path.dirname(a[0]),"format_gene_inforevistrio1.txt")
	triotarget = os.path.join(os.path.dirname(a[0]),"format_gene_info.txt")
	shutil.copy(triosource,triotarget)
	


parser = argparse.ArgumentParser()
subparsers = parser.add_subparsers(help='sub-command help')
#不包含完整家系
parser_a = subparsers.add_parser('notrio', help='can not contain family trio like proband,father,mother. this means that we only have two of three or one of three')
parser_a.add_argument('-input', help='delimited list input,format_gene_info.txt', type=lambda s: [item for item in s.split(',')],required=True)
parser_a.add_argument("-name",help='delimited list sample name, should be consistent with input file', type=lambda s: [item for item in s.split(',')],required=True)
parser_a.add_argument("-svcf",help='delimited list vcf file for single sample, should be consistent with input file ', type=lambda s: [item for item in s.split(',')],required=True)
parser_a.add_argument("-exonmode",help='exonmode', type=str)
#parser.add_argument("-mvcf",help='vcf file for the whole family members', type=str)
parser_a.set_defaults(func=finetune)

#包含完整家系
parser_b = subparsers.add_parser('trio', help='contain family trio like proband,father,mother. this means that we have them all ')
parser_b.add_argument('-input', help='delimited list input,format_gene_info.txt', type=lambda s: [item for item in s.split(',')],required=True)
parser_b.add_argument("-name",help='delimited list sample name, should be consistent with input file', type=lambda s: [item for item in s.split(',')],required=True)
parser_b.add_argument("-svcf",help='delimited list vcf file for single sample, should be consistent with input file ', type=lambda s: [item for item in s.split(',')],required=True)
parser_b.add_argument("-mvcf",help='vcf file for the whole family members', type=str,required=True)
parser_b.add_argument("-exonmode",help='exonmode', type=str)
parser_b.set_defaults(func=finetunetrio)
args = parser.parse_args()
#执行函数功能
args.func(args)
