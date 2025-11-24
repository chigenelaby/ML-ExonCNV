import argparse
import gzip
import json
import sys
import os
from collections import defaultdict, OrderedDict
import time
import re

parser = argparse.ArgumentParser(description='summary total reads, depth, covearage')
parser.add_argument('-outdir', help='out dir', required=True)
parser.add_argument('-prefix', help='output file prefix', required=True)
parser.add_argument('-capture', help='get how many capture reads such as capture.txt', default="-")
parser.add_argument('-dup', help='get how many duplication such as dup_metrics.txt', default="-")
parser.add_argument('-Insert', help='search for insert size like Insert_metrics.txt', default="-")
parser.add_argument('-json', help='json file produced by fastp like fastp.json', default="-")
parser.add_argument('-map_ratio', help='map file contains ratio of mapping,PE ratio as map_metrics.txt', default="-")
parser.add_argument('-wgs_depth', help='like chrALL.regions.bed.gz', default="-")
parser.add_argument('-wgs_cov', help='acquire covearage in different depth like chrALL.thresholds.bed.gz', default="-")
parser.add_argument('-gene_depth', help='like gene.merged.regions.bed.gz', default="-")
parser.add_argument('-gene_cov', help='acquire covearage in different depth like gene.merged.thresholds.bed.gz', default="-")

# target-region list
parser.add_argument('-target_depth', help='like NT01_target.regions.bed.gz', default="-")
parser.add_argument('-target_cov', help='acquire covearage in different depth like NT01_target.thresholds.bed.gz', default="-")
parser.add_argument('-target_bed', help='like NT01_target.bed', default="-")
args = parser.parse_args()

# time_s = time.time()


def summary_control_quality(a, b, c, d, e):
    reads_gc = OrderedDict()

    # 1. fastp.json
    if(a == "-"):
        reads_gc['raw reads'] = "-"
        reads_gc['raw data(G)'] = "-"
        reads_gc['raw GC%'] = "-"
        reads_gc['raw Q20%'] = "-"
        reads_gc['raw Q30%'] = "-"
        reads_gc['clean reads'] = "-"
        reads_gc['clean data(G)'] = "-"
        reads_gc['clean data GC%'] = "-"
        reads_gc['clean Q20%'] = "-"
        reads_gc['clean Q30%'] = "-"
    else:
        fp = open(a, 'r')
        reads_bases_gc = json.load(fp)
        reads_gc['raw reads'] = reads_bases_gc['summary']['before_filtering']['total_reads']
        reads_gc['raw data(G)'] = f"{reads_bases_gc['summary']['before_filtering']['total_bases']/1E+9:.2f}"
        reads_gc['raw GC%'] = f"{reads_bases_gc['summary']['before_filtering']['gc_content']*100:.2f}"
        reads_gc['raw Q20%'] = f"{reads_bases_gc['summary']['before_filtering']['q20_rate']*100:.2f}"
        reads_gc['raw Q30%'] = f"{reads_bases_gc['summary']['before_filtering']['q30_rate']*100:.2f}"
        reads_gc['clean reads'] = reads_bases_gc['summary']['after_filtering']['total_reads']
        reads_gc['clean data(G)'] = f"{reads_bases_gc['summary']['after_filtering']['total_bases']/1E+9:.2f}"
        reads_gc['clean data GC%'] = f"{reads_bases_gc['summary']['after_filtering']['gc_content']*100:.2f}"
        reads_gc['clean Q20%'] = f"{reads_bases_gc['summary']['after_filtering']['q20_rate']*100:.2f}"
        reads_gc['clean Q30%'] = f"{reads_bases_gc['summary']['after_filtering']['q30_rate']*100:.2f}"
        fp.close()

    # 2. map_metrics.txt
    PF_READS_ALIGNED = 0
    if(b == "-"):
        reads_gc['map%'] = "-"
        reads_gc['PE%'] = "-"
    else:
        with open(b, "r") as read_map:
            for i, lines in enumerate(read_map):
                line = lines.strip().split()
                # if(i == 0):
                #     continue
                if len(line)>17 and line[0] == "PAIR":
                    #picard/gtx/sentieon/gatk
                    mapping_ratio = float(line[6])
                    reads_gc['map%'] = f"{mapping_ratio*100:.2f}"
                    PF_READS_ALIGNED = int(line[5])
                    PE_ratio = float(line[6]) * float(line[17]) * 100
                    reads_gc['PE%'] = f"{PE_ratio:.2f}"
                elif re.search(r"mapped\s+\(", lines): 
                    #samtools/sambamba flagstat
                    reads_gc['map%'] = float(re.sub(r".*\(|%.*","", lines))
                    PF_READS_ALIGNED = int(re.sub(r"\s+\+.*", "", lines))
                elif re.search("singletons", lines):
                    #samtools/sambamba flagstat
                    se_ratio = re.sub(r".*\(|%.*","",lines)
                    reads_gc['PE%'] = f"{float(reads_gc['map%'])- float(se_ratio):.2f}"

    # 5. capture.txt
    if (e == "-" or b == "-"):
        reads_gc['捕获效率'] = "-"
    else:
        with open(e, 'r') as reads_capture:
            capture_num = 0
            for line in reads_capture:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                capture_num = int(line)
            reads_gc['捕获效率'] = f"{capture_num/PF_READS_ALIGNED:.4f}"

    # 3.dup_metrics.txt
    if(c == "-"):
        reads_gc["dup%"] = "-"
    else:
        with open(c, "r") as read_dup:
            for i, lines in enumerate(read_dup):
                line = lines.strip().split("\t")
                if(i == 0):
                    continue
                if(i == 2) and len(line) >=8: #sentieon
                    reads_gc["dup%"] = f"{float(line[8])*100:.2f}"
                    break
                if len(line) >=1 and line[0] == "percent_duplication": #gtx
                    reads_gc["dup%"] = f"{float(line[1])*100:.2f}"
                    break
                if i==4 and len(line) >=8: #bamsormadup
                    reads_gc["dup%"] = f"{float(line[7])*100:.2f}"
                    break

    # 4. Insert_metrics.txt
    if(d == "-"):
        reads_gc['插入片段'] = "-"
    else:
        with open(d, "r") as read_Insert:
            Insertsize = 0
            for i, lines in enumerate(read_Insert):
                line = lines.strip().split("\t")
                if "MEAN_INSERT_SIZE" in line:
                    Insertsize=line.index("MEAN_INSERT_SIZE")
                    continue
                if len(line) >18 and not lines.startswith("MEDIAN"): #diff soft diff line
                    reads_gc['插入片段'] =  f"{float(line[Insertsize]):.2f}"
                    break
    return reads_gc


def Gene_stat(gd, g):
    gene_gc = defaultdict(dict)

    # 5.1. depth file gene
    if(gd == "-"):
        pass
    else:
        with gzip.open(gd, 'rt') as gene_dep:
            for lines in gene_dep:
                _, start, end, gene, depth = lines.strip().split("\t")
                region = int(end) - int(start)
                gene_gc[gene]['tot_region'] = gene_gc[gene].get('tot_region', 0) + region
                gene_gc[gene]['tot_depth'] = gene_gc[gene].get('tot_depth', 0) + float(depth) * region

    # 6.1 covearage file gene_metrics
    if g == "-":
        pass
    else:
        with gzip.open(g, "rt") as read_cov:
            for i, lines in enumerate(read_cov):
                line = lines.strip().split("\t")
                # #chrom	start	end	region	1X	2X	3X	5X	10X	15X	20X	25X	30X	40X	50X
                if(i == 0):
                    continue
                gene = line[3]
                gene_gc[gene]['len'] = gene_gc[gene].get('len', 0) + (int(line[2]) - int(line[1]))
                gene_gc[gene]['cover1X'] = gene_gc[gene].get('cover1X', 0) + int(line[4])
                gene_gc[gene]['cover10X'] = gene_gc[gene].get('cover10X', 0) + int(line[8])
    return gene_gc


def sample_summary(e, f, reads_gc):

    # 5. depth file sample_summary
    if(e == "-"):
        reads_gc['平均深度'] = "-"
    else:
        with gzip.open(e, 'rt') as wgs_depth:
            tot_region = 0
            tot_depth = 0
            for lines in wgs_depth:
                line=lines.strip().split("\t")
                if len(line) == 4 :
                    _, start, end, depth = lines.strip().split("\t")
                    tot_region += int(end) - int(start)
                    tot_depth += float(depth) * ( int(end) - int(start) )
                else:
                    _, start, end, region, depth = lines.strip().split("\t")
                    tot_region += int(region)
                    tot_depth += float(depth) * int(region)
            reads_gc['平均深度'] = f"{tot_depth / tot_region:.2f}"

    # 6. covearage file wgs_metrics
    if f == "-":
        reads_gc['覆盖度1X'] = "-"
        reads_gc['覆盖度5X'] = "-"
        reads_gc['覆盖度10X'] = "-"
        reads_gc['覆盖度15X'] = "-"
        reads_gc['覆盖度20X'] = "-"
        reads_gc['覆盖度30X'] = "-"
        reads_gc['覆盖度40X'] = "-"
        reads_gc['覆盖度50X'] = "-"
    else:
        cover_dict = defaultdict(dict)
        with gzip.open(f, "rt") as read_cov:
            title_list = []
            for i, lines in enumerate(read_cov):
                line = lines.strip().split("\t")
                if(i == 0):
                    title_list = line[4:15]
                    continue
                for j, cover in enumerate(title_list):
                    index = j + 4
                    cover_dict[cover]['cover'] = cover_dict[cover].get('cover', 0) + int(line[index])
                    if line[3] != "unknown":
                        cover_dict[cover]['len'] = cover_dict[cover].get('len', 0) + int(line[3])
                    else:
                        cover_dict[cover]['len'] = cover_dict[cover].get('len', 0) + int(line[2]) - int(line[1])
            for cc in cover_dict.keys():
                new_key = '覆盖度' + cc
                reads_gc[new_key] = f"{cover_dict[cc]['cover'] / int(cover_dict[cc]['len']):.4f}"


def write_data(reads_gc, gene_gc, h):

    # output QC summary
    with open(f"{h}_summary", "w") as out:
        # k = 0
        # t = 0
        header=[]
        values=[]
        for key,value in reads_gc.items():
            header.append(key)
            values.append(str(value))
        out.write("\t".join(header) + "\n")
        out.write("\t".join(values) + "\n")
        #for key in reads_gc.keys():
         #   k = k + 1
          #  if(k < len(reads_gc)):
          #      out.write(key + "\t")
          #  else:
          #      out.write(key + "\n")
        #for key, value in reads_gc.items():
        #    t = t + 1
        #    if(t < len(reads_gc)):
        #        out.write(str(value) + "\t")
        #    else:
        #        out.write(str(value) + "\n")

    if gene_gc:
        # output gene stat
        with open(f"{h}_gene", "w") as out:
            title = "基因名\t基因平均深度\t覆盖度(1X)\t覆盖度(10X)"
            out.write(title + '\n')
            for gene, value in sorted(gene_gc.items()):
                outstr = f"{gene}\t{float(value['tot_depth'])/int(value['tot_region']):.4f}\t{value['cover1X']/value['len']:.4f}\t{value['cover10X']/value['len']:.4f}"
                out.write(outstr + "\n")


def target_bed(infile, target_gene_dict):
    # target_bed
    with open(infile, "r") as tar_bed:
        for lines in tar_bed:
            line = lines.strip().split("\t")
            if lines.startswith('#'):
                continue
            pos = tuple(line[:4])
            target_gene_dict[pos].append(line[4])


def target_depth(infile, target_gc, target_gene_dict,gene_gc, reads_gc):
    # target depth file
    with gzip.open(infile, 'rt') as tar_depth:
        tot_depth = 0
        tot_region = 0
        for lines in tar_depth:
            chrom, start, end, exon, depth = lines.strip().split("\t")
            pos = (chrom, start, end, exon)
            region_len = int(end) - int(start)
            target_gc[pos]['region_len'] = target_gc[pos].get('regeion_len', 0) + region_len
            target_gc[pos]['mean_depth'] = depth
            target_gc[pos]['tot_bases'] = int(float(depth) * region_len)
            tot_depth += float(depth) * int(region_len)
            tot_region += region_len
            if len(target_gene_dict) == 0:
                continue
            for gene in target_gene_dict[pos]:
                # gene = target_gene_dict[pos]
                gene_gc[gene]['tot_region'] = gene_gc[gene].get('tot_region', 0) + region_len
                gene_gc[gene]['tot_depth'] = gene_gc[gene].get('tot_depth', 0) + float(depth) * region_len
        reads_gc['平均深度'] = f"{tot_depth / tot_region:.2f}"

def target_cov(infile, target_gc, target_gene_dict, gene_gc, reads_gc):
    # target coverage file
    with gzip.open(infile, 'rt') as tar_cov:
        title_list = []
        cover_dict = defaultdict(dict)
        for i, lines in enumerate(tar_cov):
            chrom, start, end, exon, *cover_ = lines.strip().split("\t")
            if i == 0:
                title_list = cover_
                continue
            pos = (chrom, start, end, exon)
            region = int(end) - int(start)
            for j, cover in enumerate(title_list):
                cover_dict[cover]['cover'] = cover_dict[cover].get('cover', 0) + int(cover_[j])
                cover_dict[cover]['len'] = cover_dict[cover].get('len', 0) + int(region)
            target_gc[pos]['coverage_tot_len'] = cover_[0] + '/' + str(region)
            target_gc[pos]['coverage'] = f"{int(cover_[0]) / region:.4f}"
            if len(target_gene_dict) == 0:
                continue
            for gene in target_gene_dict[pos]:
            # gene = target_gene_dict[pos]
                gene_gc[gene]['len'] = gene_gc[gene].get('len', 0) + region
                gene_gc[gene]['cover1X'] = gene_gc[gene].get('cover1X', 0) + int(cover_[0])
                gene_gc[gene]['cover10X'] = gene_gc[gene].get('cover10X', 0) + int(cover_[4])
        for cc in cover_dict.keys():
            new_key = '覆盖度' + cc
            reads_gc[new_key] = f"{cover_dict[cc]['cover'] / int(cover_dict[cc]['len']):.4f}"


os.makedirs(args.outdir, exist_ok=True)
gene_gc = defaultdict(dict)
target_gene_dict = defaultdict(list)
target_gc = defaultdict(dict)
reads_gc = summary_control_quality(args.json, args.map_ratio, args.dup, args.Insert, args.capture)
if args.wgs_depth != '-' and args.wgs_cov != '-':
    sample_summary(args.wgs_depth, args.wgs_cov, reads_gc)

if args.gene_cov != '-' and args.gene_depth != '-':
    gene_gc = Gene_stat(args.gene_depth, args.gene_cov)

if args.target_bed != '-':
    target_bed(args.target_bed, target_gene_dict)

if args.target_depth != '-':
    target_depth(args.target_depth, target_gc, target_gene_dict, gene_gc, reads_gc)

if args.target_cov != '-':
    target_cov(args.target_cov, target_gc, target_gene_dict, gene_gc, reads_gc)

write_data(reads_gc, gene_gc, f"{args.outdir}/{args.prefix}")

# if args.target_bed != '-':
#     with open(f"{args.outdir}/{args.prefix}_target", 'w') as o:
#         titile = f"染色体\t区间起始\t区间终止\t外显子信息\t基因名\t区间总覆盖碱基数/区间总长\t区间平均深度\t覆盖长度/总长度\t覆盖度(1X)"
#         o.write(titile + '\n')
#         gene_stat_dict = defaultdict(dict)
#         for pos, stat in target_gc.items():
#             tot_bases_len = str(stat['tot_bases']) + '/' + str(stat['region_len'])
#             chrom_pos = '\t'.join(pos)
#             outstr = f"{chrom_pos}\t{target_gene_dict[pos]}\t{tot_bases_len}\t{stat['mean_depth']}\t{stat['coverage_tot_len']}\t{stat['coverage']}"
#             o.write(outstr + '\n')

if args.target_bed != '-':
    target_bedfile = args.target_bed
    output_targetfile = f"{args.outdir}/{args.prefix}_target"
    with open(target_bedfile) as inF, open(output_targetfile, 'w') as o:
        # titile = f"染色体\t区间起始\t区间终止\t外显子信息\t基因名\t区间总覆盖碱基数/区间总长\t区间平均深度\t覆盖长度/总长度\t覆盖度(1X)"
        # o.write(titile + '\n')

        # gene_stat_dict = defaultdict(dict)

        for line in inF:
        # for pos, stat in target_gc.items():
            line_parts = line.strip().split("\t")
            pos = tuple(line_parts[0:4])
            stat = target_gc[pos]

            tot_bases_len = str(stat['tot_bases']) + '/' + str(stat['region_len'])

            # chrom_pos = '\t'.join(pos)
            tmpline1 = line.strip()
            # outstr = f"{chrom_pos}\t{target_gene_dict[pos]}\t{tot_bases_len}\t{stat['mean_depth']}\t{stat['coverage_tot_len']}\t{stat['coverage']}"
            outstr = f"{tmpline1}\t{tot_bases_len}\t{stat['mean_depth']}\t{stat['coverage_tot_len']}\t{stat['coverage']}"
            o.write(outstr + '\n')



# # summary
# reads_gc = summary_control_quality(args.json, args.map_ratio, args.dup, args.Insert, args.capture)

# # wgs & gene
# if args.wgs_depth != '-' and args.wgs_cov != '-' and args.gene_cov != '-' and args.gene_depth != '-' and args.target_cov == '-' and args.target_bed == '-' and args.target_depth == '-':
#     sample_summary(args.wgs_depth, args.wgs_cov, reads_gc, 1)
#     gene_gc = Gene_stat(args.gene_depth, args.gene_cov)
#     write_data(reads_gc, gene_gc, f"{args.outdir}/{args.prefix}")
# # only wgs
# elif args.wgs_depth != '-' and args.wgs_cov != '-' and args.gene_cov == '-' and args.gene_depth == '-' and args.target_cov == '-' and args.target_bed == '-' and args.target_depth == '-':
#     sample_summary(args.wgs_depth, args.wgs_cov, reads_gc, 1)
#     write_data(reads_gc, 0, f"{args.outdir}/{args.prefix}")
# # only gene
# elif args.wgs_depth == '-' and args.wgs_cov == '-' and args.gene_cov != '-' and args.gene_depth != '-' and args.target_cov == '-' and args.target_bed == '-' and args.target_depth == '-':
#     sample_summary(args.gene_depth, args.gene_cov, reads_gc, 2)
#     gene_gc = Gene_stat(args.gene_depth, args.gene_cov)
#     write_data(reads_gc, gene_gc, f"{args.outdir}/{args.prefix}")
# # only target
# elif args.target_cov != '-' and args.target_bed != '-' and args.target_depth != '-' and args.wgs_depth == '-' and args.wgs_cov == '-' and args.gene_cov == '-' and args.gene_depth == '-':
#     Stat_target(args.target_bed, args.target_depth, args.target_cov, f"{args.outdir}/{args.prefix}", reads_gc)
# else:
#     print("请检查深度文件和覆盖度文件是否同时存在\n当计算TARGET统计结果时，需要额外参数target_bed\n")


# # print(time.time()-time_s)


