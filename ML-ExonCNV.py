# -*- coding: utf-8 -*-
# Author:: yangsh，heh，housy，yinzh
# Created Date: 2025-06-30
# Last Modified Date: 2025-06-30
# Version: 1.0
# Company: Chigene (Beijing) Translational Medical Research Center Co. Ltd.
# E-mail: yangshuanghao@zhuanhuayixue.org
# License: Apache License 2.0

import glob
from collections import defaultdict
import subprocess   
import os
import gzip
import logging
import fire
import shutil

#Execute tasks in parallel by a thread pool
def multiple_run(func,tasks,maxnum = 3):
    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
    from concurrent.futures import ThreadPoolExecutor
    with ThreadPoolExecutor(max_workers=maxnum) as executor:
        jobs = [executor.submit(func, *task) for task in tasks]
        for job in jobs:
            try:
                job.result()  # waiting mission complete
            except Exception as e:
                logging.error(f"An error occurred: {e}")

#clean intermediate files
def delete_tmpfiles(path_list):
    for path in path_list:
        if os.path.isfile(path):
            os.remove(path)
        elif os.path.isdir(path):
            shutil.rmtree(path)
        else:
            print(f"Path not found: {path}")

# fold80 requires both depth and stat to be fully processed
def getfold80(path):
    depthfile = f'{path}/depthcov/depthcov.mosdepth.region.dist.txt'
    f80cmd = f"grep total {depthfile} | awk '$3>=0.8' | head -1 | cut -f2"
    f80 = float(subprocess.run(f80cmd, shell=True, stdout=subprocess.PIPE, text=True).stdout.strip('\n'))
    statfile = f'{path}/Stat/stat_summary'
    with open(statfile, 'r') as fo:
        for x in fo:
            if x.startswith('raw'):
                continue
            linelist = x.strip('\n').split('\t')
            avgdepth = float(linelist[15])
    fold80 = round(avgdepth / f80, 4)
    with open(f'{path}/Stat/fold80.txt', 'w') as fw:
        print('FOLD80',file = fw)
        print(fold80, file=fw)

def make_handle(name,sex,outdir):
    samplepath = f'{outdir}/{name}'
    if not os.path.exists(samplepath):
        os.makedirs(samplepath)
    handlepath = f'{samplepath}/handle.txt'
    sexdir = {'XX':'FM','XY':'M'}
    with open(handlepath, 'w') as fw:
        print(f'{name}\tNAME\t{sexdir[sex]}\t{name}\t{name}\tNT\tproband\tNAME\tphenotype\t3\tn',file = fw)


def make_vcf(name,vcf,outdir):
    script_dir = os.path.dirname(os.path.abspath(__file__))
    vcfpath = f'{outdir}/{name}/Variants/variants.vcf.gz'
    if not os.path.exists(f'{outdir}/{name}/Variants'):
        os.makedirs(f'{outdir}/{name}/Variants')
    cmd1 = f'ln -sf {vcf} {vcfpath}'
    # print(f'---Running  {cmd1}')
    subprocess.run(cmd1, shell=True, text=True, stdout=subprocess.PIPE)
    cmd2 = f'python {script_dir}/scripts/vcf2table-1.3.py -S -vcf {vcf} -table {outdir}/{name}/Variants/variants.table_complete.txt'
    # print(f'---Running  {cmd2}')
    subprocess.run(cmd2, shell=True, text=True, stdout=subprocess.PIPE)
    cmd3 = f'cut -f1-5,7 {outdir}/{name}/Variants/variants.table_complete.txt > {outdir}/{name}/Variants/variants.table.txt'
    # print(f'---Running  {cmd3}')
    subprocess.run(cmd3, shell=True, text=True, stdout=subprocess.PIPE)

def get_dep_bed(name, bam, bed, ref, outdir,prename):
    output = f'{outdir}/{name}/{prename}/{prename}'
    if not os.path.exists(f'{outdir}/{name}/{prename}/'):
        os.makedirs(f'{outdir}/{name}/{prename}/')
    if not os.path.exists(f'{outdir}/{name}/Stat'):
        os.makedirs(f'{outdir}/{name}/Stat')
    mosdepth_cmd = f"mosdepth -x -n -t 3 -F 3844 -T 1,2,3,5,10,15,20,25,30,40,50 -b {bed} -f {ref} {output} {bam}"       #needs mosdepth
    # print(f'---Running  {mosdepth_cmd}')
    subprocess.run(mosdepth_cmd, shell=True, text=True, stdout=subprocess.PIPE)

#baseQC
def get_qcstat(bed,outdir):
    script_dir = os.path.dirname(os.path.abspath(__file__))
    # print(f"{outdir}/depthcov/depthcov.regions.bed.gz")
    cmd = f'python {script_dir}/scripts/Control_quality.py -target_depth {outdir}/depthcov/depthcov.regions.bed.gz -target_cov {outdir}/depthcov/depthcov.thresholds.bed.gz -target_bed {bed} -prefix stat -outdir {outdir}/Stat'     
    statsum,stattarget = f'{outdir}/Stat/stat_summary',f'{outdir}/Stat/stat_target'
    # print(f'---Running  {cmd}')
    subprocess.run(cmd, shell=True, text=True, stdout=subprocess.PIPE) ## hsy

#GC correction
def gc_correction(bed,outdir):
    script_dir = os.path.dirname(os.path.abspath(__file__))
    cmd = f'Rscript {script_dir}/scripts/gc-norm.R {outdir}/Stat/stat_target {bed} {outdir}/Stat/stat_target_gc'     #needs R
    # print(f'---Running  {cmd}')
    subprocess.run(cmd, shell=True, text=True, stdout=subprocess.PIPE)  ## hsy
    return f'{outdir}/Stat/stat_target_gc'

# Remove low-coverage or unstable coverage regions (i.e. low-quality regions)
def get_cleanbed(sampletable, outdir, sampcov, oribed, newbed):
    regiondic = defaultdict(int)
    samplenum = 0
    XXdic = defaultdict(int)
    XYdic = defaultdict(int)
    XXnum = 0
    XYnum = 0
    with open(sampletable, 'r') as fo:
        for x in fo:
            samplenum += 1
            line = x.strip('\n')
            linelist = line.split('\t')
            name, sex, bam,vcf = linelist
            if sex == 'XX':
                XXnum += 1
            elif sex == 'XY':
                XYnum += 1
            samplebed = glob.glob(f'{outdir}/{name}/rawdepthcov/rawdepthcov.regions.bed.gz')[0]
            with gzip.open(samplebed, 'rt') as fo:
                for x in fo:
                    line = x.strip('\n')
                    linelist = line.split('\t')
                    chrom, start, end = linelist[0:3]
                    dep = float(linelist[-1])
                    if chrom != 'chrX' and chrom != 'chrY':
                        regiondic[(chrom, start, end)] += 1
                    else:
                        if sex == 'XX':
                            XXdic[(chrom, start, end)] += 1
                        elif sex == 'XY':
                            XYdic[(chrom, start, end)] += 1
    with open(oribed, 'r') as fo, open(newbed, 'w') as fw:
        for x in fo:
            line = x.strip('\n')
            linelist = line.split('\t')
            chrom, start, end = linelist[0:3]
            if chrom != 'chrX' and chrom != 'chrY':
                if (chrom, start, end) in regiondic:
                    if regiondic[(chrom, start, end)] >= samplenum * float(sampcov):
                        print(line, file=fw)
                    else:
                        print(f'Ignoring {(chrom, start, end)} as not enough samples ({regiondic[(chrom, start, end)]}/{samplenum})')
            elif chrom == 'chrX':
                if (chrom, start, end) in XXdic or (chrom, start, end) in XYdic:
                    if XXdic[(chrom, start, end)] >= XXnum * float(sampcov) or XYdic[(chrom, start, end)] >= XYnum * float(sampcov):
                # if (chrom, start, end) in XXdic and (chrom, start, end) in XYdic:
                    # if XXdic[(chrom, start, end)] >= XXnum * float(sampcov) and XYdic[(chrom, start, end)] >= XYnum * float(sampcov):
                        print(line, file=fw)
                    else:
                        print(f'Ignoring {(chrom, start, end)} as not enough samples ({regiondic[(chrom, start, end)]}/{samplenum})')
            elif chrom == 'chrY':
                if (chrom, start, end) in XYdic:
                    if XYdic[(chrom, start, end)] >= XYnum * float(sampcov):
                        print(line, file=fw)
                    else:
                        print(f'Ignoring {(chrom, start, end)} as not enough samples ({regiondic[(chrom, start, end)]}/{samplenum})')

#depth info
def archivedb(samplepath,sex, database,tablename):
    script_dir = os.path.dirname(os.path.abspath(__file__))
    cmd = f'python {script_dir}/scripts/mk_deplib.py archivedb --samplepath {samplepath} --sex {sex} --database {database} --tablename {tablename}'
    # print(f'---Running  {cmd}')
    subprocess.run(cmd, shell=True, text=True, stdout=subprocess.PIPE)

# Build control library for each sample 
def getdeplib(samplepath,sex, control_save_path, database,tablename, control_num=30, corr_default = 0.9):
    script_dir = os.path.dirname(os.path.abspath(__file__))
    libcode = samplepath.split('/')[-1]
    cmd = f'python {script_dir}/scripts/mk_deplib_new.py getdeplib --samplepath {samplepath} --sex {sex} --control_save_path {control_save_path} --database {database} --tablename {tablename} --control_num {control_num} --corr_default {corr_default}'
    # print(f'---Running  {cmd}')
    subprocess.run(cmd, shell=True, text=True, stdout=subprocess.PIPE)

# Build frequency database
def mk_freqdb(dirfile,output):
    script_dir = os.path.dirname(os.path.abspath(__file__))
    cmd1 = f'python {script_dir}/src/exon_pipeline/apps/freq/build_freq_data.py build {dirfile} {output}raw'
    # print(f'---Running  {cmd1}')
    subprocess.run(cmd1, shell=True, text=True, stdout=subprocess.PIPE)
    cmd2 = f'python {script_dir}/src/exon_pipeline/tools/mod_wes_freq.py {output}raw {output}'
    # print(f'---Running  {cmd2}')
    subprocess.run(cmd2, shell=True, text=True, stdout=subprocess.PIPE)

def wescnv1(samplepath,cnvconfig,sex):
    script_dir = os.path.dirname(os.path.abspath(__file__))
    libcode = samplepath.split('/')[-1]
    statpath = f'{samplepath}/Stat'
    respath = f'{samplepath}/ExonResultV3'
    if not os.path.exists(respath):
        os.makedirs(respath)
    cmd0 = f'python {script_dir}/src/exon_pipeline/core/contrast.py --contrast_build_file {samplepath}/exon_lib1.txt --config-file {cnvconfig} --output {respath}/{sex}.txt '
    # print(f'---Running  {cmd0}')
    subprocess.run(cmd0, shell=True, text=True, stdout=subprocess.PIPE) ## hsy
    
    cmd1 = f'python {script_dir}/src/exon_pipeline/forward/make_exons.py --wkcode {libcode} --analyse_path {statpath} --contrast_file {respath}/{sex}_contrast.pic --config-file {cnvconfig} --sex {sex} --outdir {respath}'
    # print(f'---Running  {cmd1}')
    subprocess.run(cmd1, shell=True, text=True, stdout=subprocess.PIPE) 
    
    cmd2 = f'python {script_dir}/src/exon_pipeline/forward/make_flags.py --wkcode {libcode} --analyse_path {respath} --config-file {cnvconfig}' 
    # print(f'---Running  {cmd2}')
    subprocess.run(cmd2, shell=True, text=True, stdout=subprocess.PIPE)

def wescnv2(samplepath,cnvconfig,sex,freq_db,bam,version,tmpfile = False):
    script_dir = os.path.dirname(os.path.abspath(__file__))
    libcode = samplepath.split('/')[-1]
    respath = f'{samplepath}/ExonResultV3'
    if version == 'hg19':
        reference = f"{script_dir}/docs/hg19.fa"
    elif version == 'hg38':
        reference = f"{script_dir}/docs/hg38.fa"
    if not os.path.exists(respath):
        os.makedirs(respath)

    cmd3 = f'python {script_dir}/src/exon_pipeline/forward/get_freq_reliability.py --wkcode {libcode} --analyse_path {respath} --config-file {cnvconfig} --sex {sex} --freq_db_file {freq_db}'
    # print(f'---Running  {cmd3}')
    subprocess.run(cmd3, shell=True, text=True, stdout=subprocess.PIPE)

    cmd4 = f'python {script_dir}/src/exon_pipeline/forward/make_exons_details.py --wkcode {libcode} --analyse_path {respath} --config-file {cnvconfig} --outdir {respath}'
    # print(f'---Running  {cmd4}')
    subprocess.run(cmd4, shell=True, text=True, stdout=subprocess.PIPE)
    
    cmd5 = f'python -m exon_pipeline qc {libcode} {respath} {cnvconfig}'
    # print(f'---Running  {cmd5}')
    subprocess.run(cmd5, shell=True, text=True, stdout=subprocess.PIPE)
    
    cmd6 = f'python -m exon_pipeline qc_karyotype {libcode} {respath} {cnvconfig}'
    # print(f'---Running  {cmd6}')
    subprocess.run(cmd6, shell=True, text=True, stdout=subprocess.PIPE)
    
    cmd7 = f'python -m exon_pipeline format {libcode} {respath} {cnvconfig}'
    # print(f'---Running  {cmd7}')
    subprocess.run(cmd7, shell=True, text=True, stdout=subprocess.PIPE)
    
    #manta (python2)
    cmd8 = f'/opt/conda/envs/py27/bin/python2 {script_dir}/src/manta/manta-1.6.0//bin/configManta.py --bam {bam} --referenceFasta {reference} --runDir {samplepath}/Manta --exome'
    # print(f'---Running  {cmd8}')
    subprocess.run(cmd8, shell=True, text=True, stdout=subprocess.PIPE)
    cmd9 = f'/opt/conda/envs/py27/bin/python2 {samplepath}/Manta/runWorkflow.py -j 17 --quiet'
    # print(f'---Running  {cmd9}')
    subprocess.run(cmd9, shell=True, text=True, stdout=subprocess.PIPE)
    cmd10 = f'/opt/conda/envs/py27/bin/python2 {script_dir}/src/manta/convertInversion.py samtools {reference} {samplepath}/Manta/results/variants/diploidSV.vcf.gz > {samplepath}/Manta/diploidSV.vcf'
    # print(f'---Running  {cmd10}')
    subprocess.run(cmd10, shell=True, text=True, stdout=subprocess.PIPE)
    
    #mergemanta
    cmd11 = f'python {script_dir}/src/manta/merge-manta.py {samplepath}/Manta/diploidSV.vcf {samplepath}/ExonResultV3/exons_details.txt {samplepath}/ExonResultV3/qc_karyotype.txt {samplepath}/ExonResultV3/format_gene_info_pre.txt {samplepath}/ExonResultV3/format_gene_info-manta.txt --sam {libcode} --db_exon {script_dir}/docs/EXON.longestNM.bed --db_freq {script_dir}/src/manta/wesmanta.db'
    # print(f'---Running  {cmd11}')
    subprocess.run(cmd11, shell=True, text=True, stdout=subprocess.PIPE)
  
    #mos
    cmd12 = f'python {script_dir}/scripts/find-mos.py {samplepath}/ExonResultV3/exons_details.txt {samplepath}/ExonResultV3/mos.txt -v {samplepath}/Variants/variants.vcf.gz'
    # print(f'---Running  {cmd12}')
    subprocess.run(cmd12, shell=True, text=True, stdout=subprocess.PIPE)
    cmd13 = f'python {script_dir}/scripts/fix-sex.py {samplepath}/ExonResultV3/qc_karyotype.txt {samplepath}/ExonResultV3/mos.txt {samplepath}/ExonResultV3/qc_karyotype_mos.txt'
    # print(f'---Running  {cmd13}')
    subprocess.run(cmd13, shell=True, text=True, stdout=subprocess.PIPE)
    
    #reliability by ML
    cmd14 = f'/workspace/sklearn/bin/python {script_dir}/src/ExonReliability/finetunlossgainrevis.py notrio -input {samplepath}/ExonResultV3/format_gene_info-manta.txt -name {libcode} -exonmode 2 -svcf {samplepath}/Variants/variants.vcf.gz'
    # print(f'---Running  {cmd14}')
    subprocess.run(cmd14, shell=True, text=True, stdout=subprocess.PIPE)
    
    cmd15 = f'python {script_dir}/src/ExonReliability/Reportloci/reportlocibynewalgorithm.py {samplepath}/ExonResultV3/format_gene_info.txt {samplepath}/handle.txt {libcode} --filter_info=True'
    # print(f'---Running  {cmd15}')
    subprocess.run(cmd15, shell=True, text=True, stdout=subprocess.PIPE)
    
    #delete tmpfiles
    tmplist = [os.path.join(samplepath,'Variants'),os.path.join(samplepath,'exon_lib1.txt')]
    reslist = ['diploidSV.vcf','format_gene_info.txt','format_cnv_info_sample.txt','format_cnv_info_all_sample.txt',
'qc_info.txt','qc_karyotype_mos.txt','qc_karyotype.txt']
    exon_path = os.path.join(samplepath, "ExonResultV3", "*")
    manta_path = os.path.join(samplepath, "Manta", "*")
    for path in list(glob.glob(str(exon_path))) + list(glob.glob(str(manta_path))):
        filename = os.path.basename(path) 
        if filename not in reslist:
            tmplist.append(path)
    if not tmpfile:
        delete_tmpfiles(tmplist)

def stats_single(name,sex,bam,vcf, version,bed, outdir,tmpfile = False):       #NO filter (use newbed,depthcov)
    script_dir = os.path.dirname(os.path.abspath(__file__))     
    samplepath = f'{outdir}/{name}'
    if not os.path.exists(samplepath):
        os.makedirs(samplepath)
    make_handle(name,sex,outdir)
    make_vcf(name,vcf,outdir)
    if version == 'hg19':
        reference = f"{script_dir}/docs/hg19.fa"
    elif version == 'hg38':
        reference = f"{script_dir}/docs/hg38.fa"
    newbed = bed
    get_dep_bed(name, bam, newbed, reference, outdir,'depthcov')
    samplepath = f'{outdir}/{name}'
    get_qcstat(newbed,samplepath)
    getfold80(samplepath)
    gcbed = f'{newbed}_gc'
    nuc_command = f"bedtools nuc -fi {reference} -bed {newbed} >{gcbed}"
    subprocess.run(nuc_command, shell=True, capture_output=True, text=True, stdout=subprocess.PIPE)
    #GC correction
    depfile = f'{samplepath}/Stat/stat_target_gc'
    gc_correction(gcbed,samplepath)
    database = f"{script_dir}/docs/dep_contrast.db"
    depcontrast_tablename = f"wes_depcontrast_{version}"
    archivedb(samplepath,sex, database,depcontrast_tablename)
    tmplist = [os.path.join(samplepath,'depthcov')]
    if not tmpfile:
        delete_tmpfiles(tmplist)

def precnv_single(name,sex,bam,vcf,version,bed, outdir,tmpfile = False):
    stats_single(name,sex,bam,vcf, version,bed, outdir,tmpfile = tmpfile)
    cnvconfig = f'{outdir}/wescnv_final.ini'
    samplepath = f'{outdir}/{name}'
    depcontrastfile = f'{samplepath}/exon_lib1.txt'
    script_dir = os.path.dirname(os.path.abspath(__file__))     
    database = f"{script_dir}/docs/dep_contrast.db"
    depcontrast_tablename = f"wes_depcontrast_{version}"
    getdeplib(samplepath,sex, depcontrastfile, database,depcontrast_tablename, control_num=30, corr_default = 0.9)
    wescnv1(samplepath,cnvconfig,sex)

def single_analysis(name,sex,bam,vcf,bed,outdir,version,freqdb = '',tmpfile = False):
    precnv_single(name,sex,bam,vcf, version,bed, outdir,tmpfile = tmpfile)
    CNV_calling(name,sex,bam,outdir, version,freqdb = freqdb,tmpfile = tmpfile)

# Module 1: Obtain basic stats required for sample analysis (output: sample/Stat/stat_target_gc)
def get_stats(inputfile,outdir,version,sampcov = 0.8,max_threads = 3,min_dep = 50,tmpfile = False):
    script_dir = os.path.dirname(os.path.abspath(__file__))    
    if version == 'hg19':
        reference = f"{script_dir}/docs/hg19.fa"
        regionfile = f"{script_dir}/docs/target_functionnm.bed"
    elif version == 'hg38':
        reference = f"{script_dir}/docs/hg38.fa"
        regionfile = f"{script_dir}/docs/hg38_target.bed"
    # get paths
    cmd = f"cp {script_dir}/docs/dep_contrast.db {outdir}/dep_contrast.db && cp {script_dir}/docs/Rngc.bed {outdir}/Rngc.bed"
    subprocess.run(cmd, shell=True, text=True, stdout=subprocess.PIPE)
    database = f"{outdir}/dep_contrast.db"
    # database = f"{script_dir}/docs/dep_contrast.db"
    rawcnvconfig = f"{script_dir}/docs/wescnv.ini"
    nt_target_bed = f"{script_dir}/docs/target_functionnm.bed"
    telomere_bed = f"{script_dir}/docs/Tolomere1.bed"
    gene_att = f"{script_dir}/docs/gene_attentions.txt"
    gene_att1 = f"{script_dir}/docs/gene_attentions.txt_1"
    gene_att2 = f"{script_dir}/docs/gene_attentions.txt_2"
    low_qua_bed = f"{script_dir}/docs/LowQualityRegion.bed"
    repeat_obj_file = f"{script_dir}/docs/exons_max_repeat_new.pic"
    unbalance_exons_file = f"{script_dir}/docs/D-unbal-reg.txt"
    gene_filterfalse_file = f"{script_dir}/docs/gene_filterfalse.txt"
    vaf_db = f"{script_dir}/docs/s50_test.pic"
   
    inputname = inputfile.split('/')[-1].split('.')[0]
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    #generate config
    cnvconfig = f'{outdir}/wescnv_final.ini'
    with open(rawcnvconfig, 'r') as fo, open(cnvconfig,'w') as fw:
        for x in fo:
            line = x.strip('\n')
            if not line:
                continue
            print(line,file = fw)
        print('[base]',file = fw)
        print(f'nt_target = {nt_target_bed}',file = fw)
        print(f'tolomere_file = {telomere_bed}',file = fw)
        print(f'gene_attentions_file = {gene_att}',file = fw)
        print(f'gene_attentions_file1 = {gene_att1}',file = fw)
        print(f'gene_attentions_file2 = {gene_att2}',file = fw)
        print(f'low_coverage_file = {low_qua_bed}',file = fw)
        print(f'repeat_obj_file = {repeat_obj_file}',file = fw)
        print(f'unbalance_exons_file = {unbalance_exons_file}',file = fw)
        print(f'gene_filterfalse_file = {gene_filterfalse_file}',file = fw)
        print(f'vaf_db = {vaf_db}',file = fw)
        
    samplelist = []
    tmplist = []
    mosdep_missionlist1 = []
    mosdep_missionlist2 = []
    vcf_missionlist = []
    handle_missionlist = []
    stats_missionlist = []
    newbed = f"{outdir}/{inputname}_{version}.bed"
    
    with open(inputfile, 'r') as fo:
        for x in fo:
            line = x.strip('\n')
            linelist = line.split('\t')
            name, sex, bam, vcf = linelist
            samplelist.append((name,sex,bam,vcf))
            vcf_missionlist.append((name,vcf,outdir))
            handle_missionlist.append((name,sex,outdir))
            mosdep_missionlist1.append((name, bam, regionfile, reference, outdir,'rawdepthcov'))
            mosdep_missionlist2.append((name, bam, newbed, reference, outdir,'depthcov'))
            samplepath = f'{outdir}/{name}'
            stats_missionlist.append((newbed,samplepath))
            tmplist.append(os.path.join(outdir, name,'rawdepthcov'))
            tmplist.append(os.path.join(outdir, name,'depthcov'))
    # print(mosdep_missionlist1)
    multiple_run(make_handle,handle_missionlist,maxnum = int(max_threads))
    multiple_run(make_vcf,vcf_missionlist,maxnum = int(max_threads))
    # #rawmosdep
    multiple_run(get_dep_bed,mosdep_missionlist1,maxnum = int(max_threads))  ## hsy
    # #filter out low-quality region
    get_cleanbed(inputfile, outdir, sampcov, regionfile, newbed)
    # print(f'newbed: {newbed}')
    gcbed = f'{newbed}_gc'
    nuc_command = f"bedtools nuc -fi {reference} -bed {newbed} >{gcbed}"
    subprocess.run(nuc_command, shell=True, capture_output=True, text=True, stdout=subprocess.PIPE)  ## hsy
    
    #mosdepth for new bed
    # multiple_run(get_dep_bed,mosdep_missionlist2,maxnum = int(max_threads))
    
    #get final stats
    multiple_run(get_qcstat,stats_missionlist,maxnum = int(max_threads))
    sampleqcdic = {}
    sampledepdic = {}
    ignoredlist = []
    for (name,sex,bam,vcf) in samplelist:
        samplepath = f'{outdir}/{name}'
        statsum = f'{samplepath}/Stat/stat_summary'
        stattarget = f'{samplepath}/Stat/stat_target'
        sampleqcdic[(name,sex,bam,vcf)] = statsum
        sampledepdic[(name,sex,bam,vcf)] = stattarget
        getfold80(samplepath)
        #filter out samples with average depth less than min_dep
        with open(sampleqcdic[(name,sex,bam,vcf)], 'r') as fo:
            for x in fo:
                if x.startswith('raw'):
                    continue
                line = x.strip('\n')
                linelist = line.split('\t')
                if float(linelist[15]) < min_dep:
                    # print(f'Ignoring {name} {sex} {bam} which has average depth  {linelist[15]} in {newbed}')
                    ignoredlist.append((name,sex,bam,vcf))
                    continue
    
    #filter out low-quality samples
    filtered_samples = [item for item in samplelist if item not in ignoredlist]
    filtered_sample_file = f'{outdir}/sample_used_{inputname}.txt'
    XXnum = 0
    XYnum = 0
    gc_missionlist = []
    with open(inputfile, 'r') as fo, open(filtered_sample_file, 'w') as fw:
        for x in fo:
            line = x.strip('\n')
            linelist = line.split('\t')
            name, sex, bam,vcf = linelist[0:4]
            if (name,sex,bam,vcf) in filtered_samples:
                print(line, file=fw)
                if sex == 'XX':
                    XXnum += 1
                else:
                    XYnum += 1
            samplepath = f'{outdir}/{name}'
            #GC correction
            depfile = f'{samplepath}/Stat/stat_target_gc'
            gc_missionlist.append((gcbed,samplepath))
    
    print(f'Used {len(filtered_samples)} (XX:{XXnum};XY:{XYnum}) from {len(samplelist)} samples')
    
    #Run GC correction
    multiple_run(gc_correction,gc_missionlist,maxnum = int(max_threads ))
    
    # record depth info
    depcontrast_tablename = f"wes_depcontrast_{version}"
    archmission_list = []
    for (name,sex,bam,vcf) in filtered_samples:
        samplepath = f'{outdir}/{name}'
        archmission_list.append((samplepath,sex, database,depcontrast_tablename))
    multiple_run(archivedb,archmission_list,maxnum = int(max_threads))
    #delete tmpfiles
    if not tmpfile:
        delete_tmpfiles(tmplist)
    return filtered_sample_file

#Module2: build control database (based on sample stats)
def get_cnv_db(inputfile,outdir,version,filtered_sample_file,max_threads = 3,tmpfile = False):
    script_dir = os.path.dirname(os.path.abspath(__file__)) 
    database = f"{outdir}/dep_contrast.db"
    inputname = inputfile.split('/')[-1].split('.')[0]
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    cnvconfig = f'{outdir}/wescnv_final.ini'
    depcontrast_tablename = f"wes_depcontrast_{version}"
    XXdirlist = []
    XYdirlist = []
    getdepdb_missionlist = []
    wescnv_missonlist1 = []
    filtered_samples = []
    tmplist = []
    with open(filtered_sample_file,'r') as fo:
        for x in fo:
            line = x.strip('\n')
            if not line:
                continue
            name,sex,bam,vcf = line.split('\t')
            filtered_samples.append((name,sex,bam,vcf))
    for (name,sex,bam,vcf) in filtered_samples:
        samplepath = f'{outdir}/{name}'
        depcontrastfile = f'{samplepath}/exon_lib1.txt'
        getdepdb_missionlist.append((samplepath,sex, depcontrastfile, database,depcontrast_tablename,30, 0.9))
        wescnv_missonlist1.append((samplepath,cnvconfig,sex))
        dirinfo = f'{name}\t{sex}\t{samplepath}/ExonResultV3'
        if sex == 'XX':
            XXdirlist.append(dirinfo)
        else:
            XYdirlist.append(dirinfo)
    
    multiple_run(getdeplib,getdepdb_missionlist,maxnum = int(max_threads))
    multiple_run(wescnv1,wescnv_missonlist1,maxnum = int(max_threads))

    XXdirfile = f'{outdir}/XX_dir_{inputname}.txt'
    with open(f'{outdir}/XX_dir_{inputname}.txt', 'w') as fw:
        for dir in XXdirlist:
            print(dir, file=fw)
    XYdirfile = f'{outdir}/XY_dir_{inputname}.txt'
    with open(f'{outdir}/XY_dir_{inputname}.txt', 'w') as fw:
        for dir in XYdirlist:
            print(dir, file=fw)
    
    XXfreqdb = f'{outdir}/XX_freqdb.pic'
    mk_freqdb(XXdirfile,XXfreqdb)
    XYfreqdb = f'{outdir}/XY_freqdb.pic'
    mk_freqdb(XYdirfile,XYfreqdb)
    tmplist = [depcontrastfile,XXdirfile,XYdirfile,f'{XXfreqdb}raw',f'{XYfreqdb}raw',]
    #delete tmpfiles
    if not tmpfile:
        delete_tmpfiles(tmplist)

#Module3: get wescnv result (based on control database)
#Single sample
def CNV_detect(name,sex,bam,outdir, version,freqdb = '',tmpfile = False):
    samplepath = f'{outdir}/{name}'
    cnvconfig = f'{outdir}/wescnv_final.ini'
    if not freqdb:
        freqdb = f'{outdir}/{sex}_freqdb.pic'
    if not glob.glob(freqdb):
        assert False,f'frequency database {freqdb} does not exist'
    wescnv2(samplepath,cnvconfig,sex,freqdb,bam,version,tmpfile = tmpfile)


def CNV_calling(inputfile,outdir,version,freqdb = '',tmpfile = False, regions = '',max_threads = 5):
    script_dir = os.path.dirname(os.path.abspath(__file__))
    cnvmissionlist = []
    companamissionlist = []
    with open(inputfile, 'r') as fo:
        for x in fo:
            line = x.strip('\n')
            linelist = line.split('\t')
            name, sex, bam, vcf = linelist
            samplepath = f'{outdir}/{name}'
            requestfilelist = ['exon_lib1.txt','ExonResultV3','Stat','Variants']
            stats = True
            for request in requestfilelist:
                if not glob.glob(f'{samplepath}/{request}'):
                    stats = False
            if stats:
                cnvmissionlist.append((name,sex,bam,outdir, version,freqdb,tmpfile))
                # CNV_detect(name,sex,bam,outdir, version,freqdb,tmpfile = tmpfile)
            else:
                if not regions:
                    inputname = inputfile.split('/')[-1].split('.')[0]
                    newbed = f"{outdir}/{inputname}_{version}.bed"
                    if glob.glob(newbed):
                        regions = newbed
                    else:
                        regions = f"{script_dir}/docs/{version}_target.bed"
                companamissionlist.append((name,sex,bam,vcf,regions,outdir,version,freqdb,tmpfile))
                # single_analysis(name,sex,bam,vcf,regions,outdir,version,freqdb = freqdb,tmpfile = tmpfile)
    if cnvmissionlist:
        multiple_run(CNV_detect,cnvmissionlist,maxnum = int(max_threads))
    if companamissionlist:
        multiple_run(single_analysis,companamissionlist,maxnum = int(max_threads))
        

#Get control database for cnv calling (combine Module1 and Module2)
def control_train(inputfile,outdir,version,sampcov = 0.8,max_threads = 3,min_dep = 50,tmpfile = False):
    filtered_sample_file = get_stats(inputfile,outdir,version,sampcov = sampcov,max_threads = max_threads,min_dep = 50,tmpfile = tmpfile)
    get_cnv_db(inputfile,outdir,version,filtered_sample_file,max_threads = max_threads,tmpfile = tmpfile)
def main(hidden=True):
    functions = {
        "control-build": control_train,
        "CNV_calling": CNV_calling,
        "single_analysis": single_analysis,
    }

    if not hidden:
        functions.update({
            "get_stats": get_stats,
            "get_cnv_db": get_cnv_db,
            "precnv_single": precnv_single,
            "stats_single": stats_single
        })
    fire.Fire(functions)

if __name__ == "__main__":
    main()
    # main(hidden = False)
