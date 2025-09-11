<!--- TOP OF README.MD --->
<div align="center">

# ⚠️  NON-COMMERCIAL LICENSE  ⚠️

**This repository is licensed under *Creative Commons Attribution-NonCommercial 4.0 International* (CC BY-NC 4.0).**  
You may freely use, share and modify the material **for any non-commercial purpose** provided you give appropriate credit.

**ANY COMMERCIAL USE, DISTRIBUTION, OR SAAS DEPLOYMENT IS *PROHIBITED* WITHOUT A SEPARATE WRITTEN AGREEMENT.**  
Unauthorized commercial usage violates copyright and patent rights and will be pursued in court.

Contact `yshbioinfo@126.com` for commercial licensing options.

</div>

---


# ML-ExonCNV  
*A machine learning-based tool for exon-level Copy Number Variation (CNV) detection from DNA sequencing data.*  



## Table of Contents
- [Disclaimer](#disclaimer)    
- [Installation](#installation)  
- [Input/Output](#inputoutput)  
- [Usage](#usage)  


## Disclaimer
This software is made available for research purposes only. The developers make no warranties regarding the correctness, performance, or suitability of the software for any purpose. Use of the software is entirely at the user’s own risk, and the authors shall not be held liable for any direct or indirect damages arising from its use.


## Installation  
#### Docker 
```bash
docker pull chigenelaby/mlexoncnv:1.0.0

docker run -it --rm --name containter_name \
-v /absolute/host/path/inputdir:/absolute/container/path/inputdir \
-v /absolute/host/path/outputdir:/absolute/container/path/outputdir \
chigenelaby/mlexoncnv:1.0.0 \
python /workspace/ML-exon/ML-ExonCNV.py control-build INPUTFILE OUTDIR VERSION REFERENCE <flags> 


docker run -it --rm --name containter_name \
-v /absolute/host/path/inputdir:/absolute/container/path/inputdir \
-v /absolute/host/path/outputdir:/absolute/container/path/outputdir \
chigenelaby/mlexoncnv:1.0.0 \
python /workspace/ML-exon/ML-ExonCNV.py CNV_calling INPUTFILE OUTDIR VERSION REFERENCE <flags> 

```

## Input/Output
### Input
Sample information file (TSV)​ with the following form which indicates paths of bam files and vcf.gz files of samples with their names and genders (XX or XY). 

| samplename | XX/XY | /absolute/container/path/inputdir/bamfile | /absolute/container/path/inputdir/vcffile |
| ---------- | ----- | ----------------------------------------- | ----------------------------------------- |
| ALP294     | XX    | /data/ml_exon_data/ALP294.bam             | /data/ml_exon_data/ALP294.variants.vcf.gz |



### Output
For each processed sample, the following output files are generated in the sample's directory:

`ExonResultV3/format_gene_info.txt`: the processed exon-level CNVs results. 

`ExonResultV3/format_cnv_info_all_sample.txt`: the processed wes-CNVs results.

`ExonResultV3/mos.txt`: the processed mosaic CNVs results.

`ExonResultV3/qc_info.txt`: technical variability 

`ExonResultV3/qc_karyotype.txt`: karyotype results

`ExonResultV3/qc_karyotype_mos.txt`: mosaic karyotype results


#### File Descriptions

- `ExonResultV3/format_gene_info.txt` includes the following columns:


| #chr |   start   |    end    | gene_name | freq | all_freq |             infos             | exon_cnv_foldchange | exon_cnv_type | exon_cnv_reliability |
| :--: | :-------: | :-------: | :-------: | :--: | :------: | :---------------------------: | :-----------------: | :-----------: | :------------------: |
| chr1 | 43424304  | 43424322  |  SLC2A1   | 0.0  |   0.0    | {"report_tag": "Report_loci"} |       1.4927        |     gain1     | relatively_reliable  |
| chr1 | 169565363 | 169572449 |   SELP    | 0.05 |   0.05   |     {"report_tag": "NA"}      |       0.6967        |     loss1     |      unreliable      |

chr: Chromosome

start: Start position 

end: End position

gene_name: Gene name

freq: Frequency of specific exon-level CNVs

all_freq: Frequency of all types of exon-level CNVs

infos: Whether the exon-level CNV is reported. Report loci identified as reliable and/or low-frequency.

exon_cnv_foldchange: Fold change of exon-level CNVs

exon_cnv_type: Type of exon-level CNVs

exon_cnv_reliability: Reliability of exon-level CNVs



- `ExonResultV3/format_cnv_info_all_sample.txt` includes the following columns:


| chr  |   start   |    end    | foldchange | cnv_type  |
| :--: | :-------: | :-------: | :--------: | :---: |
| chr1 | 147599111 | 148202568 |   1.5187   | gain1 |

chr: Chromosome

start: Start position 

end: End position

foldchange: Fold change of wes-CNVs

cnv_type: Type of wes-CNVs


- `ExonResultV3/mos.txt` includes the following columns:

| chr  |   start   |    end    | mos_type | foldchange  | prt_cell  |
| :--: | :-------: | :-------: | :--------: | :---: | :---: |
| chr1 | 147599111 | 148202568 |   gainmos   | 1.24 | 0.35 |

chr: Chromosome

start: Start position

end: End position

mos_type: Type of mosaic CNVs

foldchange: Fold change of mosaic CNVs

prt_cell: Chimerism Rate​


## Usage
```
NAME
    ML-ExonCNV.py

SYNOPSIS
    ML-ExonCNV.py COMMAND

COMMANDS
    COMMAND is one of the following:

     control-build

     CNV_calling

     single_analysis
```

### control-build

```
ML-ExonCNV.py control-build INPUTFILE OUTDIR VERSION REFERENCE <flags>
  optional flags:        --sampcov | --max_threads | --min_dep
```
*INPUTFILE*:  Sample information file(TSV) as mentioned before(one sample per line, minimum 11 same-sex samples required for the performance of algorithm)

*OUTDIR*:  Output directory  

*VERSION*: Reference genome version (hg19)

*REFERENCE*: Reference genome

Optional:

*--sampcov*: X (float, 0-1): Retain the interval if ≥ < X > % (default: 0.8) of samples have a minimum depth (*--min_dep*, default: 50) in this region.

*--min_dep*:  Minimum average depth (default: 50)  of samples in the genomic intervals defined as above.

*--max_threads*: Maximum parallel threads (default: 3)  

Main Outcomes:
1.  **Processed sample list** which includes all samples used with the original input file format.
2.  **Filtered target regions**  
    
3.  **Gender-specific frequency databases**
  
       -  XX_freqdb.pic (female reference)
          
       -  XY_freqdb.pic (male reference)
    
4.  **Per-sample output directories**
    -   Named according to sample names from input file
        
    -   Each contains sample-level analysis

### CNV_calling
```
ML-ExonCNV.py CNV_calling INPUTFILE OUTDIR VERSION REFERENCE <flags>
  optional flags:        --freqdb | --regions | --max_threads
```
*INPUTFILE*, *OUTDIR*, *VERSION* , and *REFERENCE*: same as in *control-build* step

*--freqdb*: accepts a custom frequency database file. If unspecified, the tool defaults to the database generated during the *control-build* step

*--regions*: accepts a custom target regions file. If unspecified, the tool defaults to the *filtered target regions*  generated during the *control-build* step

This step will produce the final results (format_gene_info.txt, format_cnv_info_all_sample.txt, etc.) of the samples given in INPUTFILE.

### single_analysis
```
ML-ExonCNV.py single_analysis NAME SEX BAM VCF BED OUTDIR VERSION REFERENCE <flags>
  optional flags:        --freqdb
```
The inputs here are basically the same as above, name/sex/bam/vcf is the infomation given by a row from the input file, bed is the *filtered target file* from *control-build* or any bed you want, and *--freqdb* is defaulted by the one generated under outdir by *control-build*.
