__author__ = 'yzh'
__version__ = 1.0
__date__ = 2024.2
import os
import sys
import fire
import glob
import math
import sqlite3
import subprocess
import pandas as pd
import re
from sqlalchemy import create_engine
#对照库归档制作
def archivedb(samplepath,sex, database,tablename):
    # 指定数据库
    conn = sqlite3.connect(f"{database}")
    cursor = conn.cursor()
    # 创建user表:表中的字段  "文库号，批次号，芯片号，送检性别，分析性别，表型，路径"
    cursor.execute(f'''
    CREATE TABLE IF NOT EXISTS {tablename} (
        ID TEXT UNIQUE,
        sex VARCHAR(2), 
        save_path TEXT
    )
    ''')
    libcode = samplepath.split("/")[-1]
    if os.path.exists(f"{samplepath}/Stat/stat_target_gc"):
        save_path_ = f"{samplepath}/Stat"
        save_info = (libcode,sex,save_path_)
        cursor.execute(f'''
        INSERT OR REPLACE INTO {tablename} (ID,sex,save_path) VALUES (?,?,?)''',save_info)
        conn.commit()
        cursor.close()
        conn.close()
    else:
        print(f"核查{libcode}的stat_target_gc是否存在")
        sys.exit()
    
    
#计算相关性系数做每个样本自己的对照库
def getdeplib(samplepath,sex, control_save_path, database,tablename, control_num=25, corr_default = 0.9):
    query_sql = f'SELECT * FROM {tablename}'
    engine = create_engine(f'sqlite:///{database}')
    control_lib_data = pd.read_sql_query(query_sql,engine)
    control_num = int(control_num) 
    save_info_ = dict()
    low_default_threshold = dict()
    # 移除输入文库号计算相关系数
    random_samples = None
    libcode = samplepath.split("/")[-1]
    control_lib_data = control_lib_data[(control_lib_data["ID"]!=libcode) & (control_lib_data['sex']==sex)]
    # contral_lib_data = contral_lib_data[(contral_lib_data["library"]!=libcode) & (contral_lib_data['analyse_sex']==sex)]
    if control_lib_data.shape[0]<=2000:
        random_samples = control_lib_data["save_path"].sample(n=control_lib_data.shape[0])
    else:
        random_samples = control_lib_data["save_path"].sample(n=2000)
    if random_samples is not None:
        random_list = random_samples.tolist()
        if not os.path.exists(f"{samplepath}/Stat/stat_target_gc"):
            print(f"文库号{libcode}需要的stat_target_gc文件不存在")
            sys.exit()
        ori_depfile = pd.read_csv(f"{samplepath}/Stat/stat_target_gc", sep="\t",
                            names=["chr", "start", "end", "exon", "gene", "E", "F", "G", "H"])
        # 按染色体分组并随机抽样
        # 抽样时只保留前三列位置信息
        ori_depfile_regions = ori_depfile[["chr", "start", "end"]]
        grouped_id = ori_depfile_regions.groupby("chr")
        selected_sampleregions = pd.DataFrame(columns=ori_depfile_regions.columns)
        #对每条染色体的区间抽样 汇总bed
        for (group_name, group_data) in grouped_id:
            random_samples = group_data.sample(n=50, random_state=42)
            selected_sampleregions = pd.concat([selected_sampleregions, random_samples])
        #用上面的区间抽样结果获取原文件dep数据
        original_dep_sample = pd.merge(ori_depfile, selected_sampleregions, on=["chr", "start", "end"])
        #从数据库中挑选对照组
        for path_rmdup_file in random_list:
            if not os.path.exists(path_rmdup_file+"/stat_target_gc"):
                continue
            file_candidate = pd.read_csv(path_rmdup_file+"/stat_target_gc", sep="\t",
                                            names=["chr", "start", "end", "exon", "gene", "E", "F", "G", "H"])
            try:
                candidate_sampleregions = pd.merge(file_candidate, selected_sampleregions, on=["chr", "start", "end"])
            except ValueError as e:
                #合并出问题的文件
                continue
            save_libcode = path_rmdup_file.split('/')[-2].split("_")[0]
            #用样本数据的区间与原数据的区间进行相关系数计算
            corr_ = candidate_sampleregions["F"].corr(original_dep_sample["F"])
            # print(save_libcode,corr_)
            # if float(corr_) > 80:
            #     print("candidate_sampleregions:")
            #     print(candidate_sampleregions.head())
            #     print("\nori_depfile:")
            #     print(original_dep_sample.head())
            
            if math.isnan(corr_):  # 排除字典中的nan值
                continue
            
            if save_libcode == libcode:
                continue
            low_default_threshold[f"{save_libcode}@@@{path_rmdup_file}"] = corr_
            if corr_ < corr_default:
                continue
            save_info_[f"{save_libcode}@@@{path_rmdup_file}"] = corr_
            # print(len(save_info_),save_info_)
            if len(save_info_) == control_num:
                with open(f"{control_save_path}", "w") as f1:
                    for _ in save_info_.keys():
                        f1.write("\t".join(_.split("@@@")) + "\n")
                break
        if 10 <= len(save_info_)<control_num:
            use_control = save_info_
            with open(f"{control_save_path}", "w") as f1:
                for _ in use_control.keys():
                    f1.write("\t".join(_.split("@@@")) + "\n")
        elif len(save_info_)<10:
            if len(low_default_threshold) >=10:
                low_default_threshold_sort = sorted(low_default_threshold.items(), key=lambda t: t[1], reverse=True)
                with open(f"{control_save_path}", "w") as f1:
                    for key_ in low_default_threshold_sort[:11]:
                        f1.write("\t".join(key_[0].split("@@@")) + "\n")
            else:
                print(f"计算相关系数选取的对照数目小于10个({len(save_info_)})退出")
                sys.exit()
    else:
        print("计算相关系数，移除输入文库号的数据库空!")
        sys.exit()

if __name__ == '__main__':
    fire.Fire()
