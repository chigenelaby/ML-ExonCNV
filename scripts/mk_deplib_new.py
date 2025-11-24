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
from multiprocessing import Pool, Manager, cpu_count
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


def process_file(opts):
    path_rmdup_file, selected_samples, merged_data_original = opts
    save_library_num = path_rmdup_file.split('/')[-2].split("_")[0]
    if not os.path.exists(path_rmdup_file + "/stat_target"):
        print(path_rmdup_file, "无文件")
        return f"{save_library_num}@@@{path_rmdup_file}", -3
#    save_library_num = path_rmdup_file.split('/')[-2].split("_")[0] 
    try:
        file_need_contrala = pd.read_csv(path_rmdup_file + "/stat_target", sep="\t",
                                       names=["chr", "start", "end", "exon", "gene", "E", "F", "G", "H"],usecols=["chr", "start", "end", "E", "F", "G", "H"])
        file_need_contral = file_need_contrala.drop_duplicates(subset=["chr", "start", "end", "E", "F", "G", "H"], keep='first')
        merged_data_find_contral = pd.merge(file_need_contral, selected_samples, on=["chr", "start", "end"])

    except ValueError as e:
        print(path_rmdup_file, "合并问题")
        #exit()
        return f"{save_library_num}@@@{path_rmdup_file}", -1
    merged = pd.merge( merged_data_find_contral, merged_data_original,
                       on=["chr", "start", "end"], #, "gene"],  # 根据实际对齐需求选择关键列
                       suffixes=("_find", "_original")      # 避免列名冲突
                      ) 
#    # 筛选两列 F 均 >10 的行
    filtered = merged
    filtered['F_find'] = pd.to_numeric(filtered['F_find'], errors='coerce')
    filtered['F_original'] = pd.to_numeric(filtered['F_original'], errors='coerce')
    filtered = filtered.dropna(subset=['F_find', 'F_original'])
    try:
        corr_ = filtered["F_find"].corr(filtered["F_original"])
##    corr_ = merged_data_find_contral["F"].corr(merged_data_original["F"])
    #    print(path_rmdup_file, corr_)
    
        if math.isnan(corr_):
           return f"{save_library_num}@@@{path_rmdup_file}", 0
        else:
           return f"{save_library_num}@@@{path_rmdup_file}", corr_
    except ValueError as e:
        print(path_rmdup_file, "合并数据问题")
        return f"{save_library_num}@@@{path_rmdup_file}", -1

    
#计算相关性系数做每个样本自己的对照库
def getdeplib(samplepath,sex, control_save_path, database,tablename, control_num=30, corr_default = 0.9):
    query_sql = f'SELECT * FROM {tablename}'
    engine = create_engine(f'sqlite:///{database}')
    union_field = os.path.join(os.path.dirname(database), "Rngc.bed")
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
    # print(random_samples)
    if random_samples is not None:
        random_list = random_samples.tolist()
        if not os.path.exists(f"{samplepath}/Stat/stat_target"):
            print(f"文库号{libcode}需要的stat_target_gc文件不存在")
            sys.exit()
        file_originala = pd.read_csv(f"{samplepath}/Stat/stat_target", sep="\t",
                                        names=["chr", "start", "end", "exon", "gene", "E", "F", "G", "H"], usecols=["chr", "start", "end", "E", "F", "G", "H"])

        file_original = file_originala.drop_duplicates(subset=["chr", "start", "end", "E", "F", "G", "H"], keep='first')
        file_id = pd.read_csv(f"{union_field}", sep="\t",
                                names=["chr", "start", "end", "GC"])  # 最小共有区域
        bins = [i/10 for i in range(0, 11)]
        labels = [f"{i/10}-{(i+1)/10}" for i in range(10)]
        file_id['GC_bin'] = pd.cut(file_id['GC'], bins=bins, labels=labels, include_lowest=True)
        #sample_size = min(1000, len(group))
        #randindex=group.sample(n=sample_size, random_state=42)
        grouped_id = file_id.groupby("GC_bin")
        selected_samples = pd.DataFrame(columns=file_id.columns)
        for group_name, group_data in grouped_id:
                sample_size = min(1000, int(len(group_data)*0.9))
                if '0.4' in group_name or '0.5' in group_name or '0.6' in group_name:
                    sample_size = min(500, int(len(group_data)*0.2))
                elif sample_size<500:
                    sample_size = len(group_data)
                # print(group_name,sample_size)
                random_samples = group_data.sample(n=sample_size, random_state=42)
                selected_samples = pd.concat([selected_samples, random_samples])
        merged_data_original = pd.merge(file_original, selected_samples[["chr", "start", "end"]], on=["chr", "start", "end"])

        # print(len(random_list),len(merged_data_original))
        #exit()
        with Manager() as manager:
            low_default_threshold = manager.dict()
            save_info_ = manager.dict()
            with Pool(15) as pool:
                    results = pool.map(process_file, [(path, selected_samples[["chr", "start", "end"]],merged_data_original) for path in random_list])
    #                 exit()
                    pool.close()
                    pool.join()
            max_corr=sorted(results, key=lambda x:x[1],reverse=True)[0]
            # print(len(results))
            m1, m2, m3=0, 0, 0
            for result in results:
                    if result is None:
                        continue
                    key, corr = result

                    low_default_threshold[key] = corr
                    if corr < corr_default: #or corr < max_corr[1]-0.015:
                        continue
                    save_info_[key] = corr

            # print(len(save_info_), m1, m2, m3)
            if len(save_info_) >= control_num:
                sorted_dict = dict(sorted(save_info_.items(), key=lambda item: item[1], reverse=True)[:control_num])
                with open(f"{control_save_path}", "w") as f1:
                    for _ in sorted_dict.keys():
                        # print(_,save_info_[_])
                        f1.write("\t".join(_.split("@@@")) + "\n")
            elif 10 <=len(save_info_)<control_num:
                use_contral = save_info_
                with open(f"{control_save_path}", "w") as f1:
                    for _ in use_contral.keys():
                        # print(_,save_info_[_])
                        f1.write("\t".join(_.split("@@@")) + "\n")
            elif len(save_info_)<10:
                if len(low_default_threshold) >11:
                    low_default_threshold_sort = sorted(low_default_threshold.items(), key=lambda t: t[1], reverse=True)
                    with open(f"{control_save_path}", "w") as f1:
                        for key_ in low_default_threshold_sort[:11]:
                            # print(key_[0],low_default_threshold[key_[0]])
                            f1.write("\t".join(key_[0].split("@@@")) + "\n")
                else:
                    print("计算相关系数选取的对照数目小于10个！退出")
                    sys.exit()
            else:
                print("计算相关系数，移除输入文库号的数据库空!")
                sys.exit()


if __name__ == '__main__':
    fire.Fire()
