import pandas as pd
import joblib
import os
import numpy as np
import sklearn
# import fire
script_dir = os.path.dirname(os.path.abspath(__file__))


def singleexonfinetune(file1):
    # V5
    test_model = "model_part6"
    # test_model = "test2"
    # test_model = "test"
    # high_threshold = 0.3525
    # low_threshold = 0.25
    # V3
    # test_model = "model6"
    # high_threshold = 0.66
    # low_threshold = 0.25
    # V4
    # test_model = "model6"
    # high_threshold = 0.8
    # low_threshold = 0.49
    outfile = os.path.join(os.path.dirname(file1),"singleexonfinetune.txt")
    X_test = pd.read_csv(file1, sep='\t', keep_default_na=False)
    selected_columns = ['#chr', 'start', 'end', 'contrast.cv']
    selected_data = X_test[selected_columns].values

    # 1. 按diff值拆分X_test为三部分
    X_test_loss_ = X_test[X_test['diff'] < 0.75].copy()
    X_test_wild_ = X_test[(X_test['diff'] >= 0.75) & (X_test['diff'] < 1.25)].copy()
    X_test_gain_ = X_test[X_test['diff'] >= 1.25].copy()

    # 2. 加载三个不同模型并设置不同阈值
    model_loss_ = joblib.load(f'{script_dir}/best_pipeline_{test_model}_loss.pkl')
    model_wild_ = joblib.load(f'{script_dir}/best_pipeline_{test_model}_wild.pkl')
    model_gain_ = joblib.load(f'{script_dir}/best_pipeline_{test_model}_gain.pkl')

    # model_loss_ = joblib.load(f'{script_dir}/best_pipeline_{test_model}.pkl')
    # model_wild_ = joblib.load(f'{script_dir}/best_pipeline_{test_model}.pkl')
    # model_gain_ = joblib.load(f'{script_dir}/best_pipeline_{test_model}.pkl')

    # 设置不同阈值
    # model6_muti
    # thresholds = {
    #     'loss_': {'high': 0.3678, 'low': 0.25},  # diff < 0.75 的阈值
    #     'wild_': {'high': 0.0687, 'low': 0.05},  # 0.75 <= diff < 1.25 的阈值
    #     'gain_': {'high': 0.7633, 'low': 0.45}  # diff >= 1.25 的阈值
    # }
    # thresholds = {
    #     'loss_': {'high': 0.816, 'low': 0.25},  # diff < 0.75 的阈值
    #     'wild_': {'high': 0.816, 'low': 0.25},  # 0.75 <= diff < 1.25 的阈值
    #     'gain_': {'high': 0.9, 'low': 0.25}  # diff >= 1.25 的阈值
    # }
    # model6_muti_2

    # model_part5
    # thresholds = {
    #     'loss_': {'high': 0.3, 'low': 0.1558},  # diff < 0.75 的阈值
    #     'wild_': {'high': 0.2, 'low': 0.0944},  # 0.75 <= diff < 1.25 的阈值
    #     'gain_': {'high': 0.95, 'low': 0.9178}  # diff >= 1.25 的阈值
    # }
    # model_part6
    thresholds = {
        'loss_': {'high': 0.3, 'low': 0.1558},  # diff < 0.75 的阈值
        'wild_': {'high': 0.8, 'low': 0.4063},  # 0.75 <= diff < 1.25 的阈值
        'gain_': {'high': 0.95, 'low': 0.9037}  # diff >= 1.25 的阈值
    }



    # 3. 对各部分进行预测和分类
    def predict_and_classify(df, model, high_thresh, low_thresh):
        y_proba = model.predict_proba(df)[:, 1]
        reliability = np.where(
            y_proba >= high_thresh, "特别可靠",
            np.where(
                y_proba <= low_thresh, "不可靠",
                "相对可靠"
            )
        )
        return y_proba, reliability

    # 对三部分分别处理
    X_test_loss_['_y_proba'], reliability_loss_ = predict_and_classify(
        X_test_loss_, model_loss_, thresholds['loss_']['high'], thresholds['loss_']['low'])
    X_test_wild_['_y_proba'], reliability_wild_ = predict_and_classify(
        X_test_wild_, model_wild_, thresholds['wild_']['high'], thresholds['wild_']['low'])
    X_test_gain_['_y_proba'], reliability_gain_ = predict_and_classify(
        X_test_gain_, model_gain_, thresholds['gain_']['high'], thresholds['gain_']['low'])

    # 4. 合并三部分数据
    merged_df = pd.concat([
        X_test_loss_.assign(new_result=reliability_loss_),
        X_test_wild_.assign(new_result=reliability_wild_),
        X_test_gain_.assign(new_result=reliability_gain_)
    ]).sort_index()

    # 5. 按原逻辑处理列拆分和合并
    result_pos = merged_df.columns.get_loc('_result')
    before_result = merged_df.iloc[:, :result_pos]
    after_result = merged_df.iloc[:, result_pos+1:].drop(columns=['new_result', '_y_proba'], errors='ignore')
    original_result = merged_df['_result'].rename('原始可靠性标注')

    # 6. 创建最终DataFrame
    df_new = pd.concat([
        before_result,
        pd.Series(merged_df['new_result'], name='_result', index=merged_df.index),
        after_result,
        original_result,
        pd.Series(merged_df['_y_proba'], name='_y_proba', index=merged_df.index)
    ], axis=1)

    # 7. 输出结果
    df_new.to_csv(outfile, index=False, sep="\t")
    return selected_data.tolist()


##############################################将xgboost输出的初步结果进行校正
# 单外显子判定时：
# 1、loss: exon_num=1或总长<300， cv>0.15, |zscore|<2.5, _y_proba<0.5 特别可靠/相对可靠强制为不可靠， loss2除外
#    gain: exon_num=1或总长<300， cv>0.15或|zscore|<2.5, _y_proba<0.95 特别可靠/相对可靠强制为不可靠

# 2、loss： diff<0.55, cv<0.1 , |zscore|>4, y_proba>0.1  不可靠强制为相对可靠
# 3、gain： diff>1.45, cv<0.1 , |zscore|>5, y_proba>0.8  不可靠强制为相对可靠

# 4、exon_num==1 或总长<300 特别可靠强制为相对可靠

# thresholds = {
#         'loss_': {'high': 0.3, 'low': 0.1558},  # diff < 0.75 的阈值
#         'wild_': {'high': 0.8, 'low': 0.4063},  # 0.75 <= diff < 1.25 的阈值
#         'gain_': {'high': 0.95, 'low': 0.9037}  # diff >= 1.25 的阈值
#     }


# 20250529 对xgboost初步可靠性结果进行校正
def fmt_xgboost_relative(singleexonfinetune, output_file):
    with open(singleexonfinetune) as fo, open(output_file, "w") as fw:
        for line in fo:
            line = line.strip("\n").split("\t")
            if line[0] == "libId":
                print("\t".join(line), file=fw)
                keys = ["chrnameJoin", "_result", "diff", "z_score", "contrast.cv", "_y_proba", "_tag", "exon_num"]
                need_info_indx = [line.index(k) for k in keys]
                chrnameJoin_idx, _result_idx, diff_idx, z_score_idx, cv_idx, _y_proba_idx, tag_idx, exon_num_idx = need_info_indx
            else:
                need_info = [line[a] for a in need_info_indx]
                chrnameJoin, _result, diff, z_score, cv, _y_proba, tag, exon_num = need_info
                chrom, start, end = chrnameJoin.split(":")[0], chrnameJoin.split(":")[1].split("-")[0], chrnameJoin.split(":")[1].split("-")[1]
                # if (chrnameJoin == "chr1:891302-891393"):
                #     print(chrnameJoin, _result, diff, z_score, cv, _y_proba, tag, exon_num)
                # loss： diff<0.55, cv<0.1 , |zscore|>4, y_proba>0.1  不可靠强制为相对可靠
                if "loss" in tag and _result == "不可靠" and float(_y_proba) > 0.1 and abs(float(z_score)) > 4 and float(cv) < 0.1 and float(diff) < 0.55:
                    line[_result_idx] = "相对可靠"
                
                # gain： diff>1.45, cv<0.1 , |zscore|>5, y_proba>0.8  不可靠强制为相对可靠
                if "gain" in tag and _result == "不可靠" and float(_y_proba) > 0.8 and abs(float(z_score)) > 5 and float(cv) < 0.1 and float(diff) > 1.45:
                    line[_result_idx] = "相对可靠"

                # loss: exon_num=1或总长<300， cv>0.15, |zscore|<2.5, _y_proba<0.5 特别可靠/相对可靠强制为不可靠， loss2除外
                # gain: exon_num=1或总长<300， cv>0.15或|zscore|<2.5, _y_proba<0.95 特别可靠/相对可靠强制为不可靠
                if line[1]=='chr1:120057260-120057460': print(exon_num, _result, tag, _y_proba, z_score, cv)
                if (int(exon_num) == 1 or int(end) - int(start) < 300) and _result != "不可靠":
                    if "loss1" in tag and float(_y_proba) < 0.5 and abs(float(z_score)) < 2.5 and float(cv) > 0.15:
                        line[_result_idx] = "不可靠"
                    # if "gain" in tag and float(_y_proba) < 0.95 and (abs(float(z_score)) < 2.5 or float(cv) > 0.15):
                    if "gain" in tag and (abs(float(z_score)) < 2.5 or float(cv) > 0.15):
                        line[_result_idx] = "不可靠"
                
                # exon_num==1 或总长<300 特别可靠强制为相对可靠
                if (int(exon_num) == 1 or int(end) - int(start) < 300) and line[_result_idx] == "特别可靠" and tag != "loss2":
                    line[_result_idx] = "相对可靠"
                print("\t".join(line), file=fw)


# if __name__ == "__main__":
#     fire.Fire(fmt_xgboost_relative)
