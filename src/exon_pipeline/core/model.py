#-*-coding:utf-8-*-
"""对照库离群算法模型"""

import numpy as np

class ContrastModelBase:
    """对照库离群算法模型
    输入应该是np.array 数组
    """
    name = ''
    description = ''
    def get_normal(self, arr):
        """获取正常库"""
        pass


class ContrastModel1(ContrastModelBase):
    """对照库离群算法模型1"""
    name = '百分位离群检验法'
    description = """IQR = Q3 - Q1; (arr >= (Q1 - 1.5 * IQR)) & (arr <= (Q3 + 1.5 * IQR))"""
    def get_normal(self, arr):
        """获取正常库"""
        try:
            Q1 = np.percentile(arr, 25)
            Q3 = np.percentile(arr, 75)
            IQR = Q3 - Q1
            res = arr[(arr >= (Q1 - 1.5 * IQR)) & (arr <= (Q3 + 1.5 * IQR))]
        except:
            res = arr
        return res

class ContrastModel2(ContrastModelBase):
    """对照库离群算法模型1"""
    name = '拉依达法(k=3)'
    description = """|x_out-x_mean|>k*std"""
    def get_normal(self, arr):
        """获取正常库"""
        if not len(arr):
            return arr
        try:
            mean = arr.mean()
            std = arr.std(ddof=1)
            res = arr[abs(arr - mean) / std < 3]
        except:
            res = arr
        return res


class ContrastModel3(ContrastModelBase):
    """对照库离群算法模型1"""
    name = '拉依达法(k=2.5)'
    description = """|x_out-x_mean|>k*std"""
    def get_normal(self, arr):
        """获取正常库"""
        if not len(arr):
            return arr
        try:
            mean = arr.mean()
            std = arr.std(ddof=1)
            res = arr[abs(arr - mean) / std < 2.5]
        except:
            res = arr
        return res

models_dict = {'model1': ContrastModel1(), 
               'model2': ContrastModel2(), 
               'model3': ContrastModel3()}