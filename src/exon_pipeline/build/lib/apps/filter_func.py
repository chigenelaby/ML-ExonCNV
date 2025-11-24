#-*-coding:utf-8-*-
"""提供过滤的标准"""
from operator import itemgetter
from intervaltools.core import IntervalBase, BucketIndexIntervalList
from exon_pipeline.utils.database_store import get_index

def filter_NA_key(conv):
    """过滤NA"""
    if conv.tag == 'NA':
        return False
    return True

def filter_key(conv):
    """
    过滤key
    True 不过滤
    False 过滤
    """
    # if conv.chr == "chrX":
    #     print(conv.chr, conv.start, conv.end, conv.tag, conv.flag_diff.classify, conv.flag_is_cnv.classify)

    if conv.tag == 'NA':
        return False
    if conv.flag_diff.classify == False:
        return False
    if conv.tag == 'loss2':
        return True
    if conv.flag_is_cnv.classify:
        return True
    if hasattr(conv, "score") and conv.score > 0:
        return True
    #
    # flag_contrast = (len([ex for ex in conv.exons if ex.contrast.cv<0.2])/len(conv.exons)>0.5)
    # 如果区段恰好包含两个exon，则原阈值需要两个exon均满足条件
    flag_contrast = (len([ex for ex in conv.exons if ex.contrast.cv<0.2])/len(conv.exons)>=0.5)

    # if conv.chr == "chrX":
    #     print(flag_contrast, [ex.contrast.cv for ex in conv.exons], [ex.flag_outlier for ex in conv.exons], len(conv.exons))
    
    if len(conv.exons) < 3:
        if conv.freq < 0.05 and flag_contrast:
            print(conv, "1111")
            # return True
    if len(conv.exons) > 2 and flag_contrast:
        print(conv, "2222")
        # return True
    # return False
    return True



class FilterConditionBase:
    """过滤条件基类"""
    def check(self, conv):
        return True

class FilterConditionList(FilterConditionBase):
    """过滤条件基类"""
    def __init__(self, list_path):
        self.list_convs = get_index(list_path, key=itemgetter(0, 1, 2))
    
    def check(self, conv):
        if self.list_convs.find_one(conv, action='intersect'):
            return True
        return False

class FilterConditionAttentions(FilterConditionList):
    """gene白名单
    __init__ list_path
    """
    def check(self, conv):
        if self.list_convs.find_one(conv, action='intersect'):
            return True
        return False

## by housy
class FilterConditionAttentions_condition(FilterConditionList):
    """gene白名单2档：保留基因白名单，并且长度>=3个外显子的conv
    __init__ list_path
    """
    def check(self, conv):
        if self.list_convs.find_one(conv, action='intersect') and len(conv.exons) > 2 :
            return True
        return False

class FilterConditionListFilterFalse(FilterConditionList):
    """gene黑名单"""
    def intersect_length(self, obj1, obj2):
        """相交长度"""
        return min(obj1.end, obj2.end) - max(obj1.start, obj2.start) + 1
    
    def check(self, conv):
        """相交长度比例占目标区间0.5以上"""
        cvs = self.list_convs.find_all(conv, action='intersect')
        if not cvs:
            return True
        intersect_length = sum(map(lambda x:self.intersect_length(conv, x), cvs))
        return intersect_length / len(conv) <= 0.5
    
    def check_gene(self, conv):
        if self.list_convs.find_one(conv, action='intersect'):
            return False
        return True


class ConditionSet:
    def __init__(self, base=True):
        self._conditions = []
        self._models = []
        self._filter_models = []
        self.base = base
    
    def _add_condition(self, condition, model='instance'):
        assert model in ['instance', 'function']
        self._conditions.append(condition)
        self._models.append(model)
    
    def add_filter_condition(self, condition, model='instance'):
        """过滤标准，True为不过滤"""
        self._add_condition(condition, model)
        self._filter_models.append(True)
    
    def add_attention_condition(self, condition, model='instance'):
        """保留数据"""
        self._add_condition(condition, model)
        self._filter_models.append(False)
    
    def check(self, *argv, **keywords):
        for condition, modle, filter_model in zip(self._conditions, self._models, self._filter_models):
            if modle == 'instance':
                condition = condition.check
            v = condition(*argv, **keywords)
            if filter_model:
                if not v:
                    return False
            if not filter_model:
                if v:
                    return True
        return self.base


