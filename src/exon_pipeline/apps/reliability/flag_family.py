#-*-coding:utf-8-*-

"""获得家系flag
"""

from intervaltools.core import *
from mod_tools.mod_data_IO import *


FAMILY_INTERSECT_LIMIT = 0.5

class ConvFamily:
    def __init__(self, family_path=None):
        self.family_path = family_path
        self.family_data = self.load()
    
    def load(self):
        if self.family_path is None:
            family_data = BucketIndexIntervalList([])
        else:
            family_data = IO_pickle.load_pickle(self.family_path)
        return family_data
    
    def is_intersection(self, obj1, obj2):
        """相交比例>=0.5 FAMILY_INTERSECT_LIMIT"""
        return ((min(obj1.end, obj2.end) - max(obj1.start, obj2.start) + 1) / (obj1.end - obj1.start + 1)) >= FAMILY_INTERSECT_LIMIT
    
    def is_same_tag(self, conv1, conv2):
        """tag类型相同"""
        return conv1.tag[:4] == conv2.tag[:4]
    
    def get_flag(self, conv_obj):
        """
        flag True 家系相交且tag同
        flag None 无家系数据
        flag False 家系数据无支持数据
        """
        if self.family_path is None:
            return None
        intersect_exons = self.family_data.find_all(conv_obj, action='intersect')
        flag = True
        for ex in intersect_exons:
            if self.is_intersection(conv_obj, ex) and self.is_same_tag(conv_obj, ex):
                break
        else:
            flag = False
        return flag

class FlagCNVFamily:
    def __init__(self, family_path1=None, family_path2=None):
        self.family_obj1 = ConvFamily(family_path1)
        self.family_obj2 = ConvFamily(family_path2)
    
    def get_flag(self, conv_obj):
        """
        flag True 家系相交且tag同
        flag None 无家系数据
        flag False 家系数据无支持数据
        """
        flag1 = self.family_obj1.get_flag(conv_obj)
        flag2 = self.family_obj2.get_flag(conv_obj)
        if flag1 == flag2 == None:
            return None
        elif flag1 or flag2:
            return True
        else:
            return False