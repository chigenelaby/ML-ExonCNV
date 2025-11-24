#-*-coding:utf-8-*-

"""
家系类
"""

from intervaltools.core import *
from mod_tools.mod_data_IO import *


FAMILY_INTERSECT_LIMIT = 0.5

class ConvFamilyBase:
    def is_intersection(self, obj1, obj2):
        """相交比例>=0.5 FAMILY_INTERSECT_LIMIT"""
        return ((min(obj1.end, obj2.end) - max(obj1.start, obj2.start) + 1) / (obj1.end - obj1.start + 1)) >= FAMILY_INTERSECT_LIMIT
    
    def is_same_tag(self, tag1, tag2):
        """tag类型相同"""
        return tag1[:4] == tag2[:4]
    
    def get_intersect_conv(self, conv_obj, tag , convs):
        """
        获取在家系convs中与conv_obj能校正的区间
        """
        intersect_convs = convs.find_all(conv_obj, action='intersect')
        for conv in intersect_convs:
            if self.is_intersection(conv_obj, conv) and self.is_same_tag(tag, conv.tag):
                return conv
        return None



