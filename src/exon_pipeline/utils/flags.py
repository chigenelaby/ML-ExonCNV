#-*-coding:utf-8-*-
"""定义flag的元类"""

class FlagInfoBase:
    """
    标签信息类, 只包含了标签的信息
    flag_type 标签类型
    classify 具体标签类型
    show_string 展示字段
    comments 备注
    """
    def __init__(self, flag_type, classify, show_string, comments=''):
        self.flag_type = flag_type
        self.classify = classify
        self.show_string = show_string
        self.comments = comments
    
    def __eq__(self, other):
        return self.flag_type == other.flag_type and self.classify == other.classify
    
    def __str__(self):
        return "%s(type=%s, classify=%s)"%(self.__class__.__name__, self.flag_type, self.classify)
    
    def __repr__(self):
        return "%s(type=%s, classify=%s)"%(self.__class__.__name__, self.flag_type, self.classify)
    


class FlagScoreBase:
    """
    标签统计类定义
    flag_type标签类
    classify 具体标签类型
    score 计分
    show_string 展示字段
    flag 关联标签信息FlagInfoBase
    """
    def __init__(self, flag_type, classify, score, show_string, flag=None):
        self.flag_type = flag_type
        self.classify = classify
        self.score = score
        self.show_string = show_string
        self.flag = flag

    def __eq__(self, other):
        return self.flag_type == other.flag_type and self.classify == other.classify
    
    def __str__(self):
        return "%s(type=%s, classify=%s)"%(self.__class__.__name__, self.flag_type, self.classify)
    
    def __repr__(self):
        return "%s(type=%s, classify=%s)"%(self.__class__.__name__, self.flag_type, self.classify)

