#-*-coding:utf-8-*-

import os
import logging

logger = logging.getLogger(__name__)

class Handler:
    def _successor(self, successor):
        self.successor = successor
    
    def handle(self, *argv):
        return self.successor.handle(*argv)

class DeltaDiff:
    """(delta_diff)diff值差的相关类。
    假设：缺失重复的准确性与diff值差有关，并受到对照库深度影响
    注：所有区间范围都是上闭下包
    默认深度为100
    """
    _loss1_base = 0.5
    _loss2_base = 0
    _gain1_base = 1.5
    _gain2_base = 2
    lower_limit = 0
    upper_limit = 2
    _epsilon = 10 ** (-9)
    default_depth = 100
    
    def __init__(self, is_initialized=True):
        self.tags = []
        self.tag_condition = {}
        self._tag_obj = {}
        if is_initialized:
            self.initialize()
    
    def initialize(self):
        self.set_tag('loss1', self._loss1_base)
        self.set_tag('loss2', self._loss2_base)
        self.set_tag('gain1', self._gain1_base)
        self.set_tag('gain2', self._gain2_base)

    def limit_fomat(self, value):
        if value > self.upper_limit:
            value = self.upper_limit
        elif value < self.lower_limit:
            raise ValueError('diff value need be greater than %s'%self.lower_limit)
        return value
    

    def set_tag(self, tag, base_vale=None):
        setattr(self, tag, tag)
        setattr(self, '%s_base'%tag, base_vale)
        self.tags.append(tag)
    
    def _get_interval_func(self, interval=(None, None), limit=(None, None), return_value=1):
        """区间值判定函数"""
        lower, upper = interval
        lower_limit, upper_limit = limit
        if lower is None:
            lower = lower_limit - self._epsilon
        if upper is None:
            upper = upper_limit
        def _func(diff_value):
            if lower < diff_value and diff_value <= upper:
                return return_value
            else:
                return 0
        return _func
    
    def set_depth_ddiff(self, tag, depth_interval=(None, None), diff_interval=(None, None), return_value=1, is_uniqueness=True):
        """设置变量范围条件
        tag：缺失重复标签
        depth_interval 深度范围
        diff_interval diff值范围
        return_value 满足范围时的取值
        is_uniqueness 默认为True, 即depth_interval为唯一范围。
                      如果为False，则允许同一深度下有不同的标准返回值。这时深度范围应当不能省略。
        可多次添加，每次添加按elif逻辑判定。
        """
        assert tag in self.tags
        depth_limit = (0, 10 ** 6)
        diff_limit = (self.lower_limit, self.upper_limit)
        depth_func = self._get_interval_func(depth_interval, depth_limit)
        diff_func = self._get_interval_func(diff_interval, diff_limit)
        class ConcreteHandler(Handler):
            def handle(self, depth, diff_value):
                if depth_func(depth):
                    if diff_func(diff_value):
                        return return_value
                    elif is_uniqueness:
                        return 0
                    else:
                        try:
                            return self.successor.handle(depth, diff_value)
                        except:
                            return 0
                else:
                    try:
                        return self.successor.handle(depth, diff_value)
                    except:
                        return 0
        
        obj = self.tag_condition.setdefault(tag, Handler())
        while True:
            try:
                suc = obj.successor
                obj = suc
            except:
                obj._successor(ConcreteHandler())
                break

    def forward(self, tag, diff_value, depth=None):
        """返回diff设置的值"""
        if depth is None:
            depth = self.default_depth
        diff_value = self.limit_fomat(diff_value)
        condition = self.tag_condition.get(tag)
        return condition.handle(depth, diff_value)

diff_condition = DeltaDiff()
# diff 不可靠条件
# diff_condition.set_depth_ddiff('loss2', diff_interval=(0.15, None), return_value=-1, is_uniqueness=False)
diff_condition.set_depth_ddiff('loss2', diff_interval=(0.2, None), return_value=-1, is_uniqueness=False)
# diff_condition.set_depth_ddiff('loss1', diff_interval=(None, 0.15), return_value=-2, is_uniqueness=False)
diff_condition.set_depth_ddiff('loss1', diff_interval=(None, 0.1), return_value=-2, is_uniqueness=False)
diff_condition.set_depth_ddiff('loss1', diff_interval=(0.7, 1.3), return_value=-1, is_uniqueness=False)
diff_condition.set_depth_ddiff('loss1', diff_interval=(1.3, None), return_value=-2, is_uniqueness=False)
diff_condition.set_depth_ddiff('gain1', diff_interval=(0.7, 1.3), return_value=-1, is_uniqueness=False)
diff_condition.set_depth_ddiff('gain1', diff_interval=(None, 0.7), return_value=-2, is_uniqueness=False)
# diff_condition.set_depth_ddiff('gain1', diff_interval=(1.7, None), return_value=-1, is_uniqueness=False)
diff_condition.set_depth_ddiff('gain1', diff_interval=(1.8, None), return_value=-1, is_uniqueness=False)
diff_condition.set_depth_ddiff('gain2', diff_interval=(0.7, 1.8), return_value=-1, is_uniqueness=False)
diff_condition.set_depth_ddiff('gain2', diff_interval=(None, 0.7), return_value=-2, is_uniqueness=False)
# diff 可靠条件
diff_condition.set_depth_ddiff('loss2', diff_interval=(None, 0.03), return_value=2, is_uniqueness=False)
# diff_condition.set_depth_ddiff('loss1', diff_interval=(0.40, 0.55), return_value=1, is_uniqueness=False)
diff_condition.set_depth_ddiff('loss2', diff_interval=(0.03, 0.1), return_value=1, is_uniqueness=False)
diff_condition.set_depth_ddiff('loss1', diff_interval=(0.1, 0.55), return_value=1, is_uniqueness=False)

# diff_condition.set_depth_ddiff('gain1', diff_interval=(1.45, 1.6), return_value=1, is_uniqueness=False)
diff_condition.set_depth_ddiff('gain1', diff_interval=(1.45, 1.8), return_value=1, is_uniqueness=False)
# diff_condition.set_depth_ddiff('gain2', diff_interval=(1.85, None), return_value=1, is_uniqueness=False)
diff_condition.set_depth_ddiff('gain2', diff_interval=(1.8, None), return_value=1, is_uniqueness=False)
