#-*-coding:utf-8-*-
"""外显子类"""

from intervaltools.core import IntervalBase, BucketIndexIntervalList
from exon_pipeline.utils.database_store import DatabaseStore

class IntervalDepth(IntervalBase):
    def __init__(self, chr, start, end, depth):
        super().__init__(chr, start, end)
        self.depth = self.format_factory(float)(depth)


class ExonDatabaseStore(DatabaseStore):
    def get_data(self, chr, start, end):
        """获取数据"""
        self.initialize()
        tag_interval = self.base(chr, start, end, )
        res_intervals = []
        for intervals in self._store:
            match_intervals = intervals.find_one(tag_interval, action='be_included')
            res_intervals.append(match_intervals)
        return res_intervals
    
    def get_depths(self, chr, start, end):
        res = super().get_data(chr, start, end)
        return [i.depth for i in res]


class IntervalContrast(IntervalBase):
    def __init__(self, chr, start, end, 
                 mean=0, cv=0, std=0, size=0, 
                 min_critical=0, max_critical=0, dps_mean=0):
        super().__init__(chr, start, end)
        self.mean = self.format_factory(float)(mean)
        self.cv = self.format_factory(float)(cv)
        self.std = self.format_factory(float)(std)
        self.size = self.format_factory(int)(size)
        self.min_critical = self.format_factory(float)(min_critical)
        self.max_critical = self.format_factory(float)(max_critical)
        self.dps_mean = self.format_factory(float)(dps_mean)