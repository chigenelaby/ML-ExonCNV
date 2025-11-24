#-*-coding:utf-8-*-
"""
ExonConv.utils.load_qc 模块
"""
from mod_tools.mod_data_IO import *

class QCFormatError(Exception):
    pass

class LoadQCreport:
    def __init__(self, qc_report_path):
        self.qc_report_path = qc_report_path
        self.qc_info_dict = self.get_info_map()
    
    def get_info_map(self):
        try:
            infos = load_file(self.qc_report_path, '\t', '#')
            info_map = dict(filter(lambda x:len(x)==2, infos))
        except:
            raise QCFormatError
        return info_map
    
    def _check_easy(self, value):
        """简易检查性别字段"""
        sexuality_flag_key = 'YES'
        return sexuality_flag_key in value
    
    def get_sexuality(self):
        """rtype: 'XX' or 'XY'"""
        sexuality_key1 = 'sex_test:XY'
        sexuality_key2 = 'sex_test:XX'
        
        for sexuality_key in [sexuality_key1, sexuality_key2]:
            try:
                value = self.qc_info_dict[sexuality_key]
                break
            except KeyError:
                continue
        else:
            raise QCFormatError
        
        if not self._check_easy(value):
            raise QCFormatError(sexuality_key)
        
        sexuality = sexuality_key.split(':')[-1]
        
        return sexuality

def get_sexuality_from_qc_report(qc_report_path):
    """获取性别字段从qc文件。如果文件错误引发QCFormatError"""
    lq = LoadQCreport(qc_report_path)
    return lq.get_sexuality()


