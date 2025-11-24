# 版本参数

import re
from operator import itemgetter
from exon_pipeline.utils.sample_tools import get_info, get_info_reverse
from exon_pipeline.utils.load_qc import get_sexuality_from_qc_report
from mod_tools.mod_data_IO import *

def _get_sexuality(sex_path):
    data = load_file(sex_path, '\t', '#')
    for i in data:
        if i[0] == 'Final_sex:':
            assert i[1] in {'XY', 'XX'}
            return i[1]
    raise Exception('性别异常')


# exon_pipeline.apps.vaf
vcf_itemgetter = itemgetter(0,1,2,5)
vcf_regex = re.compile('.+/.+=(.+)()')

# exon_pipeline.core.contrast
fun_contrast = get_info_reverse
STAT_DEPTH_STR = '平均深度'


get_sexuality = _get_sexuality


# 不均衡系数
UNBALANCE_LIMIT1 = 0.02
UNBALANCE_LIMIT2 = 0.08