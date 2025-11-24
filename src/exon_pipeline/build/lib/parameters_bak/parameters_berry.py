# 版本参数

import re
from operator import itemgetter
from exon_pipeline.utils.sample_tools import get_info, get_info_reverse
from exon_pipeline.utils.load_qc import get_sexuality_from_qc_report


# exon_pipeline.apps.vaf
vcf_itemgetter = itemgetter(0,1,2,8)
vcf_regex = re.compile('.+/.+=(.+?);(.+?);')

# exon_pipeline.core.contrast
fun_contrast = get_info
STAT_DEPTH_STR = 'Average_rmdupdepth_point'


get_sexuality = get_sexuality_from_qc_report

# 不均衡系数
UNBALANCE_LIMIT1 = 0.03
UNBALANCE_LIMIT2 = 0.06