import fire
from mod_tools.mod_data_IO import *

def mod_db(freq_pic, new_freq_pic):
	a = IO_pickle.load_pickle(freq_pic)
	for att in ['exons', 'flag_cnv_family', 'flag_cnvseq', 'flag_is_cnv', 'flag_length', 'flag_reliable_exon', 'flag_vaf', 'is_male_sex_chr', 'score', 'diff', 'flag_diff']:
		for i in a:
			for j in i:
				if hasattr(j, att):
					delattr(j, att)
	IO_pickle.write_pickle(a, new_freq_pic)


if __name__ == "__main__":
	fire.Fire(mod_db)

