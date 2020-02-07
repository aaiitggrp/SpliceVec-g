from more_itertools import unique_everseen
import os
def remove_duplicates(input_file):
	cur_dir=os.getcwd()
	with open(cur_dir+'/corpus/'+input_file,'r') as f, open(cur_dir+'/corpus/'+'unique_intron26.tsv','w') as out_file:
	   out_file.writelines(unique_everseen(f))
	return ('unique_intron26.tsv')
