import random
import os
def get_false_cords_subset(len_pos_seq,input_filename):
	cur_dir=os.getcwd()
	sampler=random.sample(range(1, 1920000), len_pos_seq)
	sampler=set(sampler)

	with open(cur_dir+'/corpus/'+input_filename,'r') as two,open(cur_dir+'/corpus/'+'actual_negative_intron_26_subset.tsv','w') as fillu:
	    count=0
	    for line in two:
	        count+=1
	        if(count in sampler):
	            line=line.strip()
	            fillu.write(line+'\n')
	return ('actual_negative_intron_26_subset.tsv')
