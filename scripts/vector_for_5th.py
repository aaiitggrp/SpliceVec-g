import numpy as np
import pandas as pd
import gensim
import tensorflow as tf
import os

def get_vectors_for_single_file(filename,vec_choice,model_name):
	model=gensim.models.Word2Vec.load(model_name)
	vectors=model.syn0
	words=model.index2word
	cur_dir=os.getcwd()
	array_introns=np.empty((0,100)) #should be 100
	with open(cur_dir+"/corpus/"+filename,"r") as infile:
	    for line in infile:
	        line = line.strip()
	        if line:
		        cha1 = 0
		        sw = 3
		        mang = np.zeros(100)
		        while cha1 <= len(line) - sw:
		            seq = line[cha1:cha1 + sw]
		            if seq in words:
		                ind1=words.index(seq)
		                mang = mang + vectors[ind1]
		            cha1 = cha1 + 1
		        print(array_introns.shape)
		        if(vec_choice==2):
		            mang=mang/cha1
		        array_introns=np.vstack((array_introns,mang[None,:]))

	return array_introns		        