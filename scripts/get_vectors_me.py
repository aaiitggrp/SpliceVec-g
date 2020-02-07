import gensim
import numpy as np
import os
import sys

def generate_vectors(model_name,list_of_files,vec_choice):
	cur_dir=os.getcwd()
	len_of_files=list()
	files_generated_name=list()
	model=gensim.models.Word2Vec.load(cur_dir+'/data/'+model_name)
	vectors=model.syn0
	words=model.index2word
	#get vectors_helper returns two list;first contains the size of true and false data ;
	# 2nd one contains the name of the vector file generated 
	def get_vectors_helper(filename1,filename2):
		array_introns=np.empty((0,100)) #should be 100
		with open(cur_dir+"/corpus/"+filename1,"r") as infile:
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
		            #if(array_introns.shape[0]>=5612):
		                #break

		len1=array_introns.shape[0]
		with open(cur_dir+"/corpus/"+filename2,"r") as infile:
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
		                    mang +=vectors[ind1]
		                cha1 += 1
		            if (vec_choice==2):
		            	mang = mang / cha1
		            array_introns = np.vstack((array_introns, mang[None, :]))
		            print(array_introns.shape)
		            #if (array_introns.shape[0] >= 11224):
		             #   break
		len2=array_introns.shape[0]-len1
		print(array_introns.shape)
		if (count==1):
			output_filename=["all-chr_3mer_non-coding_withn(100_5_5_8_hs_iter20)_v26_" + "train " + "_40stream.nparray",]
		else:
			output_filename=["all-chr_3mer_non-coding_withn(100_5_5_8_hs_iter20)_v26_" + "test " + "_40stream.nparray",]
		with open(output_filename[0],"wb") as outfile:
		    np.save(outfile,array_introns)
		#output_filename=["all-chr_3mer_non-coding_withn(100_5_5_8_hs_iter20)_v26_test_40stream.nparray",]
		return [len1,len2,],output_filename








	count=1
	if(len(list_of_files)==2):
		len_of_files,files_generated_name=get_vectors_helper(list_of_files[0],list_of_files[1])
	else:
		temp1,temp2=get_vectors_helper(list_of_files[0],list_of_files[1])
		len_of_files.append(temp1[0],temp1[1])
		files_generated_name.append(temp2[0])
		del temp1[:];del temp2[:]
		count=2
		temp1,temp2=get_vectors_helper(list_of_files[2],list_of_files[3])
		len_of_files.append(temp1[0],temp1[1])
		files_generated_name.append(temp2[0])
		del temp1[:];del temp2[:]

	return len_of_files,list_of_files


