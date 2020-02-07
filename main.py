import os
import subprocess
import sys
cur_dir=os.getcwd()
sys.path.insert(0,cur_dir+'/scripts')
from split_chr import splt_chr_code
from break_scr import break_into_lines
from inp import break_into_3mers
from run_gen import generate_model
from scr1 import start_and_end_cord
from scr2 import remove_duplicates
from extract30bp import extract_threshold_bp
from slice_intron_positives import get_pos_seq
from false_data_coords import get_false_intron_cords
from subtract_script import get_actual_false_cords
from subset_negative_data import get_false_cords_subset
from slice_intron_negatives import get_false_intron_seq
from get_vectors_me import generate_vectors
from train_classifier_me import train_classifier_fun
from newop import test_on_new_data
#import all the scripts required
cur_dir=os.getcwd()
prev_model='all-chr_3mer_non-coding_withn(200_5_5_8_hs_iter20).model'
sys.path.insert(0,cur_dir+'/scripts')
def options():
	print("Choose from one of the following:")
	print("1 if u want to use given corpus for SpliceVec-g.")
	print("2 if u want to use your own corpus.")
	print("3 if you have true_data and false_data file.")
	print("4 if you have train_true_data, train_false_data, test_true_data and test_false_data")
	print("5 if you want to use the pretrained model for testing on new data set")

options()
choice=int(input())
print("Enter 1 for generating vectors using summation and 2 for using average:")
vec_choice=int(input())
if(choice==1):
	temp='GRCh38.p10.genome.fa'
	splt_str='>c'
elif(choice==2):
	print("Enter the filename:")
	temp=input()
	print("Enter the splitting character that separates one chromosome from another:")
	splt_str=input()
elif(choice==3 or choice==4):
	filename=[i for i in input("enter the filenames seperated by space:").split()]
#filename has names of file that have different chromosome sequence or name of true and false data file
if(choice==1 or choice ==2):
	seprate_chr_files=splt_chr_code(temp,splt_str)#call to splt_chr.py
	print("separate files for each chromosome is created")
	filename=break_into_lines(seprate_chr_files)#call to break.py change to break_scr.py
	print("Files created containing chromosomes broken into lines")
#now we will break the corpus files into 3mers

three_mers=break_into_3mers(filename) #call to inp.py
print("Files containing 3mers is created")
#model_name=generate_model()  #call to run_gen.py
model_name=prev_model  
print("model has been trained")
#now generate true and false data for case 1 and 2
if(choice==1 or choice==2):
	if(choice==1):
		gencode_filename='gencode.v26.annotation.gff3'
	elif(choice==2):
		gencode_filename=input("Enter the full name of gencode file:")
	new_dir=cur_dir+'/corpus'
	os.chdir(new_dir)
	command="cat "+gencode_filename+ " | awk '/\sexon\s|\stranscript\s|\sgene\s/ && /gene_type=protein_coding/' > gencode_protein_coding.gff"
	os.system(str(command))
	os.chdir(cur_dir)
	print("shell script done!!")
	intron_cord=start_and_end_cord('gencode_protein_coding.gff')  #call to scr1.py
	print("Start and end cordinates of introns generated!!")
	intron_cord=remove_duplicates(intron_cord) #call to scr2.py       
	print("Duplicate introns removed!!")
	intron_cord=extract_threshold_bp(intron_cord) #call to extract3bp.py      
	print("Entrons with length greater than 30bp is created!!")
	len_pos,positive_seq=get_pos_seq(intron_cord) #call to slice_intron_positives.py        
	print("Entronic sequences from the chromosomes is extracted")
	#false data generation
	negative_intron=get_false_intron_cords(seprate_chr_files) #call to false_data_coords.py 
	print("False intron cordinates is generated")
	negative_intron=get_actual_false_cords(intron_cord,negative_intron) #call to subtractscript.py 
	print("True cordinates are removed from false cords!!")
	actual_negative_seq=get_false_cords_subset(len_pos,negative_intron) #call to subset_negative_data.py
	print("number of negative cords reduced to number of positive cords")
	actual_negative_seq=get_false_intron_seq(actual_negative_seq) #call to slice_intron_negatives.py   <------------------
	print("false entron samples extracted from all chromosomes")
	#actual_negative_seq and positive_seq are true and false data for case 1 and 2
	files_for_vectors=[positive_seq,actual_negative_seq]


if(choice==3):
	files_for_vectors=filename

#call get_vectors_me.py for getting vectors
len_of_datasets,name_of_vector_files=generate_vectors(model_name,files_for_vectors,vec_choice)
print("Vectors are generated!!!")
saved_file_name=train_classifier_fun(len_of_datasets,name_of_vector_files) # call to train_classifier_me.py
print("classifier trained and saved")
if (choice==5):
	printf ("Enter the name of files in the format given in readme:")
	filename=[i for i in input.split()]
	test_on_new_data(filename,vec_choice,prev_model,'saved_weights')
print("-------------End of splice-vec-g programme!!!!-------------------")
