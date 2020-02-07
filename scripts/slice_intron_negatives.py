import numpy as np
import os

def get_false_intron_seq(input_filename):
    def complementer(string1):
        string2=str()
        complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A','N':'N'}
        for char in string1:
            string2+=complement[char]
        return string2
    output_filename='actual_negative_intron_26_subset.txt'
    cur_dir=os.getcwd()
    ar = np.genfromtxt(cur_dir+'/corpus/'+input_filename, dtype= 'U')
    num_rows, num_cols = ar.shape
    x = 0
    up_down_stream = 40
    fout = open(cur_dir+'/corpus/'+output_filename, "w")
    chr_name_change=1
    chr_name = ar[0][3]
    sense_chr=ar[0][2]
    while x < num_rows:
        if(chr_name != ar[x][3]):
            chr_name_change = 1
            chr_name = ar[x][3]
            sense_chr = ar[x][2]
        elif(sense_chr != ar[x][2]):
            sense_chr=ar[x][2]
        else:
            chr_name_change = 0
        if chr_name_change == 1 or x==0:
            filename =  chr_name
            fin = open(cur_dir+'/corpus/'+filename, "r")
        fin.seek(int(ar[x][0])-up_down_stream)
        string1= fin.read(int(ar[x][1])-int(ar[x][0])+2+(2*up_down_stream))
        if sense_chr=='+':
            fout.write(string1+'\n')
        else:
            string1=complementer(string1)
            fout.write(string1[::-1] + '\n')
        x = x + 1
    fin.close()
    fout.close()
    return output_filename
    	
    		

