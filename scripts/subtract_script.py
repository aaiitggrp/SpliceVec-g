import os
             #true introns                               #false introns
def get_actual_false_cords(unique_intron,negative_seq):
    cur_dir=os.getcwd()
    out_file_tobe_returned=negative_seq
    with open(cur_dir+'/corpus/'+unique_intron,'r') as one, open(cur_dir+'/corpus/'+negative_seq,'r') as two:
        set1=set()
        set2=set()
        for line in one:
            line=line.strip()
            line=line.split('\t')
            set1.add(tuple(line))
        for line in two:
            line=line.strip()
            line=line.split('\t')
            set2.add(tuple(line))
        print(len(set1))
        print(len(set2))
        set3=set2-set1
        print(len(set3))
        if len(set2) > len(set3):
            print("file created for true negative introns")
            three = open(cur_dir+'/corpus/'+'actual_negative_intron_26.tsv','w')
            out_file_tobe_returned='actual_negative_intron_26.tsv'
            for line in set3:
                three.write(line +'\n')
        else:
            print("no false negative cases found")
    	#rename negative_intron_26.tsv to actual_negative_intron_26.tsv
    return out_file_tobe_returned