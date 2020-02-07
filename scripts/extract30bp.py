import os
def extract_threshold_bp(infilename):
    cur_dir=os.getcwd()
    with open(cur_dir+'/corpus/'+infilename,'r') as f1, open(cur_dir+'/corpus/'+'unique_26_intron_30bp.tsv','w') as outfile:
        c3=0
        for line in f1:
            line1=line.split('\t')
            c1=int(line1[0])
            c2=int(line1[1])
            if(c2-c1>30):
                outfile.write(line)
            if(c3<c2-c1):
                c3=c2-c1

        print(str(c3))
    return ('unique_26_intron_30bp.tsv')

