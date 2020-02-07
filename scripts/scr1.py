import os
def start_and_end_cord(filename):
    cur_dir=os.getcwd()
    name_of_output_file='intron26.tsv'
    with open(cur_dir+'/corpus/'+filename) as infi, open(cur_dir+'/corpus/'+'intron26.tsv','w') as outfile:
        list_one=list()
        for line in infi:
            tet=line.split("\t")
            if(tet[2]=='transcript'):
                if(list_one):
                    l2=sorted(list_one)
                    for w in range(len(l2)-1):
                        outfile.write(str(l2[w][1]+1)+'\t'+str(l2[w+1][0]-1)+'\t'+str(l2[w][2])+'\t'+str(l2[w][3])+'\n')
                    list_one=list()
            if(tet[2]=='exon'):
                list_one.append((int(tet[3]),int(tet[4]),str(tet[6]),str(tet[0])))
        if (list_one):
            l2 = sorted(list_one)
            for w in range(len(l2) - 1):
                outfile.write(
                    str(l2[w][1] + 1) + '\t' + str(l2[w + 1][0] - 1) + '\t' + str(l2[w][2]) + '\t' + str(l2[w][3]) + '\n')
            list_one = list()
    return name_of_output_file
