import os
def splt_chr_code(temp,splt_str):
    cur_dir=os.getcwd()
    f1=open(cur_dir+'/corpus/'+temp,'r')
    new_list=list()
    chro=1
    for line in f1:
        if(line[0:2]==splt_str[0:2]):
            if(chro==23):
                fout = open(cur_dir+'/corpus/'+'chrX' + '.fa', 'w')
                new_list.append('chrX' + '.fa')
            elif(chro==24):
                fout = open(cur_dir+'/corpus/'+'chrY' + '.fa', 'w')
                new_list.append('chrY' + '.fa')
            else:
                fout = open(cur_dir+'/corpus/'+'chr' + str(chro) + '.fa', 'w')
                new_list.append('chr' + str(chro) + '.fa')
            chro+=1
           # if(chro<=24):
                #new_list.append('chr' + str(chro) + '.fa')
            continue
        elif(line[0:2]=='>G'):
            break
        else:
            line=line.strip()
            line=line.upper()
            fout.write(line)
    return new_list
