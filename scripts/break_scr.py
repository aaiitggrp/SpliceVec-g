import os
def break_into_lines(seprate_chr_files):
    #iterlist=[i+1 for i in range(22)]
    #iterlist.extend(['X','Y'])
    files_in_lines=list()
    cur_dir=os.getcwd()
    for i in seprate_chr_files:
        #filename = 'chr' + str(i) + '.fa'
        outfile=str(i)+'.split'
        files_in_lines.append(outfile)
        with open(cur_dir+'/corpus/'+str(i),'r') as inf, open(cur_dir+'/corpus/'+outfile,'w') as outf:
            string=str()
            while(True):
                buf=inf.read(1800)
                if not buf:
                    break
                buf=buf.upper()
                outf.write(buf+'\n')
    return files_in_lines
