import random
import os
def get_false_intron_cords(filename_list):
    cur_dir=os.getcwd()
    def count_char(filename):
        f1=open(cur_dir+"/corpus/"+filename,'r')
        char=0
        for line in f1:
            char+=len(line)
        return char
    cur_dir=os.getcwd()
    intron_min=30
    intron_max=1240200
    #iterlist=[i+1 for i in range(22)]
    #iterlist.extend(['X','Y'])
    dictoflen={}
    for i in filename_list:
        #filename='chr'+str(i)+'.fa'
        filename=i
        dictoflen[filename]=count_char(filename)

    for i in filename_list:
        #filename = 'chr' + str(i) + '.fa'
        filename=i
        with open(cur_dir+'/corpus/'+filename,'r') as infi,open(cur_dir+'/corpus/'+'negative_intron_26.tsv','a') as outfile:
            setofneg = set()
            random.seed(int(dictoflen[filename]))
            while(len(setofneg)<40000):
                seekpos=random.randint(0,int(dictoflen[filename])-1-intron_max-10000)
                infi.seek(seekpos)
                string1=infi.read(10000)
                start=string1.find('GT')
                if(start != -1):
                    startpos=seekpos+start
                    infi.seek(startpos+intron_min)
                    string2=infi.read(intron_max-intron_min)
                    seekendpos=random.randint(0,intron_max-intron_min)
                    stoppos=string2[seekendpos:].find('AG')
                    if(stoppos != -1):
                        stoppos+=startpos+intron_min+seekendpos
                        setofneg.add((startpos,stoppos))
                else:
                    continue
            for ele in setofneg:
                outfile.write(str(ele[0])+'\t'+str(ele[1])+'\t'+'+'+'\t'+str(i)+'\n')
        return ('negative_intron_26.tsv')




