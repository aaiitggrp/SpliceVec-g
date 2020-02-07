import os
import time
 # actual  code which i used to prepare input
# for chrc in range(22):
#     chrn='../all_chr/'+'chr'+str(chrc+1)+'.fa'
#     with open(chrn) as file1:
#         count=0
#         iter=True
#         outfile='chr'+str(chrc+1)+'.fa'+'.split'
#         out=open(outfile,'w')
#         for line in file1:
#             if iter:
#                 iter=False
#                 continue
#             line=line.strip()
#             line=line.replace('N','')
#             line = line.replace('n', '')
#             if(len(line)>0):
#                 out.write(line.upper())
#                 count+=1
#                 if count>10:
#                     count=0
#                     out.write('\n')
# end1=time.time()
# print('phase 1 complete...')
# print("elapsed time= "+str((end1-start)/60))



#this was used to make sliding window also for 6 and 9
def break_into_3mers(filename):
    #iterlist=[i+1 for i in range(22)]
    #iterlist.extend(['X','Y'])
    output_list=list()
    cur_dir=os.getcwd()
    sw=0
    if not os.path.isdir(cur_dir+'/corpus/split_chr_3mers'):
        os.makedirs(cur_dir+'/corpus/split_chr_3mers')
    if not os.path.isdir(cur_dir+'/data'):
        os.makedirs(cur_dir+'/data')
    for i in range(1): #range can be increased to generate 6,9 etc k-mers
        start = time.time()
        sw+=3
        for chrc in filename:
            #chrn = 'chr' + str(chrc) + '.fa'
            infile=chrc
            with open(cur_dir+'/corpus/'+infile) as file1:
                outf=infile+'.'+str(sw)+'mer'
                output_list.append(outf)
                outfile=open(cur_dir+'/corpus/split_chr_3mers/'+outf,'w')
                count=0
                for line in file1:
                    line=line.strip()
                    if line:
                        cha1=0
                        mang=""
                        while cha1<=len(line)-sw:
                            mang+=line[cha1:cha1+sw]
                            mang+=' '
                            cha1+=1
                        mang.strip()
                        outfile.write(mang+'\n')
        end2=time.time()
        #print('phase 2 complete...')
        #print("elapsed time= "+str((end2-end1)/60))
        print("total time elapsed= "+str((end2-start)/60))
    return output_list
