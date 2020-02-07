import gensim
import numpy as np
import os
model = gensim.models.Word2Vec.load('all-chr_3mer_non-coding_withn(250_5_5_8_hs_iter20).model') #<------ should be 200_5
vectors=model.syn0
words=model.index2word
array_introns=np.empty((0,250)) #should be 100
with open("/home/d.aparajita/word2vec_spliceVec/test_true.txt","r") as infile:
    for line in infile:
        line = line.strip()
        if line:
            cha1 = 0
            sw = 3
            mang = np.zeros(250)
            while cha1 <= len(line) - sw:
                seq = line[cha1:cha1 + sw]
                if seq in words:
                    ind1=words.index(seq)
                    mang = mang + vectors[ind1]
                cha1 = cha1 + 1
            print(array_introns.shape)
            mang=mang/cha1
            array_introns=np.vstack((array_introns,mang[None,:]))
            if(array_introns.shape[0]>=5612):
                break

with open("/home/d.aparajita/word2vec_spliceVec/test_false.txt","r") as infile:
    for line in infile:
        line = line.strip()
        if line:
            cha1 = 0
            sw = 3
            mang = np.zeros(250)
            while cha1 <= len(line) - sw:
                seq = line[cha1:cha1 + sw]
                if seq in words:
                    ind1=words.index(seq)
                    mang +=vectors[ind1]
                cha1 += 1
            mang = mang / cha1
            array_introns = np.vstack((array_introns, mang[None, :]))
            print(array_introns.shape)
            if (array_introns.shape[0] >= 11224):
                break
print(array_introns.shape)
with open("all-chr_3mer_non-coding_withn(250_5_5_8_hs_iter20)_v26_test_40stream.nparray","wb") as outfile:
    np.save(outfile,array_introns)
