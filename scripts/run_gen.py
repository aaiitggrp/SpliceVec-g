# import gensim
#
# sentences=list()
# with open('small.txt') as file1:
#     for line in file1:
#
#         sentences.append(line.split())
#
#
# model = gensim.models.Word2Vec(sentences, size=100, window=5, min_count=5, workers=4)
# model.save("super.model")

#
# import gensim
#
# model=gensim.models.Word2Vec.load("super.model")
# print(model.wv["AAA"])
# print(model.wv.most_similar(positive=["AAA"]))
import os
import gensim
import time
def generate_model():
    class MySentences(object):
        def __init__(self, dirname):
            self.dirname = dirname

        def __iter__(self):
            for fname in os.listdir(self.dirname):
                for line in open(os.path.join(self.dirname, fname)):
                    yield line.split()
    ww=0
    cur_dir=os.getcwd()
    for i in range(1): #can increase this for different window lengths
        ww+=5
        start=time.time()
        sentences = MySentences(cur_dir+'/corpus/split_chr_3mers')  # a memory-friendly iterator
        model = gensim.models.Word2Vec(sentences,size=100,window=ww,min_count=5,workers=8,negative=0,hs=1,iter=20)
        model.save(cur_dir+"/data/"+"all-chr_3mer_non-coding_withn(100_"+str(ww)+"_5_8_hs_iter20).model")
        end=time.time()
        print("elapsed: "+str((end-start)/60))
    model_name="all-chr_3mer_non-coding_withn(100_"+str(ww)+"_5_8_hs_iter20).model"
    return model_name
