import numpy as np
import os
## total_x made
total_x=np.load("all-chr_3mer_non-coding_withn(100_5_5_8_hs_iter15).intron_avg.nparray")

# total_y made
labels=np.array([0,1])
labels=np.tile(labels,(293889,1))  #<--------------  this number is not hardcoded replace
l21=np.array([1,0])
l21=np.tile(l21,(293889,1))  #<------------------ same here	
total_y=np.concatenate((labels,l21))


#####for i in range(len(total_x)):
#####	print(total_x[i],total_y[i])



# #no of elements dataset
# n_elem = len(total_x)
# print(n_elem)
# train_start = 0
# # 70% of actual data
# main_div = (n_elem*7)//10
# test_start = main_div
# test_end = n_elem-1
# # 80% to train
# train_end = (main_div*8)//10
# dev_start = train_end
# dev_end = main_div


# print(n_elem,dev_start,test_start)


# #permutate
# indices = np.random.permutation(n_elem)
# print(indices)
# np.savetxt("indices.txt",indices,fmt='%d')
# training_idx,dev_idx,test_idx = indices[0:dev_start], indices[dev_start:test_start], indices[test_start:]
# np.savetxt("training_idx.txt",training_idx,fmt='%d')
# np.savetxt("dev_idx.txt",dev_idx,fmt='%d')
# np.savetxt("test_idx.txt",test_idx,fmt='%d')
# #print(training_idx)
# #print(dev_idx)
# #print(test_idx)


#following implementation needed
#------------------------------------------------------------------------------------------------------------------------------------------------
#if the train test dev text files are not present, then 
	#if train test and dev idx.txt are not present, then create by above code. Else just load the idx to create train test dev text files.
#else no need of any partitioning. Data is already partitioned.

f3=open("training_idx.txt",'r')
f4=open("dev_idx.txt",'r')
f5=open("test_idx.txt",'r')

training_idx=np.loadtxt(f3,dtype="int64")
dev_idx=np.loadtxt(f4,dtype="int64")
test_idx=np.loadtxt(f5,dtype="int64")


train_x,dev_x,test_x = total_x[training_idx,:], total_x[dev_idx,:],total_x[test_idx,:]
train_y,dev_y,test_y = total_y[training_idx,:], total_y[dev_idx,:],total_y[test_idx,:]

np.savetxt("train_x.txt",train_x,fmt='%f')
np.savetxt("train_y.txt",train_y,fmt='%d')
np.savetxt("test_x.txt",test_x,fmt='%f')
np.savetxt("test_y.txt",test_y,fmt='%d')
np.savetxt("dev_x.txt",dev_x,fmt='%f')
np.savetxt("dev_y.txt",dev_y,fmt='%d')

count = 0
cou=0
for i in range(len(train_y)):
	if train_y[i][0]==1:
		count=count+1
	if train_y[i][0]==0:
		cou=cou+1
print('true', cou, 'false', count, 'out of',i+1, len(train_y),'training samples')

count = 0
cou=0
for i in range(len(dev_y)):
	if dev_y[i][0]==1:
		count=count+1
	if dev_y[i][0]==0:
		cou=cou+1
print('true', cou, 'false', count, 'out of',i+1,len(dev_y), 'dev samples')

count = 0
cou=0
for i in range(len(test_y)):
	if test_y[i][0]==1:
		count=count+1
	if test_y[i][0]==0:
		cou=cou+1
print('true', cou, 'false', count, 'out of',i+1,len(test_y),'test samples')
