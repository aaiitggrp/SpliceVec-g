import tensorflow as tf
import pickle 
import time
import  numpy as np
import sklearn as sk
from sklearn.model_selection import train_test_split
import os
import sys
def train_classifier_fun(len_of_files,name_of_nparray_files):
	saver=tf.train.Saver()
	cur_dir=os.getcwd()
	if (len(len_of_files)==4):
		train_x=np.load(cur_dir+"/corpus/"+name_of_nparray_files[0])
		test_x=np.load(cur_dir+"/corpus/"+name_of_nparray_files[1])
		labels=np.array([0,1])
		labels=np.tile(labels,(len_of_files[0],1))
		l21=np.array([0,1])
		l21=np.tile(l21,(len_of_files[1],1))
		train_y=np.concatenate((labels,l21))
		labels=np.array([0,1])
		labels=np.tile(labels,(len_of_files[2],1))
		l21=np.array([1,0])
		l21=np.tile(l21,(len_of_files[3],1))
		test_y=np.concatenate((labels,l21))
		train_x,_,train_y,_=train_test_split(train_x,train_y,test_size=0,random_state=0)


	else:
		total_x=np.load(cur_dir+"/corpus/"+name_of_nparray_files[0])
		labels=np.array([1,0])
		labels=np.tile(labels,(len_of_files[0],1))
		l21=np.array([0,1])
		l21=np.tile(l21,(len_of_files[1],1))
		total_y=np.concatenate((labels,l21))
		n_elem=len_of_files[0]+len_of_files[1]
		indices=np.random.permutation(n_elem)
		main_div=(n_elem*7)//10
		test_start=main_div
		test_end=n_elem-1
		#80% to train
		train_end=(main_div*8)//10
		dev_start=train_end
		dev_end=main_div
		training_idx,dev_idx,test_idx = indices[0:dev_start], indices[dev_start:test_start], indices[test_start:]
		train_x,dev_x,test_x = total_x[training_idx,:], total_x[dev_idx,:],total_x[test_idx,:]
		train_y,dev_y,test_y = total_y[training_idx,:], total_y[dev_idx,:],total_y[test_idx,:]


	# train_x and test_x has been prepared...
	#now prepare nn and train model

	n_nodes_hl1 = 2500
#n_nodes_hl2 = 2500
#n_nodes_hl3 = 2500

	n_classes = 2
	batch_size = 128
	hm_epochs = 50

	x = tf.placeholder('float',name='X')
	y = tf.placeholder('float',name='Y')

	hidden_1_layer = {'f_fum':n_nodes_hl1,
	                  'weight':tf.Variable(tf.random_normal([len(train_x[0]), n_nodes_hl1])),
	                  'bias':tf.Variable(tf.random_normal([n_nodes_hl1]))}

	# hidden_2_layer = {'f_fum':n_nodes_hl2,
	#                    'weight':tf.Variable(tf.random_normal([n_nodes_hl1, n_nodes_hl2])),
	#                    'bias':tf.Variable(tf.random_normal([n_nodes_hl2]))}

	# hidden_3_layer = {'f_fum':n_nodes_hl3,
	#                  'weight':tf.Variable(tf.random_normal([n_nodes_hl2, n_nodes_hl3])),
	#                   'bias':tf.Variable(tf.random_normal([n_nodes_hl3]))}

	output_layer = {'f_fum':None,
	                'weight':tf.Variable(tf.random_normal([n_nodes_hl1, n_classes])),
	                'bias':tf.Variable(tf.random_normal([n_classes])),}


	# Nothing changes
	def neural_network_model(data):

	    l1 = tf.add(tf.matmul(data,hidden_1_layer['weight']), hidden_1_layer['bias'])
	    l1 = tf.nn.relu(l1)

	    # l2 = tf.add(tf.matmul(l1,hidden_2_layer['weight']), hidden_2_layer['bias'])
	    # l2 = tf.nn.relu(l2)

	    # l3 = tf.add(tf.matmul(l2,hidden_3_layer['weight']), hidden_3_layer['bias'])
	    # l3 = tf.nn.relu(l3)

	    output = tf.add(tf.matmul(l1,output_layer['weight']) , output_layer['bias'],name='z')

	    return output

	def train_neural_network(x):
		prediction = neural_network_model(x)
		cost = tf.reduce_mean( tf.nn.softmax_cross_entropy_with_logits(prediction,y) )
		optimizer = tf.train.AdamOptimizer(learning_rate=0.001).minimize(cost)

		with tf.Session() as sess:
			sess.run(tf.initialize_all_variables())	
		    
			for epoch in range(hm_epochs):
				epoch_loss = 0
				i=0
				while i < len(train_x):
					start = i
					end = i+batch_size
					batch_x = train_x[start:end]
					batch_y = train_y[start:end]

					_, c = sess.run([optimizer, cost], feed_dict={x: batch_x,
					                                              y: batch_y})
					epoch_loss += c
					i+=batch_size
					
				print('Epoch', epoch+1, 'completed out of',hm_epochs,'loss:',epoch_loss)
				
				#correct = tf.equal(tf.argmax(prediction, 1), tf.argmax(y, 1))
				#accuracy = tf.reduce_mean(tf.cast(correct, 'float'))
				#print('Accuracy:',accuracy.eval({x:test_x, y:test_y}))
				#print("--- %s seconds ---" % (time.time() - start_time))
			saver.save(cur_dir+'/saved_weights')
			correct = tf.equal(tf.argmax(prediction, 1), tf.argmax(y, 1))
			accuracy = tf.reduce_mean(tf.cast(correct, 'float'))
			y_p = tf.argmax(prediction, 1)
			all_y_pred = []
			print("between all_y_pred and j")
			j=0
			while j < len(test_x):
				start_test = j
				end_test = j + batch_size
				batch_test_x = test_x[start_test:end_test]
				batch_test_y = test_y[start_test:end_test]

				val_accuracy, y_pred = sess.run([accuracy, y_p], feed_dict={x:batch_test_x, y:batch_test_y})
				all_y_pred = np.concatenate([all_y_pred,y_pred])
				j+=batch_size

			#print("all y pred:",all_y_pred[0])
			#print('Accuracy:',val_accuracy)
			#print('Accuracy:',accuracy.eval({x:test_x, y:test_y}))
			y_true = np.argmax(test_y, 1)
			print('Accuracy:', sk.metrics.accuracy_score(y_true, all_y_pred))
			print('Precision:', sk.metrics.precision_score(y_true, all_y_pred))
			print('Recall:', sk.metrics.recall_score(y_true, all_y_pred))
			print('f1_score:', sk.metrics.f1_score(y_true, all_y_pred))
			print ('confusion_matrix:')
			print (sk.metrics.confusion_matrix(y_true, all_y_pred))
			print(y_true)
			print(all_y_pred)
			fpr, tpr, thresholds = sk.metrics.roc_curve(y_true, all_y_pred)
			print ( sk.metrics.roc_auc_score(y_true, all_y_pred))

		    
	train_neural_network(x)
	print("--- %s seconds ---" % (time.time() - start_time))
	return 'saved_weights'

