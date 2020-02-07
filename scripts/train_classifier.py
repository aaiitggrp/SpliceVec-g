#from create_sentiment_featuresets import create_feature_sets_and_labels
import tensorflow as tf
#from tensorflow.examples.tutorials.mnist import input_data
import pickle
import time
import numpy as np
import sklearn as sk
from sklearn.model_selection import train_test_split
import os
start_time = time.time()
train_x=np.load("all-chr_3mer_non-coding_withn(250_5_5_8_hs_iter20)_v20_train_40stream_intronic_type2.nparray")
test_x=np.load("all-chr_3mer_non-coding_withn(250_5_5_8_hs_iter20)_v26_test_40stream_intronic_type2.nparray")
labels=np.array([0,1])
labels=np.tile(labels,(290502,1))
l21=np.array([1,0])
l21=np.tile(l21,(290502,1))
train_y=np.concatenate((labels,l21))

labels=np.array([0,1])
labels=np.tile(labels,(5612,1))
l21=np.array([1,0])
l21=np.tile(l21,(5612,1))
test_y=np.concatenate((labels,l21))

print(train_x.shape)
print(train_y.shape)
print(test_x.shape)
print(test_y.shape)
# Issue was the inorrect order of variables for train test split
train_x,_,train_y,_ = train_test_split(train_x, train_y, test_size=0, random_state=0)
print(train_x.shape)
print(train_y.shape)
n_nodes_hl1 = 2500
#n_nodes_hl2 = 2500
#n_nodes_hl3 = 2500

n_classes = 2
batch_size = 128
hm_epochs = 50

x = tf.placeholder('float')
y = tf.placeholder('float')

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

    output = tf.matmul(l1,output_layer['weight']) + output_layer['bias']

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
