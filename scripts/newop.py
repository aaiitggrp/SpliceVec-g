import numpy as np
import tensorflow as tf
import os 
import gensim 
import sys
cur_dir=os.getcwd()
sys.path.insert(0,cur_dir+'/scripts')
from vector_for_5th import get_vectors_for_single_file
from get_vectors_me import generate_vectors
def test_on_new_data(filename,vec_choice,model_name,saved_weights_name):
	cur_dir=os.getcwd()
	new_g=tf.train.import_meta_graph(cur_dir+'/'+saved_weights_name+'.meta')
	sess=tf.Session()
	new_g.restore(sess,tf.train.latest_checkpoint('./'))
	cx=tf.get_default_graph()
	x=cx.get_tensor_by_name('X:0')
	y=cx.get_tensor_by_name('Y:0')
	z=cx.get_tensor_by_name('z:0')
	pred_temp = tf.equal(tf.argmax(z, 1), tf.argmax(y, 1))
	accuracy = tf.reduce_mean(tf.cast(pred_temp, "float"))
	if len(filename)==2:
		len_of_files,output_files=generate_vectors(model_name,filename,vec_choice)
		test_x=np.load(cur_dir+"/corpus/"+output_files[1])
		labels=np.array([0,1])
		labels=np.tile(labels,(len_of_files[0],1))
		l21=np.array([0,1])
		l21=np.tile(l21,(len_of_files[1],1))
		test_y=np.concatenate((labels,l21))
		print ("testing  Accuracy:", accuracy.eval({x:test_x, y:test_y},session=sess))
	else:
		test_x=get_vectors_for_single_file(filename[0],vec_choice,model_name)
		pred_y=tf.argmax(z,1)
		print ("predicted value is:",pred_y.eval({x:test_x},session=sess))

	print ("------------------prediction complete------------------------------")



