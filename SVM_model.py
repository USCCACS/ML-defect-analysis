import sys
import numpy as np
from sklearn import svm
from Ni_ML import *

train_file = sys.argv[1]
test_file = sys.argv[2]
print(train_file,test_file)
#Train ML Model
inputfile=train_file
outfile='feature/train'
#extract feature vector, lables and cooridnates from the input file.
create_input_data(inputfile,outfile)
#build classifier
Ni_model = build_classifier('feature/train_XX.npy','feature/train_YY.npy')
Ni_model.train()


#Test ML Model
#Step 1: extract feature vector, lables and cooridnates from the input file.
inputfile=test_file
outfile='feature/test'
create_input_data(inputfile,outfile)
test_X = np.load('feature/test_XX.npy')
test_Y = np.load('feature/test_YY.npy')
test_pos = np.load('feature/test_pos.npy')
test_Y = test_Y.ravel()
#Step 2: Precit lables of the test data using the build classier
test_1_predict = Ni_model.predict(test_X)
test_1_accuracy = Ni_model.accuracy(test_Y,test_1_predict)
print("Test error: ",np.mean(test_1_accuracy)*100.0)
print("Test accuracy: ",100.00-np.mean(test_1_accuracy)*100.0)
writexyz(len(test_Y),test_pos,test_1_predict,test_Y)
