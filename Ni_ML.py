import numpy as np
from sklearn import svm
from readinput import *
import random

#extract feature vector, lables and cooridnates from the input file.
# these information will be used to ML model 
def create_input_data(filename,outfile):
    Natoms,nfeature,position,property2,features=readfile(filename)
    print("Natoms  and number of features per atom",Natoms,nfeature)
    np.save(outfile+"_XX",features);
    np.save(outfile+"_YY",property2);
    total_atoms = np.concatenate((position,property2),axis=1)
    np.save(outfile+"_pos",total_atoms);


class build_classifier:
    def __init__(self,trainX,trainY):
        self.training_X = np.load(trainX)
        self.training_Y = np.load(trainY)
        self.meanval = self.training_X.mean(axis=0)
        self.stdval = self.training_X.std(axis=0)
        self.ntraining_X = (self.training_X-self.meanval)/self.stdval
        self.training_Y = self.training_Y.ravel()
    
    def train(self):
        atom_frac = [np.sum(self.training_Y == val) for val in np.unique(self.training_Y)]
        atom_select=0.1+atom_frac[2]/atom_frac
        #select training examples
        select_atom=[]
        bulk=0;surface=0;defect=0
        for ii in range(0,len(self.training_Y)):
            n_rand = random.uniform(0, 1)
            if self.training_Y[ii] == 0 and n_rand <= atom_select[0]:
                select_atom.append(True)
                bulk+=1
            elif self.training_Y[ii] == 1 and n_rand <= atom_select[1]:
                select_atom.append(True)
                surface+=1
            elif self.training_Y[ii] == 2 and n_rand <= atom_select[2]:
                select_atom.append(True)
                defect+=1
            else:
                select_atom.append(False)
        select_atom = np.array(select_atom, dtype = bool)
        utraining_Y=self.training_Y[select_atom]
        u_ntrain=self.ntraining_X[select_atom]
        print("Number of training examples: ",len(u_ntrain))

        #SVM classifier
        clf = svm.LinearSVC(penalty='l1',loss='squared_hinge', dual=False,C=0.9,
                    tol=0.00001,max_iter=50000).fit(u_ntrain, utraining_Y)
        #Training Accuracy
        y_label = clf.predict(u_ntrain)
        yy_ = self.accuracy(utraining_Y,y_label)
        print("training error: ",np.mean(yy_)*100.0)
        print("training accuracy: ",100.00-np.mean(yy_)*100.0)
        self.model = clf
    
    def predict(self,input_X):
        n_X = (input_X-self.meanval)/self.stdval
        y_label = self.model.predict(n_X)
        return y_label
    
    def accuracy(self,y_true,y_predict):
        yy_ = []
        yy_.append(1 - np.mean(y_true == y_predict))
        return yy_





