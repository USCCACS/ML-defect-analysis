import numpy as np

def readfile(filename):
    with open(filename,'r') as inputfile:
        Natoms=int(inputfile.readline().strip())
        nfeature = int(inputfile.readline().strip())
        boxsize = inputfile.readline()
        position = np.empty((Natoms, 3), dtype='float32')
        property2 = np.empty((Natoms,1), dtype='float32')
        features = np.empty((Natoms,nfeature),dtype='float32')
        ii=0
        for val in inputfile:
            val=val.strip().split()
            position[ii,:] = np.array(val[1:4]).astype('float32')
            property2[ii] = np.array(val[4]).astype('float32')
            features[ii,:] = np.array(val[5:]).astype('float32')
            ii +=1

    return Natoms,nfeature,position,property2,features


def writexyz(Natoms,position,predict,ground_t):

    with open('output.xyz','w') as outputfile:
        outputfile.write(str(Natoms) + "\n")
        outputfile.write(str(Natoms) + "\n")
        for ii in range(0,Natoms):
            outputfile.write("Ni %12.6f %12.6f %12.6f  %6d  %6d\n" % (position[ii][0], position[ii][1], position[ii][2],
                                                                 predict[ii],ground_t[ii]))











