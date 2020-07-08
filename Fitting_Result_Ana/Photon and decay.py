import h5py
import numpy as np
import os

base_path = 'Z:\\Chenghang\\7.7.20.D2OTesting\\D2O\\'
read_path = base_path + 'acquisition\\bins\\'
if not os.path.exists(base_path + 'Fitting_result\\'):
    os.mkdir(base_path + 'Fitting_result\\')
out_path = base_path + 'Fitting_result\\'

name = '488storm_0002'

txtname = name + '.txt'
Filename = read_path + name + '_mlist.hdf5'

f1 = h5py.File(Filename,'r+')
f2 = open(out_path + txtname,"w+")

Frame = 0
Molecule = 0
Real_frame = 0

for key in f1.keys():
    print('Writing' + ' ' + str(key))
    if key == 'metadata.xml':
        break
    temp = str(key)
    Real_frame = int(temp[3:])
    length_temp = f1[key]['x'].shape
    length = length_temp[0]
    for i in range(length):
        a = str(f1[key]['background'][i])
        b = str(f1[key]['height'][i])
        c = str(f1[key]['sum'][i])
        d = str(Real_frame)
        L = a + ' ' + b + ' ' + c + ' ' + d + '\n'
        Molecule = Molecule + 1
        f2.write(L)
    Frame = Frame + 1

f1.close()
f2.close()
First_L = str(Frame) + ' ' + str(Molecule) + '\n'
f2 = open(out_path + txtname,"r")
oline = f2.readlines()
oline.insert(0,First_L)

f2.close()
f2 = open(out_path + txtname,"w")
f2.writelines(oline)
f2.close()
