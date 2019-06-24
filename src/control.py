from os import system as os
import sys
import numpy as np
    

#Compile all code
os("gfortran -ffree-form -ffree-line-length-none -fdefault-real-8 -O3 *.f -J mod/")
os("gfortran -ffree-form -ffree-line-length-none -fdefault-real-8 -O3 *.f -J mod/")
os("gfortran -ffree-form -ffree-line-length-none -fdefault-real-8 -O3 *.f -J mod/")
    




def Feed_R(r):
    os("./a.out " +r)
    #os("python ../Tools/SNR.py ")








def process():
    #PSR Orbital dynamics
    os("./a.out")

#process()










r = '6.0'

Feed_R(str(r))
sys.exit()

os('rm SNR_Data.txt')
os('touch SNR_Data.txt')


rs = np.linspace(6.0,8.0,20.0)
for r in rs:
    Feed_R(str(r))



