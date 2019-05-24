from os import system as os
import sys
import numpy as np
    

#Compile all code
os("gfortran -ffree-form -ffree-line-length-none -fdefault-real-8 -O3 *.f")
os("gfortran -ffree-form -ffree-line-length-none -fdefault-real-8 -O3 *.f")
os("gfortran -ffree-form -ffree-line-length-none -fdefault-real-8 -O3 *.f")
    
def process():
    #PSR Orbital dynamics
    os("./a.out")

process()




