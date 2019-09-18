from os import system as os

os("gfortran -fcheck=all -ffree-form -fdefault-real-8 -O3 parameters.f constants.f GW.f tensors.f derivatives.f initial_conditions.f rk.f main.f -J mod/") 
os("./a.out") 
