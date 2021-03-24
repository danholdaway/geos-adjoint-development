close all
clear
clc
mydir = pwd;


A_tlm = (1e-14*rand(3)) + [ 4+2i 6 3 ; 9-9i 3 13+7i ; 1 14+3i 8];
A_adm = (1e-14*rand(3)) + A_tlm';

x_tlm = [21+4i ; 3-2i ; 9-7i];
y_tlm = A_tlm*x_tlm;

y_adm = [18-4i ; 4-8i ; 9+13i];
x_adm = A_adm*y_adm;

format long
dot1 = dot(y_tlm,y_adm)
dot2 = dot(x_tlm,x_adm)
format short

cd(mydir)