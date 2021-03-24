close all
clear
clc
mydir = pwd;


A_tlm = (1e-12*rand(3)) + [ 4 6 3 ; 9 3 13 ; 1 14 8];
A_adm = (1e-12*rand(3)) + A_tlm';

x_tlm = [21 ; 3 ; 9];
y_tlm = A_tlm*x_tlm;

y_adm = [18 ; 4 ; 9];
x_adm = A_adm*y_adm;

format long
dot1 = dot(y_tlm,y_adm)
dot2 = dot(x_tlm,x_adm)
format short

cd(mydir)