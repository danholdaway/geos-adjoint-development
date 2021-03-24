close all
clear
clc
mydir = pwd;

cd ../

tracersphereend;

figure
contourf(q')
colorbar

figure
contourf(qlrlw')
colorbar

figure
contourf(qlrlwlim')
colorbar

cd(mydir); clear mydir