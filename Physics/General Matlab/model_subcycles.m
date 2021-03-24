close all
clear
clc
mydir = pwd;

DT = 450;
NPX = 360;

numsc_1 = 6* round (abs(DT/1800.)*((4.0*NPX)/180.0) + 0.49 )
numsc_2 = 7;

DT_sc = DT/numsc_1


DT = 450;
NPX = 720;

numsc_1 = 6* round (abs(DT/1800.)*((4.0*NPX)/180.0) + 0.49 )
numsc_2 = 7;

DT_sc = DT/numsc_1

cd(mydir)