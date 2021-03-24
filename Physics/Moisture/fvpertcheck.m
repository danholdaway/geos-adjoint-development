close all
clear
clc

cd /discover/nobackup/drholdaw/btmp.25956/sens.20140202.000000/

file = 'fvpert.eta.nc4';

qi = ncread(file,'qitot');
ql = ncread(file,'qltot');