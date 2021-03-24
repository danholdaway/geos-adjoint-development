close all
clear
clc


cd /discover/nobackup/drholdaw/btmp.25956/JNorm/
file = 'Jgradf_twe.eta.nc4';

DU001 = ncread(file,'DU001');
DU002 = ncread(file,'DU002');
DU003 = ncread(file,'DU003');
DU004 = ncread(file,'DU004');
DU005 = ncread(file,'DU005');

DU001 = 0.0;
DU002 = 0.0;
DU003 = 0.0;
DU004 = 0.0;
DU005 = 0.0;

ncwrite(file,'DU001',DU001);
ncwrite(file,'DU002',DU002);
ncwrite(file,'DU003',DU003);
ncwrite(file,'DU004',DU004);
ncwrite(file,'DU005',DU005);


