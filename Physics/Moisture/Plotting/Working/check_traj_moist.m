close all
clear
clc

cd /discover/nobackup/drholdaw/tmp.20188

TS1 = ncread('iau590a.traj.M2.lcv.20120317_0000z.nc4','TS');
TS2 = ncread('iau590a.traj.M2.lcv.20120317_0000z.nc4','TS1');

cd /home/drholdaw/Lin_Moist_Physics/DAS_plotting/

contourf(TS1-TS2)