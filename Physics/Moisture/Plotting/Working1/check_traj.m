close all
clear
clc


cd /discover/nobackup/drholdaw/btmp.25956/traj/traj_free
file = 'v000_C180.traj.lcv.20140201_0145z.nc4';


DELT = ncread(file,'DELTIRD');

contourf(DELT)