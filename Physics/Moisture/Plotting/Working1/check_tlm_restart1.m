close all
clear
clc

cd /discover/nobackup/drholdaw/btmp.25956/tlmrestarts_03z/


ql = ncread('fvpert.eta.nc4','qltot');
lon = ncread('fvpert.eta.nc4','lon');
lat = ncread('fvpert.eta.nc4','lat');

contourf(lon,lat,ql(:,:,63)')