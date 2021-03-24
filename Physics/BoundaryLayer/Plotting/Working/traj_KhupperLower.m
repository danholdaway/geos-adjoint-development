close all
clc
clear

cd /discover/nobackup/drholdaw/tmp.22292/
file = 'x0011dh_a.traj.lcv.20130116_0200z.nc4';

lon = ncread(file,'lon');
lat = ncread(file,'lat');

KhL = ncread(file,'KHLOWER');
KhU = ncread(file,'KHUPPER');

diff = KhU - KhL