close all
clear
clc
mydir = pwd;

cd /discover/nobackup/drholdaw/wrk.dust/sens.20060828.000000/

file = 'fvpert.eta.nc4';
u = ncread(file,'u');
v = ncread(file,'u');
t = ncread(file,'u');
q = ncread(file,'u');
p = ncread(file,'u');
ql = ncread(file,'u');
qi = ncread(file,'u');
o3 = ncread(file,'u');

u = 10e-4 * u./u;
v = 10e-4 * u./u;
t = 10e-4 * u./u;
q = 10e-4 * u./u;
p = 10e-4 * u./u;
ql = 10e-4 * u./u;
qi = 10e-4 * u./u;
o3 = 10e-4 * u./u;


file = 'fvpert1.eta.nc4';
ncwrite(file,'u',u);
ncwrite(file,'v',v);
ncwrite(file,'tv',t);
ncwrite(file,'sphu',q);
ncwrite(file,'delp',p);
ncwrite(file,'qltot',ql);
ncwrite(file,'qitot',qi);
ncwrite(file,'ozone',o3);

cd(mydir)