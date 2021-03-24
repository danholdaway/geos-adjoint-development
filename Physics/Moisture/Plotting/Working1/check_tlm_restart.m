close all
clear
clc
mydir = pwd;

cd /discover/nobackup/drholdaw/wrk.dust/sens.20060828.000000/Restarts_DotP_Simp/
file = 'fvpert.eta.nc4';

lon = ncread(file,'lon');
lat = ncread(file,'lat');
lev = ncread(file,'lev');

u_nom = ncread(file,'u');
v_nom = ncread(file,'v');
t_nom = ncread(file,'tv');
q_nom = ncread(file,'sphu');
p_nom = ncread(file,'delp');
qi_nom = ncread(file,'qitot');
ql_nom = ncread(file,'qltot');
o3_nom = ncread(file,'ozone');


cd /discover/nobackup/drholdaw/wrk.dust/sens.20060828.000000/Restarts_DotP_Simp/
file = 'Jgradf.eta.nc4';

u_wim = ncread(file,'u');
v_wim = ncread(file,'v');
t_wim = ncread(file,'tv');
q_wim = ncread(file,'sphu');
p_wim = ncread(file,'delp');
qi_wim = ncread(file,'qitot');
ql_wim = ncread(file,'qltot');
o3_wim = ncread(file,'ozone');

cd(mydir)



figure
contourf(lon,lat,u_nom(:,:,50)','LineStyle','none')

figure
contourf(lon,lat,u_wim(:,:,50)','LineStyle','none')