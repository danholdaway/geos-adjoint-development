close all
clear
clc
mydir = pwd;


dir1 = '/discover/nobackup/drholdaw/btmp.25956/hord_tr/hord_tr1/prog_free/';
dir2 = '/discover/nobackup/drholdaw/btmp.25956/';
dir3 = '/discover/nobackup/drholdaw/btmp.25956/hord_tr/hord_tr1/prog_free/';

file1 = 'v000_C180.prog.eta.20140201_0600z.nc4';
file2 = 'v000_C180.prog.eta.20140201_0600z.nc4';
file3 = 'v000_C180.prog.eta.20140201_0600z.nc4';

cd(dir1)
file = file1;

lon = ncread(file,'lon');
lat = ncread(file,'lat');
lev = ncread(file,'lev');

u_1 = ncread(file,'u');
v_1 = ncread(file,'v');
t_1 = ncread(file,'tv');
q_1 = ncread(file,'sphu');
p_1 = ncread(file,'delp');
qi_1 = ncread(file,'qitot');
ql_1 = ncread(file,'qltot');
o3_1 = ncread(file,'ozone');

cd(dir2)
file = file2;

u_2 = ncread(file,'u');
v_2 = ncread(file,'v');
t_2 = ncread(file,'tv');
q_2 = ncread(file,'sphu');
p_2 = ncread(file,'delp');
qi_2 = ncread(file,'qitot');
ql_2 = ncread(file,'qltot');
o3_2 = ncread(file,'ozone');

cd(dir3)
file = file3;

u_3 = ncread(file,'u');
v_3 = ncread(file,'v');
t_3 = ncread(file,'tv');
q_3 = ncread(file,'sphu');
p_3 = ncread(file,'delp');
qi_3 = ncread(file,'qitot');
ql_3 = ncread(file,'qltot');
o3_3 = ncread(file,'ozone');

cd(mydir)


x_all = [u_1-u_2; v_1-v_2; t_1-t_2; q_1-q_2; p_1-p_2; o3_1-o3_2];
x_all = [qi_1-qi_2; ql_1-ql_2;];

max(abs(x_all(:)))

asd

plot_lev = 62;

a = (ql_1(:,:,plot_lev)'-ql_3(:,:,plot_lev)');
b = (ql_2(:,:,plot_lev)'-ql_3(:,:,plot_lev)');

figure
contour(lon,lat,a)
colorbar

figure
contour(lon,lat,0.1*b)
colorbar

figure
contour(lon,lat,a-0.1*b)
colorbar

