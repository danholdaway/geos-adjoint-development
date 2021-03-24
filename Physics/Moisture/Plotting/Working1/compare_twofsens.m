close all
clear
clc
mydir = pwd;


dir1 = '/discover/nobackup/drholdaw/atmp.22293/sens.20130117.000000';
dir2 = '/discover/nobackup/drholdaw/atmp.22293/sens.20130117.000000';

% dir1 = '/archive/u/drholdaw/u001_C180_PH1/prog/Y2013/M01/D08/H21/';
% dir2 = '/archive/u/drholdaw/u001_C180_PH2/prog/Y2013/M01/D08/H21/';

file1 = 'x0011dh_a.fvpert.eta.20130116_0600z_RADTEST_RADOFFM.nc4';
file2 = 'x0011dh_a.fvpert.eta.20130116_0600z_RADTEST_RADON.nc4';

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


cd(mydir)


x_all = [u_1-u_2; v_1-v_2; t_1-t_2; q_1-q_2; p_1-p_2; qi_1-qi_2; ql_1-ql_2; o3_1-o3_2];

max(abs(x_all(:)))

plot_lev = 63;

u_1l = u_1(:,:,plot_lev);
u_2l = u_2(:,:,plot_lev);
v_1l = v_1(:,:,plot_lev);
v_2l = v_2(:,:,plot_lev);
t_1l = t_1(:,:,plot_lev);
t_2l = t_2(:,:,plot_lev);
q_1l = q_1(:,:,plot_lev);
q_2l = q_2(:,:,plot_lev);
p_1l = p_1(:,:,plot_lev);
p_2l = p_2(:,:,plot_lev);
qi_1l = qi_1(:,:,plot_lev);
qi_2l = qi_2(:,:,plot_lev);
ql_1l = ql_1(:,:,plot_lev);
ql_2l = ql_2(:,:,plot_lev);
o3_1l = o3_1(:,:,plot_lev);
o3_2l = o3_2(:,:,plot_lev);

umax = max(abs([u_1l(:); u_1l(:)]));
ustep = 2*umax/20;
vmax = max(abs([v_1l(:); v_1l(:)]));
vstep = 2*vmax/20;
tmax = max(abs([t_1l(:); t_1l(:)]));
tstep = 2*tmax/20;
qmax = max(abs([q_1l(:); q_1l(:)]));
qstep = 2*qmax/20;
pmax = max(abs([p_1l(:); p_1l(:)]));
pstep = 2*pmax/20;
qimax = max(abs([qi_1l(:); qi_1l(:)]));
qistep = 2*qimax/20;
qlmax = max(abs([ql_1l(:); ql_1l(:)]));
qlstep = 2*qlmax/20;
o3max = max(abs([o3_1l(:); o3_1l(:)]));
o3step = 2*o3max/20;

figure
set(gcf,'position',[719    31   560   888])

subplot(3,1,1)
[~,h] = contourf(lon,lat,u_1l','LineStyle','none');
caxis([-umax umax])
% xlim([-65 -25])
% ylim([25 55])
colorbar
set(h,'LevelStep',ustep)

subplot(3,1,2)
[~,h] = contourf(lon,lat,u_2l','LineStyle','none');
caxis([-umax umax])
% xlim([-65 -25])
% ylim([25 55])
colorbar
set(h,'LevelStep',ustep)

subplot(3,1,3)
[~,h] = contourf(lon,lat,u_1l'-u_2l','LineStyle','none');
caxis([-umax umax])
% xlim([-65 -25])
% ylim([25 55])
colorbar
set(h,'LevelStep',ustep)


figure
set(gcf,'position',[719    31   560   888])

subplot(3,1,1)
[~,h] = contourf(lon,lat,v_1l','LineStyle','none');
caxis([-vmax vmax])
% xlim([-65 -25])
% ylim([25 55])
colorbar
set(h,'LevelStep',vstep)

subplot(3,1,2)
[~,h] = contourf(lon,lat,v_2l','LineStyle','none');
caxis([-vmax vmax])
% xlim([-65 -25])
% ylim([25 55])
colorbar
set(h,'LevelStep',vstep)

subplot(3,1,3)
[~,h] = contourf(lon,lat,v_1l'-v_2l','LineStyle','none');
caxis([-vmax vmax])
% xlim([-65 -25])
% ylim([25 55])
colorbar
set(h,'LevelStep',vstep)

figure
set(gcf,'position',[719    31   560   888])

subplot(3,1,1)
[~,h] = contourf(lon,lat,t_1l','LineStyle','none');
caxis([-tmax tmax])
% xlim([-65 -25])
% ylim([25 55])
colorbar
set(h,'LevelStep',tstep)

subplot(3,1,2)
[~,h] = contourf(lon,lat,t_2l','LineStyle','none');
caxis([-tmax tmax])
% xlim([-65 -25])
% ylim([25 55])
colorbar
set(h,'LevelStep',tstep)

subplot(3,1,3)
[~,h] = contourf(lon,lat,t_1l'-t_2l','LineStyle','none');
caxis([-tmax tmax])
% xlim([-65 -25])
% ylim([25 55])
colorbar
set(h,'LevelStep',tstep)

figure
set(gcf,'position',[719    31   560   888])

subplot(3,1,1)
[~,h] = contourf(lon,lat,q_1l','LineStyle','none');
caxis([-qmax qmax])
% xlim([-65 -25])
% ylim([25 55])
colorbar
set(h,'LevelStep',qstep)

subplot(3,1,2)
[~,h] = contourf(lon,lat,q_2l','LineStyle','none');
caxis([-qmax qmax])
% xlim([-65 -25])
% ylim([25 55])
colorbar
set(h,'LevelStep',qstep)

subplot(3,1,3)
[~,h] = contourf(lon,lat,q_1l'-q_2l','LineStyle','none');
caxis([-qmax qmax])
% xlim([-65 -25])
% ylim([25 55])
colorbar
set(h,'LevelStep',qstep)
