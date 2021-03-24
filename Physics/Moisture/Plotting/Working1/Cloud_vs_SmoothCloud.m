close all
clear
clc

%Choose the time being considered: 16th at 12z would be 6_12
% file_end = '6_0330'; % 15 mins
% file_end = '6_0400'; % 1 hours
% file_end = '6_0600'; % 3 hours
% file_end = '6_0900'; % 6 hours
% file_end = '6_1200'; % 9 hours
file_end = '6_1500'; % 12 hours
% file_end = '6_1800'; % 15 hours
% file_end = '6_2100'; % 18 hours
% file_end = '7_0000'; % 21 hours
% file_end = '7_0300'; % 24 hours

im = 576;
jm = 361;
lm = 72;

LevMin = 1;
LevMax = lm;

%Load Free (background) State.
cd /discover/nobackup/drholdaw/atmp.22292/prog/prog_free/
file = ['x0011dh_a.prog.eta.2013011',file_end,'z.nc4'];

lon = ncread(file,'lon');
lat = ncread(file,'lat');
lev = ncread(file,'lev');

u = ncread(file,'u',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
v = ncread(file,'v',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
t = ncread(file,'tv',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
q = ncread(file,'sphu',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
qi = ncread(file,'qitot',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
ql = ncread(file,'qltot',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);

%Load Smoothed Free (background) State.
cd /discover/nobackup/drholdaw/atmp.22292/prog/prog_free_smoothcloud/
file = ['x0011dh_a.prog.eta.2013011',file_end,'z.nc4'];

lon = ncread(file,'lon');
lat = ncread(file,'lat');
lev = ncread(file,'lev');

u_s = ncread(file,'u',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
v_s = ncread(file,'v',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
t_s = ncread(file,'tv',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
q_s = ncread(file,'sphu',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
qi_s = ncread(file,'qitot',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
ql_s = ncread(file,'qltot',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);



plot_level_l = 64;
plot_level_i = 44;

u_l = u(:,:,plot_level_l);
v_l = v(:,:,plot_level_l);
t_l = t(:,:,plot_level_l);
q_l = q(:,:,plot_level_l);
qi_l = qi(:,:,plot_level_i);
ql_l = ql(:,:,plot_level_l);


u_s_l = u_s(:,:,plot_level_l);
v_s_l = v_s(:,:,plot_level_l);
t_s_l = t_s(:,:,plot_level_l);
q_s_l = q_s(:,:,plot_level_l);
qi_s_l = qi_s(:,:,plot_level_i);
ql_s_l = ql_s(:,:,plot_level_l);


ql_all = [ql_l ; ql_s_l; ql_l-ql_s_l];
max_ql = max(abs(ql_all(:)));

qi_all = [qi_l ; qi_s_l; qi_l-qi_s_l];
max_qi = max(abs(qi_all(:)));


figure
set(gcf,'position',[719    52   560   867])

subplot(3,1,1)
contour(lon,lat,ql_l')
colorbar
caxis([0 max_ql])

subplot(3,1,2)
contour(lon,lat,ql_s_l')
colorbar
caxis([0 max_ql])

subplot(3,1,3)
contourf(lon,lat,ql_l'-ql_s_l','LineStyle','none')
colorbar
caxis([-max_ql max_ql])



figure
set(gcf,'position',[719    52   560   867])

subplot(3,1,1)
contour(lon,lat,qi_l')
colorbar
caxis([0 max_qi])

subplot(3,1,2)
contour(lon,lat,qi_s_l')
colorbar
caxis([0 max_qi])

subplot(3,1,3)
contourf(lon,lat,qi_l'-qi_s_l','LineStyle','none')
colorbar
caxis([-max_qi max_qi])


















