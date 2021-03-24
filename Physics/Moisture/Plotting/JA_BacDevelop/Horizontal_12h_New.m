close all
clear
clc
mydir = pwd;


%Load MATLAB topography
load coast
coast_lat = lat; clear lat
coast_lon = long; clear long

%Choose whether you want plots of pert traj (0/1) and which model level to plot 
makeplots = 1;

plot_level_ql = 63;
plot_level_qi = 44;

% 0 to Scale plots individually
% 1 to Scale plots by nonlinear
scale = 1; 

%Make center of colormap white
cmap = colormap;
cmap(32,:) = [1 1 1];
cmap(33,:) = [1 1 1];
close

%Some plotting things
fontsize = 9;
line_wid_cont = 1.0;
line_wid_det = 0.6;
grey = 0.75;

%Choose the time being considered: 16th at 12z would be 6_12
file_end = '1_1500';

im = 576;
jm = 361;
lm = 72;

LevMin = 1;
LevMax = lm;

cd /home/drholdaw/LinearisedPhysics/Inputs/
pref
 
%Load Free (background) State.
cd /discover/nobackup/drholdaw/ExperimentData/BacmeisterPaper/ModelOutput/prog/free/
file = ['v000_C180.prog.eta.2014020',file_end,'z.nc4'];

lon = ncread(file,'lon'); lon = round(lon*10000)/10000;
lat = ncread(file,'lat');

qi_free = ncread(file,'qitot',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
ql_free = ncread(file,'qltot',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);

%Load Perturbed (analysis) state
cd /discover/nobackup/drholdaw/ExperimentData/BacmeisterPaper/ModelOutput/prog/replay/
file = ['v000_C180.prog.eta.2014020',file_end,'z.nc4'];

qi_replay = ncread(file,'qitot',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
ql_replay = ncread(file,'qltot',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);

%Load TLM state to compare.
cd /discover/nobackup/drholdaw/ExperimentData/BacmeisterPaper/ModelOutput/4DVAR
file = ['v000_C180.fvpert.eta.2014020',file_end,'z_ST03z_P21000_S1.nc4'];

qi_tlm = ncread(file,'QI',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
ql_tlm = ncread(file,'QL',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);

cd(mydir)

%Compute NL perturbation trajectory.
qi_nlm = qi_replay - qi_free;
ql_nlm = ql_replay - ql_free;


%Global Plot
ilonmin0 = find(lon==-180);
ilonmax0 = find(lon==179.375);
ilatmin0 = find(lat==-90);
ilatmax0 = find(lat==90);

%South Indian Plot
ilonmin1 = find(lon==50);
ilonmax1 = find(lon==150);
ilatmin1 = find(lat==-70);
ilatmax1 = find(lat==-30);

%South Pacific Plot
ilonmin2 = find(lon==-140);
ilonmax2 = find(lon==-60);
ilatmin2 = find(lat==-70);
ilatmax2 = find(lat==-28);

%North Atlantic Plot
ilonmin3 = find(lon==-70);
ilonmax3 = find(lon==0);
ilatmin3 = find(lat==20);
ilatmax3 = find(lat==50);


%First plot 

ilonmin = ilonmin3;
ilonmax = ilonmax3;
ilatmin = ilatmin3;
ilatmax = ilatmax3;


%Compute maximum absolute values for contour limits and intervals.
maxql = 7e-4;
cint_ql = 4e-5;



%PLOT FOR QL
figure
set(gcf,'position',[1320 15 878 450])


%Second plot 

ilonmin = ilonmin1;
ilonmax = ilonmax1;
ilatmin = ilatmin1;
ilatmax = ilatmax1;

ql_a_lev = ql_nlm(ilonmin:ilonmax,ilatmin:ilatmax,plot_level_ql);
ql_b_lev = ql_tlm(ilonmin:ilonmax,ilatmin:ilatmax,plot_level_ql);


subplot(2,2,1)
contour(lon(ilonmin:ilonmax),lat(ilatmin:ilatmax),ql_a_lev','LineWidth',line_wid_cont,'LevelStep',cint_ql);
hold on;
plot(coast_lon,coast_lat,'Color',[grey grey grey])
caxis([-maxql maxql])
title('(a) q_l Nonlinear Perturbation at 850hPa','FontSize',fontsize,'FontName','TimesNewRoman');
ylabel('Latitude (\circ)','FontSize',fontsize,'FontName','TimesNewRoman')
xlabel('Longitude (\circ)','FontSize',fontsize,'FontName','TimesNewRoman')
a = get(gca,'position');
set(gca,'position',[0.1 a(2) a(3) a(4)])
set(gca,'XTick',[50 70 90 110 130 150])
colormap(cmap)
box on
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')


subplot(2,2,2)
contour(lon(ilonmin:ilonmax),lat(ilatmin:ilatmax),ql_b_lev','LineWidth',line_wid_cont,'LevelStep',cint_ql);
hold on;
plot(coast_lon,coast_lat,'Color',[grey grey grey])
caxis([-maxql maxql])
title('(b) q_l Tangent Linear Perturbation at 850hPa','FontSize',fontsize,'FontName','TimesNewRoman');
ylabel('Latitude (\circ)','FontSize',fontsize,'FontName','TimesNewRoman')
xlabel('Longitude (\circ)','FontSize',fontsize,'FontName','TimesNewRoman')
colormap(cmap)
cbar = colorbar;
b = get(gca,'position');
set(gca,'position',[0.475 b(2) a(3) b(4)])
c = get(cbar,'position');
set(cbar,'position',[0.87 c(2) c(3) 0.35])
set(gca,'XTick',[50 70 90 110 130 150])
set(get(cbar,'ylabel'),'String', '(kgkg^{-1})','FontSize',fontsize,'FontName','TimesNewRoman')

set(gca,'YAxisLocation','right')
box on   
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Australasia Plot
ilonmin1 = find(lon==100);
ilonmax1 = find(lon==max(lon));
ilatmin1 = find(lat==-40);
ilatmax1 = find(lat==20);

%Brazil Plot
ilonmin2 = find(lon==-80);
ilonmax2 = find(lon==-35);
ilatmin2 = find(lat==-25);
ilatmax2 = find(lat==7.5);

%First plot 

ilonmin = ilonmin1
ilonmax = ilonmax1
ilatmin = ilatmin1
ilatmax = ilatmax1

qi_a_lev = qi_nlm(ilonmin:ilonmax,ilatmin:ilatmax,plot_level_qi);
qi_b_lev = qi_tlm(ilonmin:ilonmax,ilatmin:ilatmax,plot_level_qi);

%Compute maximum absolute values for contour limits and intervals.
maxqi = 3.0e-4;
cint_qi = maxqi/15;





%Second plot 

ilonmin = ilonmin2;
ilonmax = ilonmax2;
ilatmin = ilatmin2;
ilatmax = ilatmax2;

qi_a_lev = qi_nlm(ilonmin:ilonmax,ilatmin:ilatmax,plot_level_qi);
qi_b_lev = qi_tlm(ilonmin:ilonmax,ilatmin:ilatmax,plot_level_qi);


subplot(2,2,3)
contour(lon(ilonmin:ilonmax),lat(ilatmin:ilatmax),qi_a_lev','LineWidth',line_wid_cont,'LevelStep',cint_qi);
hold on;
plot(coast_lon,coast_lat,'Color',[grey grey grey])
caxis([-maxqi maxqi])
title('(c) q_i Nonlinear Perturbation at 250hPa','FontSize',fontsize,'FontName','TimesNewRoman');
ylabel('Latitude (\circ)','FontSize',fontsize,'FontName','TimesNewRoman')
xlabel('Longitude (\circ)','FontSize',fontsize,'FontName','TimesNewRoman')
a = get(gca,'position');
set(gca,'position',[0.1 a(2) a(3) a(4)])
colormap(cmap)
box on
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')


subplot(2,2,4)
contour(lon(ilonmin:ilonmax),lat(ilatmin:ilatmax),qi_b_lev','LineWidth',line_wid_cont,'LevelStep',cint_qi);
hold on;
plot(coast_lon,coast_lat,'Color',[grey grey grey])
caxis([-maxqi maxqi])
title('(d) q_i Tangent Linear Perturbation at 250hPa','FontSize',fontsize,'FontName','TimesNewRoman');
ylabel('Latitude (\circ)','FontSize',fontsize,'FontName','TimesNewRoman')
xlabel('Longitude (\circ)','FontSize',fontsize,'FontName','TimesNewRoman')
colormap(cmap)
cbar = colorbar;
b = get(gca,'position');
set(gca,'position',[0.475 b(2) a(3) b(4)])
c = get(cbar,'position');
set(cbar,'position',[0.87 c(2) c(3) 0.35])
set(get(cbar,'ylabel'),'String', '(kgkg^{-1})','FontSize',fontsize,'FontName','TimesNewRoman')

set(gca,'YAxisLocation','right')
box on   
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')



