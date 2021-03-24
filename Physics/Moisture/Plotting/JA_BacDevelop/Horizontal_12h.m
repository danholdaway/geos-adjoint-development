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

lon = ncread(file,'lon');
lat = ncread(file,'lat');
lev = ncread(file,'lev');

qi_free = ncread(file,'qitot',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
ql_free = ncread(file,'qltot',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);

%Load Perturbed (analysis) st
cd /discover/nobackup/drholdaw/ExperimentData/BacmeisterPaper/ModelOutput/prog/replay/
file = ['v000_C180.prog.eta.2014020',file_end,'z.nc4'];

qi_replay = ncread(file,'qitot',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
ql_replay = ncread(file,'qltot',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);

%Load TLM state to compare.
cd /discover/nobackup/drholdaw/wrk.bac/sens.20140202.000000/
file = ['v000_C180.fvpert.eta.2014020',file_end,'z_ST03z_P21000_S1.nc4'];

qi_tlm = ncread(file,'QI',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
ql_tlm = ncread(file,'QL',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);

cd(mydir)

%Compute NL perturbation trajectory.
qi_nlm = qi_replay - qi_free;
ql_nlm = ql_replay - ql_free;

%Pick out just the level to be plotted.
ql_a_lev = ql_nlm(:,:,plot_level_ql);

%Uncomment here to plot actual TLM pert trajectories
ql_b_lev = ql_tlm(:,:,plot_level_ql);

% %Uncomment here to plot the difference between the NL and TLM pert trajectories
% t_b_lev = t_tlm(:,:,plot_level_ql) - t_nlm(:,:,plot_level_ql);
% q_b_lev = q_tlm(:,:,plot_level_ql) - q_nlm(:,:,plot_level_ql);
% ql_b_lev = ql_tlm(:,:,plot_level_ql) - ql_nlm(:,:,plot_level_ql);

%Compute maximum absolute values for contour limits and intervals.
maxql = max(abs( [ql_a_lev(:); ql_b_lev(:)] ));

cint = 20;

%PLOT FOR QL
figure
set(gcf,'position',[1320 596 878 305])

%ql 1458         596         878         281

if scale == 0 || scale == 1
    maxql = max(abs(ql_a_lev(:)));
end
cint_ql = 2*maxql/cint;

subplot(2,2,1)
contour(lon,lat,ql_a_lev','LineWidth',line_wid_cont,'LevelStep',cint_ql);
hold on;
plot(coast_lon,coast_lat,'Color',[grey grey grey])
caxis([-maxql maxql])

ti = title('(a) q_l Nonlinear Perturbation at 850hPa','FontSize',fontsize,'FontName','TimesNewRoman');
ylabel('Latitude (\circ)','FontSize',fontsize,'FontName','TimesNewRoman')
a = get(gca,'position');
set(gca,'position',[0.1 a(2) a(3) a(4)])
d = get(ti,'position');
set(ti,'position',[d(1) d(2)-15 d(3)]);
colormap(cmap)
box on
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

if scale == 0
    maxql = max(abs(ql_c_lev(:)));
end
cint_ql = 2*maxql/cint;

if scale == 1
    for i = 1:im
        for j = 1:jm
            if abs(ql_b_lev(i,j)) >= maxql
                ql_b_lev(i,j) = maxql * abs(ql_b_lev(i,j))/ql_b_lev(i,j);
            end
        end
    end         
end

subplot(2,2,2)
contour(lon,lat,ql_b_lev','LineWidth',line_wid_cont,'LevelStep',cint_ql);
hold on;
plot(coast_lon,coast_lat,'Color',[grey grey grey])
caxis([-maxql maxql])

ti = title('(b) q_l Tangent Linear Perturbation at 850hPa','FontSize',fontsize,'FontName','TimesNewRoman');
colormap(cmap)
cbar = colorbar;
b = get(gca,'position');
set(gca,'position',[0.475 b(2) a(3) b(4)])
c = get(cbar,'position');
set(cbar,'position',[0.85 c(2) c(3) c(4)])
d = get(ti,'position');
set(ti,'position',[d(1) d(2)-15 d(3)]);
set(get(cbar,'ylabel'),'String', '(kgkg^{-1})','FontSize',fontsize,'FontName','TimesNewRoman')

set(gca,'YAxisLocation','right')
box on   
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')


%Pick out just the level to be plotted.
qi_a_lev = qi_nlm(:,:,plot_level_qi);

%Uncomment here to plot actual TLM pert trajectories
qi_b_lev = qi_tlm(:,:,plot_level_qi);

%Compute maximum absolute values for contour limits and intervals.
maxqi = max(abs( [qi_a_lev(:); qi_b_lev(:)] ));

if scale == 0 || scale == 1
    maxqi = max(abs(qi_a_lev(:)));
end
cint_qi = 2*maxqi/cint;

subplot(2,2,3)
contour(lon,lat,qi_a_lev','LineWidth',line_wid_cont,'LevelStep',cint_qi);
hold on;
plot(coast_lon,coast_lat,'Color',[grey grey grey])
caxis([-maxqi maxqi])

ti = title('(c) q_i Nonlinear Perturbation at 250hPa','FontSize',fontsize,'FontName','TimesNewRoman');
xlabel('Longitude (\circ)','FontSize',fontsize,'FontName','TimesNewRoman')
ylabel('Latitude (\circ)','FontSize',fontsize,'FontName','TimesNewRoman')
a = get(gca,'position');
set(gca,'position',[0.1 a(2) a(3) a(4)])
d = get(ti,'position');
set(ti,'position',[d(1) d(2)-15 d(3)]);
colormap(cmap)
box on
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

if scale == 0
    maxqi = max(abs(qi_c_lev(:)));
end
cint_qi = 2*maxqi/cint;

if scale == 1
    for i = 1:im
        for j = 1:jm
            if abs(qi_b_lev(i,j)) >= maxqi
                qi_b_lev(i,j) = maxqi * abs(qi_b_lev(i,j))/qi_b_lev(i,j);
            end
        end
    end         
end

subplot(2,2,4)
contour(lon,lat,qi_b_lev','LineWidth',line_wid_cont,'LevelStep',cint_qi);
hold on;
plot(coast_lon,coast_lat,'Color',[grey grey grey])
caxis([-maxqi maxqi])

ti = title('(d) q_i Tangent Linear Perturbation at 250hPa','FontSize',fontsize,'FontName','TimesNewRoman');
xlabel('Longitude (\circ)','FontSize',fontsize,'FontName','TimesNewRoman')
colormap(cmap)
cbar = colorbar;
b = get(gca,'position');
set(gca,'position',[0.475 b(2) a(3) b(4)])
c = get(cbar,'position');
set(cbar,'position',[0.85 c(2) c(3) c(4)])
d = get(ti,'position');
set(ti,'position',[d(1) d(2)-15 d(3)]);
set(get(cbar,'ylabel'),'String', '(kgkg^{-1})','FontSize',fontsize,'FontName','TimesNewRoman')

set(gca,'YAxisLocation','right')
box on   
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

