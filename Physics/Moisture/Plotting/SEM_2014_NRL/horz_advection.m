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

plot_level = 56;
plot_level_qi = 40;

% 0 to Scale plots individually
% 1 to Scale plots by nonlinear
scale = 1;

%Make center of colormap white
cmap = colormap;
cmap(32,:) = [1 1 1];
cmap(33,:) = [1 1 1];
close

%Some plotting things
fontsize = 11;
line_wid_cont = 1.0;
line_wid_det = 0.6;
grey = 0.75;

%Choose the time being considered: 16th at 12z would be 6_12
file_end = '2_0000';

im = 576;
jm = 361;
lm = 72;

LevMin = 1;
LevMax = lm;

cd /home/drholdaw/LinearisedPhysics/Inputs/
pref
 
%Load Free (background) State.
cd /discover/nobackup/drholdaw/wrk.advec/Model_Output/hord_tr12/prog_free/
file = ['v000_C180.prog.eta.2014020',file_end(1:6),'z.nc4'];

lon = ncread(file,'lon');
lat = ncread(file,'lat');
lev = ncread(file,'lev');

ql_free = ncread(file,'qltot',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);

%Load Perturbed (analysis) st
cd /discover/nobackup/drholdaw/wrk.advec/Model_Output/hord_tr12/prog_replay_b_1p0001/
file = ['v000_C180.prog.eta.2014020',file_end(1:6),'z.nc4'];

ql_replay = ncread(file,'qltot',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);

%Load first TLM state to compare.
cd /discover/nobackup/drholdaw/wrk.advec/Model_Output/hord_tr12/fvpert/
file = ['v000_C180.fvpert.eta.2014020',file_end,'z_P20000.nc4'];

ql_tlm = ncread(file,'QL',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);

%Load first TLM state to compare.
cd /discover/nobackup/drholdaw/wrk.advec/Model_Output/hord_tr1/fvpert_tr12traj/
file = ['v000_C180.fvpert.eta.2014020',file_end,'z_P20000.nc4'];

ql_tlm1 = ncread(file,'QL',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);

cd(mydir)

%Compute NL perturbation trajectory.
ql_nlm = ql_replay - ql_free;

%Pick out just the level to be plotted.
ql_a_lev = ql_nlm(:,:,plot_level);
ql_b_lev = ql_tlm(:,:,plot_level);
ql_c_lev = ql_tlm1(:,:,plot_level);



%Compute maximum absolute values for contour limits and intervals.
maxql = max(abs( [ql_a_lev(:); ql_b_lev(:); ql_c_lev(:)] ));


%Set contour intervals for each variable
cint_ql = 2*maxql/10;

fontsize = 13;

figure
set(gcf,'position',[621 30 658 889])

if scale == 0 || scale == 1
    maxql = max(abs(ql_a_lev(:)));
end
cint_ql = 2*maxql/20;

subplot(3,1,1)
contour(lon,lat,ql_a_lev','LineWidth',line_wid_cont,'LevelStep',cint_ql);
hold on;
plot(coast_lon,coast_lat,'Color',[grey grey grey])
caxis([-maxql maxql])
title('NLM ql (PPM Lin Scheme)','FontSize',fontsize,'FontName','TimesNewRoman')
colormap(cmap)
colorbar;
box on
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

if scale == 0
    maxql = max(abs(ql_b_lev(:)));
end
cint_ql = 2*maxql/20;

if scale == 1
    for i = 1:im
        for j = 1:jm
            if abs(ql_b_lev(i,j)) >= maxql
                ql_b_lev(i,j) = maxql * abs(ql_b_lev(i,j))/ql_b_lev(i,j);
            end
        end
    end         
end

subplot(3,1,2)
contour(lon,lat,ql_b_lev','LineWidth',line_wid_cont,'LevelStep',cint_ql);
hold on;
plot(coast_lon,coast_lat,'Color',[grey grey grey])
caxis([-maxql maxql])
title('TLM ql (PPM Lin Scheme)','FontSize',fontsize,'FontName','TimesNewRoman')
colormap(cmap)
colorbar;
box on
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

        
if scale == 0
    maxql = max(abs(ql_c_lev(:)));
end
cint_ql = 2*maxql/20;

if scale == 1
    for i = 1:im
        for j = 1:jm
            if abs(ql_c_lev(i,j)) >= maxql
                ql_c_lev(i,j) = maxql * abs(ql_c_lev(i,j))/ql_c_lev(i,j);
            end
        end
    end         
end

subplot(3,1,3)
contour(lon,lat,ql_c_lev','LineWidth',line_wid_cont,'LevelStep',cint_ql);
hold on;
plot(coast_lon,coast_lat,'Color',[grey grey grey])
caxis([-maxql maxql])
title('TLM ql (First Order Scheme)','FontSize',fontsize,'FontName','TimesNewRoman')
colormap(cmap)
colorbar;
box on
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')



