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
plot_level = 60;
match_cax = 0;

plot_u = 0;
plot_v = 0;
plot_t = 1;
plot_q = 0;
plot_p = 0;
plot_qi = 0;
plot_ql = 0;
plot_o3 = 0;
plot_d1 = 0;
plot_d2 = 0;
plot_d3 = 0;
plot_d4 = 0;
plot_d5 = 0;

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
% file_end = '1_0100';
% file_end = '1_0300';
% file_end = '1_0600';
% file_end = '1_0900';
% file_end = '1_1200';
% file_end = '1_1500';
% file_end = '1_1800';
% file_end = '1_2100';
% file_end = '2_0000';
file_end = '2_0300';

im = 576;
jm = 361;
lm = 72;

LevMin = 1;
LevMax = lm;

cd /home/drholdaw/LinearisedPhysics/Inputs/
pref
 
%Load Free (background) State.
cd /discover/nobackup/drholdaw/btmp.25956/prog/prog_free/
file = ['v000_C180.prog.eta.2014020',file_end,'z.nc4'];

lon = ncread(file,'lon');
lat = ncread(file,'lat');
lev = ncread(file,'lev');

u_free = ncread(file,'u',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
v_free = ncread(file,'v',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
t_free = ncread(file,'tv',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
q_free = ncread(file,'sphu',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
p_free = ncread(file,'delp',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
qi_free = ncread(file,'qitot',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
ql_free = ncread(file,'qltot',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
o3_free = ncread(file,'ozone',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
d1_free = ncread(file,'DU001',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
d2_free = ncread(file,'DU002',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
d3_free = ncread(file,'DU003',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
d4_free = ncread(file,'DU004',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
d5_free = ncread(file,'DU005',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);

%Load Perturbed (analysis) st
cd /discover/nobackup/drholdaw/btmp.25956/prog/prog_replay/
file = ['v000_C180.prog.eta.2014020',file_end,'z.nc4'];

u_replay = ncread(file,'u',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
v_replay = ncread(file,'v',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
t_replay = ncread(file,'tv',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
q_replay = ncread(file,'sphu',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
p_replay = ncread(file,'delp',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
qi_replay = ncread(file,'qitot',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
ql_replay = ncread(file,'qltot',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
o3_replay = ncread(file,'ozone',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
d1_replay = ncread(file,'DU001',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
d2_replay = ncread(file,'DU002',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
d3_replay = ncread(file,'DU003',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
d4_replay = ncread(file,'DU004',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
d5_replay = ncread(file,'DU005',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);

%Load first TLM state to compare.
cd /discover/nobackup/drholdaw/btmp.25956/sens.20140202.000000/
file = ['v000_C180.fvpert.eta.2014020',file_end,'z_ST03z_DUST_P21000_S0.nc4'];
% file = 'fvpertX.eta.nc4';

u_tlm1 = ncread(file,'U',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
v_tlm1 = ncread(file,'V',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
t_tlm1 = ncread(file,'TV',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
q_tlm1 = ncread(file,'QV',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
p_tlm1 = ncread(file,'DP',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
qi_tlm1 = ncread(file,'QI',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
ql_tlm1 = ncread(file,'QL',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
o3_tlm1 = ncread(file,'O3',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
d1_tlm1 = ncread(file,'DU001',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
d2_tlm1 = ncread(file,'DU002',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
d3_tlm1 = ncread(file,'DU003',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
d4_tlm1 = ncread(file,'DU004',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
d5_tlm1 = ncread(file,'DU005',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);

% u_tlm1 = ncread(file,'u',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
% v_tlm1 = ncread(file,'v',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
% t_tlm1 = ncread(file,'tv',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
% q_tlm1 = ncread(file,'sphu',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
% p_tlm1 = ncread(file,'delp',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
% qi_tlm1 = ncread(file,'qitot',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
% ql_tlm1 = ncread(file,'qltot',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
% o3_tlm1 = ncread(file,'ozone',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
% d1_tlm1 = ncread(file,'DU001',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
% d2_tlm1 = ncread(file,'DU002',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
% d3_tlm1 = ncread(file,'DU003',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
% d4_tlm1 = ncread(file,'DU004',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
% d5_tlm1 = ncread(file,'DU005',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);

%Load second TLM state to compare.
cd /discover/nobackup/drholdaw/btmp.25956/sens.20140202.000000/
file = ['v000_C180.fvpert.eta.2014020',file_end,'z_ST03z_DUST_P22010_S0.nc4'];
% file = 'fvpertX.eta.nc4';

u_tlm2 = ncread(file,'U',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
v_tlm2 = ncread(file,'V',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
t_tlm2 = ncread(file,'TV',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
q_tlm2 = ncread(file,'QV',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
p_tlm2 = ncread(file,'DP',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
qi_tlm2 = ncread(file,'QI',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
ql_tlm2 = ncread(file,'QL',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
o3_tlm2 = ncread(file,'O3',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
d1_tlm2 = ncread(file,'DU001',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
d2_tlm2 = ncread(file,'DU002',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
d3_tlm2 = ncread(file,'DU003',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
d4_tlm2 = ncread(file,'DU004',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
d5_tlm2 = ncread(file,'DU005',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);

% u_tlm2 = ncread(file,'u',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
% v_tlm2 = ncread(file,'v',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
% t_tlm2 = ncread(file,'tv',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
% q_tlm2 = ncread(file,'sphu',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
% p_tlm2 = ncread(file,'delp',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
% qi_tlm2 = ncread(file,'qitot',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
% ql_tlm2 = ncread(file,'qltot',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
% o3_tlm2 = ncread(file,'ozone',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
% d1_tlm2 = ncread(file,'DU001',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
% d2_tlm2 = ncread(file,'DU002',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
% d3_tlm2 = ncread(file,'DU003',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
% d4_tlm2 = ncread(file,'DU004',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
% d5_tlm2 = ncread(file,'DU005',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);

cd(mydir)

%Compute NL perturbation trajectory.
u_nlm = u_replay - u_free;
v_nlm = v_replay - v_free;
t_nlm = t_replay - t_free;
q_nlm = q_replay - q_free;
p_nlm = p_replay - p_free;
qi_nlm = qi_replay - qi_free;
ql_nlm = ql_replay - ql_free;
o3_nlm = o3_replay - o3_free;
d1_nlm = d1_replay - d1_free;
d2_nlm = d2_replay - d2_free;
d3_nlm = d3_replay - d3_free;
d4_nlm = d4_replay - d4_free;
d5_nlm = d5_replay - d5_free;

%Pick out just the level to be plotted.
u_a_lev = u_nlm(:,:,plot_level);
v_a_lev = v_nlm(:,:,plot_level);
t_a_lev = t_nlm(:,:,plot_level);
q_a_lev = q_nlm(:,:,plot_level);
p_a_lev = p_nlm(:,:,plot_level);
qi_a_lev = qi_nlm(:,:,plot_level);
ql_a_lev = ql_nlm(:,:,plot_level);
o3_a_lev = o3_nlm(:,:,plot_level);
d1_a_lev = d1_nlm(:,:,plot_level);
d2_a_lev = d2_nlm(:,:,plot_level);
d3_a_lev = d3_nlm(:,:,plot_level);
d4_a_lev = d4_nlm(:,:,plot_level);
d5_a_lev = d5_nlm(:,:,plot_level);

%Uncomment here to plot actual TLM pert trajectories
u_b_lev = u_tlm1(:,:,plot_level);
v_b_lev = v_tlm1(:,:,plot_level);
t_b_lev = t_tlm1(:,:,plot_level);
q_b_lev = q_tlm1(:,:,plot_level);
p_b_lev = p_tlm1(:,:,plot_level);
qi_b_lev = qi_tlm1(:,:,plot_level);
ql_b_lev = ql_tlm1(:,:,plot_level);
o3_b_lev = o3_tlm1(:,:,plot_level);
d1_b_lev = d1_tlm1(:,:,plot_level);
d2_b_lev = d2_tlm1(:,:,plot_level);
d3_b_lev = d3_tlm1(:,:,plot_level);
d4_b_lev = d4_tlm1(:,:,plot_level);
d5_b_lev = d5_tlm1(:,:,plot_level);

u_c_lev = u_tlm2(:,:,plot_level);
v_c_lev = v_tlm2(:,:,plot_level);
t_c_lev = t_tlm2(:,:,plot_level);
q_c_lev = q_tlm2(:,:,plot_level);
p_c_lev = p_tlm2(:,:,plot_level);
qi_c_lev = qi_tlm2(:,:,plot_level);
ql_c_lev = ql_tlm2(:,:,plot_level);
o3_c_lev = o3_tlm2(:,:,plot_level);
d1_c_lev = d1_tlm2(:,:,plot_level);
d2_c_lev = d2_tlm2(:,:,plot_level);
d3_c_lev = d3_tlm2(:,:,plot_level);
d4_c_lev = d4_tlm2(:,:,plot_level);
d5_c_lev = d5_tlm2(:,:,plot_level);

% %Uncomment here to plot the difference between the NL and TLM pert trajectories
% u_b_lev = u_tlm1(:,:,plot_level) - u_nlm(:,:,plot_level);
% v_b_lev = v_tlm1(:,:,plot_level) - v_nlm(:,:,plot_level);
% t_b_lev = t_tlm1(:,:,plot_level) - t_nlm(:,:,plot_level);
% q_b_lev = q_tlm1(:,:,plot_level) - q_nlm(:,:,plot_level);
% p_b_lev = p_tlm1(:,:,plot_level) - p_nlm(:,:,plot_level);
% qi_b_lev = qi_tlm1(:,:,plot_level) - qi_nlm(:,:,plot_level);
% ql_b_lev = ql_tlm1(:,:,plot_level) - ql_nlm(:,:,plot_level);
% o3_b_lev = o3_tlm1(:,:,plot_level) - o3_nlm(:,:,plot_level);
% 
% u_c_lev = u_tlm2(:,:,plot_level) - u_nlm(:,:,plot_level);
% v_c_lev = v_tlm2(:,:,plot_level) - v_nlm(:,:,plot_level);
% t_c_lev = t_tlm2(:,:,plot_level) - t_nlm(:,:,plot_level);
% q_c_lev = q_tlm2(:,:,plot_level) - q_nlm(:,:,plot_level);
% p_c_lev = p_tlm2(:,:,plot_level) - p_nlm(:,:,plot_level);
% qi_c_lev = qi_tlm2(:,:,plot_level) - qi_nlm(:,:,plot_level);
% ql_c_lev = ql_tlm2(:,:,plot_level) - ql_nlm(:,:,plot_level);
% o3_c_lev = o3_tlm2(:,:,plot_level) - o3_nlm(:,:,plot_level);

%Compute maximum absolute values for contour limits and intervals.
maxu = max([max(abs(u_a_lev(:))) max(abs(u_b_lev(:)))  max(abs(u_c_lev(:)))]);
minu = -maxu; 
maxv = max([max(abs(v_a_lev(:))) max(abs(v_b_lev(:)))  max(abs(v_c_lev(:)))]);
minv = -maxv;
maxt = max([max(abs(t_a_lev(:))) max(abs(t_b_lev(:)))  max(abs(t_c_lev(:)))]);
mint = -maxt;
maxq = max([max(abs(q_a_lev(:))) max(abs(q_b_lev(:))) ]);
minq = -maxq;
maxp = max([max(abs(p_a_lev(:))) max(abs(p_b_lev(:)))  max(abs(p_c_lev(:)))]);
minp = -maxp; 
maxqi = max([max(abs(qi_a_lev(:))) max(abs(qi_b_lev(:))) ]);
minqi = -maxqi;
maxql = max([max(abs(ql_a_lev(:))) max(abs(qi_b_lev(:))) ]);
minql = -maxql;
maxo3 = max([max(abs(o3_a_lev(:))) max(abs(o3_b_lev(:))) max(abs(o3_c_lev(:)))]);
mino3 = -maxo3;
maxd1 = max([max(abs(d1_a_lev(:))) max(abs(d1_b_lev(:))) max(abs(d1_c_lev(:)))]);
mind1 = -maxd1;
maxd2 = max([max(abs(d2_a_lev(:))) max(abs(d2_b_lev(:))) max(abs(d2_c_lev(:)))]);
mind2 = -maxd2;
maxd3 = max([max(abs(d3_a_lev(:))) max(abs(d3_b_lev(:))) max(abs(d3_c_lev(:)))]);
mind3 = -maxd3;
maxd4 = max([max(abs(d4_a_lev(:))) max(abs(d4_b_lev(:))) max(abs(d4_c_lev(:)))]);
mind4 = -maxd4;
maxd5 = max([max(abs(d5_a_lev(:))) max(abs(d5_b_lev(:))) max(abs(d5_c_lev(:)))]);
mind5 = -maxd5;

%Set contour intervals for each variable
cint_u = 2*maxu/10;
cint_v = 2*maxv/10;
cint_t = 2*maxt/10;
cint_q = 2*maxq/30;
cint_p = 2*maxp/10;
cint_qi = 2*maxqi/20;
cint_ql = 2*maxql/10;
cint_o3 = 2*maxo3/10;
cint_d1 = 2*maxd1/10;
cint_d2 = 2*maxd2/10;
cint_d3 = 2*maxd3/10;
cint_d4 = 2*maxd4/10;
cint_d5 = 2*maxd5/10;

%Make figures of pert trajectories.
if makeplots == 1 

    if plot_u == 1
        
        figure
        set(gcf,'position',[621 30 658 889])

        subplot(3,1,1)
        contour(lon,lat,((u_a_lev))','LineWidth',line_wid_cont);
        hold on; 
        plot(coast_lon,coast_lat,'Color',[grey grey grey])
        caxis([minu maxu])
        title('u')
        colormap(cmap)
        colorbar;
        set(h,'LevelStep',cint_u);

        subplot(3,1,2)
        contour(lon,lat,((u_b_lev))','LineWidth',line_wid_cont);
        hold on; 
        plot(coast_lon,coast_lat,'Color',[grey grey grey])
        caxis([minu maxu])
        title('u')
        colormap(cmap)
        colorbar;
        set(h,'LevelStep',cint_u);

        subplot(3,1,3)
        contour(lon,lat,((u_c_lev))','LineWidth',line_wid_cont);
        hold on; 
        plot(coast_lon,coast_lat,'Color',[grey grey grey])
        caxis([minu maxu])
        title('v')
        colormap(cmap)
        colorbar;
        set(h,'LevelStep',cint_u);
    
    end
        
    if plot_v == 1
    
        figure
        set(gcf,'position',[621 30 658 889])

        subplot(3,1,1)
        contour(lon,lat,((v_a_lev))','LineWidth',line_wid_cont);
        hold on; 
        plot(coast_lon,coast_lat,'Color',[grey grey grey])
        caxis([minv maxv])
        title('v')
        colormap(cmap)
        colorbar;
        set(h,'LevelStep',cint_v);

        subplot(3,1,2)
        contour(lon,lat,((v_b_lev))','LineWidth',line_wid_cont);
        hold on; 
        plot(coast_lon,coast_lat,'Color',[grey grey grey])
        caxis([minv maxv])
        title('v')
        colormap(cmap)
        colorbar;
        set(h,'LevelStep',cint_v);

        subplot(3,1,3)
        contour(lon,lat,((v_c_lev))','LineWidth',line_wid_cont);
        hold on; 
        plot(coast_lon,coast_lat,'Color',[grey grey grey])
        caxis([minv maxv])
        title('v')
        colormap(cmap)
        colorbar;
        set(h,'LevelStep',cint_v);
    
    end
        
    if plot_t == 1

        figure
        set(gcf,'position',[621 30 658 889])

        subplot(3,1,1)
        [~,h] = contour(lon,lat,((t_a_lev))','LineWidth',line_wid_cont);
        hold on; 
        plot(coast_lon,coast_lat,'Color',[grey grey grey])
        if match_cax == 1
            caxis([mint maxt])
        else
            caxis([-max(abs(t_a_lev(:))) max(abs(t_a_lev(:)))])
        end
        title('t')
        colormap(cmap)
        colorbar;
        if match_cax == 1
            set(h,'LevelStep',cint_t);
        end
            
        subplot(3,1,2)
        [~,h] = contour(lon,lat,((t_b_lev))','LineWidth',line_wid_cont);
        hold on; 
        plot(coast_lon,coast_lat,'Color',[grey grey grey])
        if match_cax == 1
            caxis([mint maxt])
        else
            caxis([-max(abs(t_b_lev(:))) max(abs(t_b_lev(:)))])
        end
        title('t')
        colormap(cmap)
        colorbar;
        if match_cax == 1
            set(h,'LevelStep',cint_t);
        end

        subplot(3,1,3)
        [~,h] = contour(lon,lat,(min(abs(t_c_lev),max(abs(t_a_lev(:)))))'.*(abs(t_c_lev)./t_c_lev)','LineWidth',line_wid_cont);
        hold on; 
        plot(coast_lon,coast_lat,'Color',[grey grey grey])
        if match_cax == 1
            caxis([mint maxt])
        else
%             caxis([-max(abs(t_c_lev(:))) max(abs(t_c_lev(:)))])
        end
        title('t')
        colormap(cmap)
        colorbar;
        if match_cax == 1
            set(h,'LevelStep',cint_t);
        end
    
    end
        
    if plot_q == 1

        figure
        set(gcf,'position',[621 30 658 889])

        subplot(3,1,1)
        contour(lon,lat,((q_a_lev))','LineWidth',line_wid_cont);
        hold on; 
        plot(coast_lon,coast_lat,'Color',[grey grey grey])
%         caxis([-max(abs(q_a_lev(:))) max(abs(q_a_lev(:))) ])
        caxis([minq maxq])
        title('q')
        colormap(cmap)
        colorbar;
        set(h,'LevelStep',cint_q);

        subplot(3,1,2)
        contour(lon,lat,((q_b_lev))','LineWidth',line_wid_cont);
        hold on; 
        plot(coast_lon,coast_lat,'Color',[grey grey grey])
%         caxis([-max(abs(q_b_lev(:))) max(abs(q_b_lev(:))) ])
        caxis([minq maxq])
        title('q')
        colormap(cmap)
        colorbar;
        set(h,'LevelStep',cint_q);

        subplot(3,1,3)
        contour(lon,lat,((q_c_lev))','LineWidth',line_wid_cont);
        hold on; 
        plot(coast_lon,coast_lat,'Color',[grey grey grey])
%         caxis([-max(abs(q_a_lev(:))) max(abs(q_a_lev(:))) ])
        caxis([minq maxq])
        title('q')
        colormap(cmap)
        colorbar;
        set(h,'LevelStep',cint_q);
    
    end
        
    if plot_p == 1

        figure
        set(gcf,'position',[621 30 658 889])

        subplot(3,1,1)
        contour(lon,lat,((p_a_lev))','LineWidth',line_wid_cont);
        hold on; 
        plot(coast_lon,coast_lat,'Color',[grey grey grey])
        caxis([minp maxp])
        title('p')
        colormap(cmap)
        colorbar;
        set(h,'LevelStep',cint_p);

        subplot(3,1,2)
        contour(lon,lat,((p_b_lev))','LineWidth',line_wid_cont);
        hold on; 
        plot(coast_lon,coast_lat,'Color',[grey grey grey])
        caxis([minp maxp])
        title('p')
        colormap(cmap)
        colorbar;
        set(h,'LevelStep',cint_p);

        subplot(3,1,3)
        contour(lon,lat,((p_c_lev))','LineWidth',line_wid_cont);
        hold on; 
        plot(coast_lon,coast_lat,'Color',[grey grey grey])
        caxis([minp maxp])
        title('p')
        colormap(cmap)
        colorbar;
        set(h,'LevelStep',cint_p);
    
    end
        
    if plot_qi == 1

        figure
        set(gcf,'position',[621 30 658 889])

        subplot(3,1,1)
        contour(lon,lat,((qi_a_lev))','LineWidth',line_wid_cont);
        hold on; 
        plot(coast_lon,coast_lat,'Color',[grey grey grey])
%         caxis([minqi maxqi])
        title('qi')
        colormap(cmap)
        colorbar;
%         set(h,'LevelStep',cint_qi);

        subplot(3,1,2)
        contour(lon,lat,((qi_b_lev))','LineWidth',line_wid_cont);
        hold on; 
        plot(coast_lon,coast_lat,'Color',[grey grey grey])
        caxis([-max(abs(qi_b_lev(:))) max(abs(qi_b_lev(:))) ])
%         caxis([minqi maxqi])
        title('qi')
        colormap(cmap)
        colorbar;
%         set(h,'LevelStep',cint_qi);

        subplot(3,1,3)
        contour(lon,lat,((qi_c_lev))','LineWidth',line_wid_cont);
        hold on; 
        plot(coast_lon,coast_lat,'Color',[grey grey grey])
%         caxis([minqi maxqi])
        caxis([-max(abs(qi_c_lev(:))) max(abs(qi_c_lev(:))) ])
        title('qi')
        colormap(cmap)
        colorbar;
%         set(h,'LevelStep',cint_qi);
    
    end
        
    if plot_ql == 1

        figure
        set(gcf,'position',[621 30 658 889])

        subplot(3,1,1)
        contour(lon,lat,((ql_a_lev))','LineWidth',line_wid_cont);
        hold on; 
        plot(coast_lon,coast_lat,'Color',[grey grey grey])
caxis([-max(abs(ql_a_lev(:))) max(abs(ql_a_lev(:)))])
%         title('ql')
%         colormap(cmap)
        colorbar;
        set(h,'LevelStep',cint_ql);

        subplot(3,1,2)
        contour(lon,lat,((ql_b_lev))','LineWidth',line_wid_cont);
        hold on; 
        plot(coast_lon,coast_lat,'Color',[grey grey grey])
        caxis([-max(abs(ql_a_lev(:))) max(abs(ql_a_lev(:)))])  
%         title('ql')
%         colormap(cmap)
        colorbar;
        set(h,'LevelStep',cint_ql);

        subplot(3,1,3)
        contour(lon,lat,((ql_c_lev))','LineWidth',line_wid_cont);
        hold on; 
        plot(coast_lon,coast_lat,'Color',[grey grey grey])
        caxis([-max(abs(ql_a_lev(:))) max(abs(ql_a_lev(:)))])  

        %         title('ql')
%         colormap(cmap)
        colorbar;
        set(h,'LevelStep',cint_ql);
    
    end
        
    if plot_o3 == 1

        figure
        set(gcf,'position',[621 30 658 889])
    
        subplot(3,1,1)
        contour(lon,lat,((o3_a_lev))','LineWidth',line_wid_cont);
        hold on; 
        plot(coast_lon,coast_lat,'Color',[grey grey grey])
%         caxis([mino3 maxo3])
        title('o3')
%         colormap(cmap)
        colorbar;
        set(h,'LevelStep',cint_o3);
        
        subplot(3,1,2)
        contour(lon,lat,((o3_b_lev))','LineWidth',line_wid_cont);
        hold on; 
        plot(coast_lon,coast_lat,'Color',[grey grey grey])
%         caxis([mino3 maxo3])
        title('o3')
%         colormap(cmap)
        colorbar;
        set(h,'LevelStep',cint_o3);
            
        subplot(3,1,3)
        contour(lon,lat,((o3_c_lev))','LineWidth',line_wid_cont);
        hold on; 
        plot(coast_lon,coast_lat,'Color',[grey grey grey])
%         caxis([mino3 maxo3])
        title('o3')
%         colormap(cmap)
        colorbar;
        set(h,'LevelStep',cint_o3);
        
    end
    
    if plot_d1 == 1

        figure
        set(gcf,'position',[621 30 658 889])
    
        subplot(3,1,1)
        contour(lon,lat,((d1_a_lev))','LineWidth',line_wid_cont);
        hold on; 
        plot(coast_lon,coast_lat,'Color',[grey grey grey])
        caxis([-max(abs(d1_a_lev(:))) max(abs(d1_a_lev(:)))])
        caxis([mind1 maxd1])
        title('d1')
        colormap(cmap)
        colorbar;
%         set(h,'LevelStep',cint_d1);
        
        subplot(3,1,2)
        contour(lon,lat,((d1_b_lev))','LineWidth',line_wid_cont);
        hold on; 
        plot(coast_lon,coast_lat,'Color',[grey grey grey])
        caxis([-max(abs(d1_b_lev(:))) max(abs(d1_b_lev(:)))])
%         caxis([mind1 maxd1])
        title('d1')
        colormap(cmap)
        colorbar;
%         set(h,'LevelStep',cint_d1);
            
        subplot(3,1,3)
        contour(lon,lat,((d1_c_lev))','LineWidth',line_wid_cont);
        hold on; 
        plot(coast_lon,coast_lat,'Color',[grey grey grey])
        caxis([-max(abs(d1_c_lev(:))) max(abs(d1_c_lev(:)))])
%         caxis([mind1 maxd1])
        title('d1')
        colormap(cmap)
        colorbar;
%         set(h,'LevelStep',cint_d1);
        
    end
    
end


