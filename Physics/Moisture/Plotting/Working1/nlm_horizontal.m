close all
clear
clc

%Load MATLAB topography
load coast
coast_lat = lat; clear lat
coast_lon = long; clear long

%Choose whether you want plots of pert traj (0/1) and which model level to plot 
makeplots = 1;
plot_level = 50;

plot_u = 0;
plot_v = 0;
plot_t = 0;
plot_q = 1;
plot_p = 0;
plot_qi = 0;
plot_ql = 0;
plot_o3 = 0;

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
% file_end = '6_0300';
% file_end = '6_0600';
% file_end = '6_1200';
% file_end = '6_1500';
file_end = '7_0000';

cd /home/drholdaw/Lin_Moist_Physics/Inputs/
pref
 
%Load Free (background) State.
cd /discover/nobackup/drholdaw/tmp.22292/nlm_runs/free/
file = ['x0011dh_a.prog.eta.2013011',file_end(1:4),'z.nc4'];

lon = ncread(file,'lon');
lat = ncread(file,'lat');
lev = ncread(file,'lev');

u_free = ncread(file,'u');
v_free = ncread(file,'v');
t_free = ncread(file,'tv');
q_free = ncread(file,'sphu');
p_free = ncread(file,'delp');
qi_free = ncread(file,'qitot');
ql_free = ncread(file,'qltot');
o3_free = ncread(file,'ozone');

%Load Perturbed (analysis) state.
cd /discover/nobackup/drholdaw/tmp.22292/nlm_runs/replay10/
file = ['x0011dh_a.prog.eta.2013011',file_end(1:4),'z.nc4'];

u_replay10 = ncread(file,'u');
v_replay10 = ncread(file,'v');
t_replay10 = ncread(file,'tv');
q_replay10 = ncread(file,'sphu');
p_replay10 = ncread(file,'delp');
qi_replay10 = ncread(file,'qitot');
ql_replay10 = ncread(file,'qltot');
o3_replay10 = ncread(file,'ozone');

cd /discover/nobackup/drholdaw/tmp.22292/nlm_runs/replay05/
file = ['x0011dh_a.prog.eta.2013011',file_end(1:4),'z.nc4'];

u_replay05 = ncread(file,'u');
v_replay05 = ncread(file,'v');
t_replay05 = ncread(file,'tv');
q_replay05 = ncread(file,'sphu');
p_replay05 = ncread(file,'delp');
qi_replay05 = ncread(file,'qitot');
ql_replay05 = ncread(file,'qltot');
o3_replay05 = ncread(file,'ozone');

cd /discover/nobackup/drholdaw/tmp.22292/nlm_runs/replay01/
file = ['x0011dh_a.prog.eta.2013011',file_end(1:4),'z.nc4'];

u_replay01 = ncread(file,'u');
v_replay01 = ncread(file,'v');
t_replay01 = ncread(file,'tv');
q_replay01 = ncread(file,'sphu');
p_replay01 = ncread(file,'delp');
qi_replay01 = ncread(file,'qitot');
ql_replay01 = ncread(file,'qltot');
o3_replay01 = ncread(file,'ozone');

cd /home/drholdaw/Lin_Moist_Physics/DAS_plotting/Cloud' Fraction'/

%Compute NL perturbation trajectory.
u_nlm10 = u_replay10 - u_free;
v_nlm10 = v_replay10 - v_free;
t_nlm10 = t_replay10 - t_free;
q_nlm10 = q_replay10 - q_free;
p_nlm10 = p_replay10 - p_free;
qi_nlm10 = qi_replay10 - qi_free;
ql_nlm10 = ql_replay10 - ql_free;
o3_nlm10 = o3_replay10 - o3_free;

u_nlm05 = u_replay05 - u_free;
v_nlm05 = v_replay05 - v_free;
t_nlm05 = t_replay05 - t_free;
q_nlm05 = q_replay05 - q_free;
p_nlm05 = p_replay05 - p_free;
qi_nlm05 = qi_replay05 - qi_free;
ql_nlm05 = ql_replay05 - ql_free;
o3_nlm05 = o3_replay05 - o3_free;

u_nlm01 = u_replay01 - u_free;
v_nlm01 = v_replay01 - v_free;
t_nlm01 = t_replay01 - t_free;
q_nlm01 = q_replay01 - q_free;
p_nlm01 = p_replay01 - p_free;
qi_nlm01 = qi_replay01 - qi_free;
ql_nlm01 = ql_replay01 - ql_free;
o3_nlm01 = o3_replay01 - o3_free;


%Pick out just the level to be plotted.
u_a_lev = u_nlm10(:,:,plot_level);
v_a_lev = v_nlm10(:,:,plot_level);
t_a_lev = t_nlm10(:,:,plot_level);
q_a_lev = q_nlm10(:,:,plot_level);
p_a_lev = p_nlm10(:,:,plot_level);
qi_a_lev = qi_nlm10(:,:,plot_level);
ql_a_lev = ql_nlm10(:,:,plot_level);
o3_a_lev = o3_nlm10(:,:,plot_level);

%Uncomment here to plot actual TLM pert trajectories
u_b_lev = u_nlm05(:,:,plot_level)*2;
v_b_lev = v_nlm05(:,:,plot_level)*2;
t_b_lev = t_nlm05(:,:,plot_level)*2;
q_b_lev = q_nlm05(:,:,plot_level)*2;
p_b_lev = p_nlm05(:,:,plot_level)*2;
qi_b_lev = qi_nlm05(:,:,plot_level)*2;
ql_b_lev = ql_nlm05(:,:,plot_level)*2;
o3_b_lev = o3_nlm05(:,:,plot_level)*2;

u_c_lev = u_nlm01(:,:,plot_level)*10;
v_c_lev = v_nlm01(:,:,plot_level)*10;
t_c_lev = t_nlm01(:,:,plot_level)*10;
q_c_lev = q_nlm01(:,:,plot_level)*10;
p_c_lev = p_nlm01(:,:,plot_level)*10;
qi_c_lev = qi_nlm01(:,:,plot_level)*10;
ql_c_lev = ql_nlm01(:,:,plot_level)*10;
o3_c_lev = o3_nlm01(:,:,plot_level)*10;

%Compute maximum absolute values for contour limits and intervals.
maxu = max([max(abs(u_a_lev(:))) max(abs(u_b_lev(:)))  max(abs(u_a_lev(:)))]);
minu = -maxu; 
maxv = max([max(abs(v_a_lev(:))) max(abs(v_b_lev(:)))  max(abs(v_a_lev(:)))]);
minv = -maxv;
maxt = max([max(abs(t_a_lev(:))) max(abs(t_b_lev(:)))  max(abs(t_a_lev(:)))]);
mint = -maxt;
maxq = max([max(abs(q_a_lev(:))) max(abs(q_b_lev(:)))  max(abs(q_a_lev(:)))]);
minq = -maxq;
maxp = max([max(abs(p_a_lev(:))) max(abs(p_b_lev(:)))  max(abs(p_a_lev(:)))]);
minp = -maxp; 
maxqi = max([max(abs(qi_a_lev(:))) max(abs(qi_b_lev(:))) max(abs(qi_a_lev(:)))]);
minqi = -maxqi;
maxql = max([max(abs(ql_a_lev(:))) max(abs(ql_b_lev(:))) max(abs(qi_a_lev(:)))]);
minql = -maxql;
maxo3 = max([max(abs(o3_a_lev(:))) max(abs(o3_b_lev(:))) max(abs(qi_a_lev(:)))]);
mino3 = -maxo3;

%Set contour intervals for each variable
cint_u = 2*maxu/10;
cint_v = 2*maxv/10;
cint_t = 2*maxt/10;
cint_q = 2*maxq/30;
cint_p = 2*maxp/10;
cint_qi = 2*maxqi/40;
cint_ql = 2*maxql/40;
cint_o3 = 2*maxo3/10;

%Make figures of pert trajectories.
if makeplots == 1 

    if plot_u == 1
        
        figure
        set(gcf,'position',[621 30 658 889])

        subplot(3,1,1)
        [C,h] = contour(lon,lat,((u_a_lev))','LineWidth',line_wid_cont);
        hold on; 
        plot(coast_lon,coast_lat,'Color',[grey grey grey])
        caxis([minu maxu])
        title('u')
        colormap(cmap)
        colbar = colorbar;
        set(h,'LevelStep',cint_u);

        subplot(3,1,2)
        [C,h] = contour(lon,lat,((u_b_lev))','LineWidth',line_wid_cont);
        hold on; 
        plot(coast_lon,coast_lat,'Color',[grey grey grey])
        caxis([minu maxu])
        title('u')
        colormap(cmap)
        colbar = colorbar;
        set(h,'LevelStep',cint_u);

        subplot(3,1,3)
        [C,h] = contour(lon,lat,((u_c_lev))','LineWidth',line_wid_cont);
        hold on; 
        plot(coast_lon,coast_lat,'Color',[grey grey grey])
        caxis([minu maxu])
        title('v')
        colormap(cmap)
        colbar = colorbar;
        set(h,'LevelStep',cint_u);
    
    end
        
    if plot_v == 1
    
        figure
        set(gcf,'position',[621 30 658 889])

        subplot(3,1,1)
        [C,h] = contour(lon,lat,((v_a_lev))','LineWidth',line_wid_cont);
        hold on; 
        plot(coast_lon,coast_lat,'Color',[grey grey grey])
        caxis([minv maxv])
        title('v')
        colormap(cmap)
        colbar = colorbar;
        set(h,'LevelStep',cint_v);

        subplot(3,1,2)
        [C,h] = contour(lon,lat,((v_b_lev))','LineWidth',line_wid_cont);
        hold on; 
        plot(coast_lon,coast_lat,'Color',[grey grey grey])
        caxis([minv maxv])
        title('v')
        colormap(cmap)
        colbar = colorbar;
        set(h,'LevelStep',cint_v);

        subplot(3,1,3)
        [C,h] = contour(lon,lat,((v_c_lev))','LineWidth',line_wid_cont);
        hold on; 
        plot(coast_lon,coast_lat,'Color',[grey grey grey])
        caxis([minv maxv])
        title('v')
        colormap(cmap)
        colbar = colorbar;
        set(h,'LevelStep',cint_v);
    
    end
        
    if plot_t == 1

        figure
        set(gcf,'position',[621 30 658 889])

        subplot(3,1,1)
        [C,h] = contour(lon,lat,((t_a_lev))','LineWidth',line_wid_cont);
        hold on; 
        plot(coast_lon,coast_lat,'Color',[grey grey grey])
        caxis([mint maxt])
        title('t')
        colormap(cmap)
        colbar = colorbar;
        set(h,'LevelStep',cint_t);

        subplot(3,1,2)
        [C,h] = contour(lon,lat,((t_b_lev))','LineWidth',line_wid_cont);
        hold on; 
        plot(coast_lon,coast_lat,'Color',[grey grey grey])
        caxis([mint maxt])
        title('t')
        colormap(cmap)
        colbar = colorbar;
        set(h,'LevelStep',cint_t);

        subplot(3,1,3)
        [C,h] = contour(lon,lat,((t_c_lev))','LineWidth',line_wid_cont);
        hold on; 
        plot(coast_lon,coast_lat,'Color',[grey grey grey])
        caxis([mint maxt])
        title('t')
        colormap(cmap)
        colbar = colorbar;
        set(h,'LevelStep',cint_t);
    
    end
        
    if plot_q == 1

        figure
        set(gcf,'position',[621 30 658 889])

        subplot(3,1,1)
        [C,h] = contour(lon,lat,((q_a_lev))','LineWidth',line_wid_cont);
        hold on; 
        plot(coast_lon,coast_lat,'Color',[grey grey grey])
        caxis([minq maxq])
        title('q')
        colormap(cmap)
        colbar = colorbar;
        set(h,'LevelStep',cint_q);

        subplot(3,1,2)
        [C,h] = contour(lon,lat,((q_b_lev))','LineWidth',line_wid_cont);
        hold on; 
        plot(coast_lon,coast_lat,'Color',[grey grey grey])
        caxis([minq maxq])
        title('q')
        colormap(cmap)
        colbar = colorbar;
        set(h,'LevelStep',cint_q);

        subplot(3,1,3)
        [C,h] = contour(lon,lat,((q_c_lev))','LineWidth',line_wid_cont);
        hold on; 
        plot(coast_lon,coast_lat,'Color',[grey grey grey])
        caxis([minq maxq])
        title('q')
        colormap(cmap)
        colbar = colorbar;
        set(h,'LevelStep',cint_q);
    
    end
        
    if plot_p == 1

        figure
        set(gcf,'position',[621 30 658 889])

        subplot(3,1,1)
        [C,h] = contour(lon,lat,((p_a_lev))','LineWidth',line_wid_cont);
        hold on; 
        plot(coast_lon,coast_lat,'Color',[grey grey grey])
        caxis([minp maxp])
        title('p')
        colormap(cmap)
        colbar = colorbar;
        set(h,'LevelStep',cint_p);

        subplot(3,1,2)
        [C,h] = contour(lon,lat,((p_b_lev))','LineWidth',line_wid_cont);
        hold on; 
        plot(coast_lon,coast_lat,'Color',[grey grey grey])
        caxis([minp maxp])
        title('p')
        colormap(cmap)
        colbar = colorbar;
        set(h,'LevelStep',cint_p);

        subplot(3,1,3)
        [C,h] = contour(lon,lat,((p_c_lev))','LineWidth',line_wid_cont);
        hold on; 
        plot(coast_lon,coast_lat,'Color',[grey grey grey])
        caxis([minp maxp])
        title('p')
        colormap(cmap)
        colbar = colorbar;
        set(h,'LevelStep',cint_p);
    
    end
        
    if plot_qi == 1

        figure
        set(gcf,'position',[621 30 658 889])

        subplot(3,1,1)
        [C,h] = contour(lon,lat,((qi_a_lev))','LineWidth',line_wid_cont);
        hold on; 
        plot(coast_lon,coast_lat,'Color',[grey grey grey])
        caxis([minqi maxqi])
        title('qi')
        colormap(cmap)
        colbar = colorbar;
        set(h,'LevelStep',cint_qi);

        subplot(3,1,2)
        [C,h] = contour(lon,lat,((qi_b_lev))','LineWidth',line_wid_cont);
        hold on; 
        plot(coast_lon,coast_lat,'Color',[grey grey grey])
        caxis([minqi maxqi])
        title('qi')
        colormap(cmap)
        colbar = colorbar;
        set(h,'LevelStep',cint_qi);

        subplot(3,1,3)
        [C,h] = contour(lon,lat,((qi_c_lev))','LineWidth',line_wid_cont);
        hold on; 
        plot(coast_lon,coast_lat,'Color',[grey grey grey])
        caxis([minqi maxqi])
        title('qi')
        colormap(cmap)
        colbar = colorbar;
        set(h,'LevelStep',cint_qi);
    
    end
        
    if plot_ql == 1

        figure
        set(gcf,'position',[621 30 658 889])

        subplot(3,1,1)
        [C,h] = contour(lon,lat,((ql_a_lev))','LineWidth',line_wid_cont);
        hold on; 
        plot(coast_lon,coast_lat,'Color',[grey grey grey])
        caxis([minql maxql])    
        title('ql')
        colormap(cmap)
        colbar = colorbar;
        set(h,'LevelStep',cint_ql);

        subplot(3,1,2)
        [C,h] = contour(lon,lat,((ql_b_lev))','LineWidth',line_wid_cont);
        hold on; 
        plot(coast_lon,coast_lat,'Color',[grey grey grey])
        caxis([minql maxql])    
        title('ql')
        colormap(cmap)
        colbar = colorbar;
        set(h,'LevelStep',cint_ql);

        subplot(3,1,3)
        [C,h] = contour(lon,lat,((ql_c_lev))','LineWidth',line_wid_cont);
        hold on; 
        plot(coast_lon,coast_lat,'Color',[grey grey grey])
        caxis([minql maxql])    
        title('ql')
        colormap(cmap)
        colbar = colorbar;
        set(h,'LevelStep',cint_ql);
    
    end
        
    if plot_o3 == 1

        figure
        set(gcf,'position',[621 30 658 889])
    
        subplot(3,1,1)
        [C,h] = contour(lon,lat,((o3_a_lev))','LineWidth',line_wid_cont);
        hold on; 
        plot(coast_lon,coast_lat,'Color',[grey grey grey])
        caxis([mino3 maxo3])
        title('o3')
        colormap(cmap)
        colbar = colorbar;
        set(h,'LevelStep',cint_o3);
        
        subplot(3,1,2)
        [C,h] = contour(lon,lat,((o3_b_lev))','LineWidth',line_wid_cont);
        hold on; 
        plot(coast_lon,coast_lat,'Color',[grey grey grey])
        caxis([mino3 maxo3])
        title('o3')
        colormap(cmap)
        colbar = colorbar;
        set(h,'LevelStep',cint_o3);
            
        subplot(3,1,3)
        [C,h] = contour(lon,lat,((o3_c_lev))','LineWidth',line_wid_cont);
        hold on; 
        plot(coast_lon,coast_lat,'Color',[grey grey grey])
        caxis([mino3 maxo3])
        title('o3')
        colormap(cmap)
        colbar = colorbar;
        set(h,'LevelStep',cint_o3);
        
    end
    
    
end


