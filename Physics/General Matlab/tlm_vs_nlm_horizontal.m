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

plot_level = 65;
plot_level_qi = 41;

% 0 to Scale plots individually
% 1 to Scale plots by nonlinear
scale = 1;

plot_u = 0;
plot_v = 0;
plot_t = 0;
plot_q = 0;
plot_p = 0;
plot_qi = 1;
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
file = ['v000_C180.prog.eta.2014020',file_end(1:6),'z.nc4'];

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

%Load Perturbed (analysis) st
cd /discover/nobackup/drholdaw/btmp.25956/prog/prog_replay/
file = ['v000_C180.prog.eta.2014020',file_end(1:6),'z.nc4'];

u_replay = ncread(file,'u',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
v_replay = ncread(file,'v',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
t_replay = ncread(file,'tv',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
q_replay = ncread(file,'sphu',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
p_replay = ncread(file,'delp',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
qi_replay = ncread(file,'qitot',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
ql_replay = ncread(file,'qltot',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
o3_replay = ncread(file,'ozone',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);

%Load first TLM state to compare.
cd /discover/nobackup/drholdaw/wrk.bac/sens.20140202.000000/
file = ['v000_C180.fvpert.eta.2014020',file_end,'z_ST03z_P20000_S1.nc4'];

u_tlmdry = ncread(file,'U',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
v_tlmdry = ncread(file,'V',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
t_tlmdry = ncread(file,'TV',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
q_tlmdry = ncread(file,'QV',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
p_tlmdry = ncread(file,'DP',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
qi_tlmdry = ncread(file,'QI',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
ql_tlmdry = ncread(file,'QL',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
o3_tlmdry = ncread(file,'O3',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);

%Load second TLM state to compare.
cd /discover/nobackup/drholdaw/wrk.bac/sens.20140202.000000/
file = ['v000_C180.fvpert.eta.2014020',file_end,'z_ST03z_P22000_S1.nc4'];

u_tlmmoi = ncread(file,'U',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
v_tlmmoi = ncread(file,'V',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
t_tlmmoi = ncread(file,'TV',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
q_tlmmoi = ncread(file,'QV',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
p_tlmmoi = ncread(file,'DP',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
qi_tlmmoi = ncread(file,'QI',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
ql_tlmmoi = ncread(file,'QL',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
o3_tlmmoi = ncread(file,'O3',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);

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

%Pick out just the level to be plotted.
u_a_lev = u_nlm(:,:,plot_level);
v_a_lev = v_nlm(:,:,plot_level);
t_a_lev = t_nlm(:,:,plot_level);
q_a_lev = q_nlm(:,:,plot_level);
p_a_lev = p_nlm(:,:,plot_level);
qi_a_lev = qi_nlm(:,:,plot_level_qi);
ql_a_lev = ql_nlm(:,:,plot_level);
o3_a_lev = o3_nlm(:,:,plot_level);

%Uncomment here to plot actual TLM pert trajectories
u_b_lev = u_tlmdry(:,:,plot_level);
v_b_lev = v_tlmdry(:,:,plot_level);
t_b_lev = t_tlmdry(:,:,plot_level);
q_b_lev = q_tlmdry(:,:,plot_level);
p_b_lev = p_tlmdry(:,:,plot_level);
qi_b_lev = qi_tlmdry(:,:,plot_level_qi);
ql_b_lev = ql_tlmdry(:,:,plot_level);
o3_b_lev = o3_tlmdry(:,:,plot_level);

u_c_lev = u_tlmmoi(:,:,plot_level);
v_c_lev = v_tlmmoi(:,:,plot_level);
t_c_lev = t_tlmmoi(:,:,plot_level);
q_c_lev = q_tlmmoi(:,:,plot_level);
p_c_lev = p_tlmmoi(:,:,plot_level);
qi_c_lev = qi_tlmmoi(:,:,plot_level_qi);
ql_c_lev = ql_tlmmoi(:,:,plot_level);
o3_c_lev = o3_tlmmoi(:,:,plot_level);

% %Uncomment here to plot the difference between the NL and TLM pert trajectories
% u_b_lev = u_tlmdry(:,:,plot_level) - u_nlm(:,:,plot_level);
% v_b_lev = v_tlmdry(:,:,plot_level) - v_nlm(:,:,plot_level);
% t_b_lev = t_tlmdry(:,:,plot_level) - t_nlm(:,:,plot_level);
% q_b_lev = q_tlmdry(:,:,plot_level) - q_nlm(:,:,plot_level);
% p_b_lev = p_tlmdry(:,:,plot_level) - p_nlm(:,:,plot_level);
% qi_b_lev = qi_tlmdry(:,:,plot_level) - qi_nlm(:,:,plot_level);
% ql_b_lev = ql_tlmdry(:,:,plot_level) - ql_nlm(:,:,plot_level);
% o3_b_lev = o3_tlmdry(:,:,plot_level) - o3_nlm(:,:,plot_level);
% 
% u_c_lev = u_tlmmoi(:,:,plot_level) - u_nlm(:,:,plot_level);
% v_c_lev = v_tlmmoi(:,:,plot_level) - v_nlm(:,:,plot_level);
% t_c_lev = t_tlmmoi(:,:,plot_level) - t_nlm(:,:,plot_level);
% q_c_lev = q_tlmmoi(:,:,plot_level) - q_nlm(:,:,plot_level);
% p_c_lev = p_tlmmoi(:,:,plot_level) - p_nlm(:,:,plot_level);
% qi_c_lev = qi_tlmmoi(:,:,plot_level) - qi_nlm(:,:,plot_level);
% ql_c_lev = ql_tlmmoi(:,:,plot_level) - ql_nlm(:,:,plot_level);
% o3_c_lev = o3_tlmmoi(:,:,plot_level) - o3_nlm(:,:,plot_level);

%Compute maximum absolute values for contour limits and intervals.
maxu  = max(abs( [u_a_lev(:) ; u_b_lev(:) ; u_c_lev(:)]  ));
maxv  = max(abs( [v_a_lev(:) ; v_b_lev(:) ; v_c_lev(:)]  ));
maxt  = max(abs( [t_a_lev(:) ; t_b_lev(:) ; t_c_lev(:)]  ));
maxq  = max(abs( [q_a_lev(:) ; q_b_lev(:) ; q_c_lev(:)]  ));
maxp  = max(abs( [p_a_lev(:) ; p_b_lev(:) ; p_c_lev(:)]  ));
maxqi = max(abs( [qi_a_lev(:); qi_b_lev(:); qi_c_lev(:)] ));
maxql = max(abs( [ql_a_lev(:); ql_b_lev(:); ql_c_lev(:)] ));
maxo3 = max(abs( [o3_a_lev(:); o3_b_lev(:); o3_c_lev(:)] ));

%Set contour intervals for each variable
cint_u = 2*maxu/20;
cint_v = 2*maxv/20;
cint_t = 2*maxt/20;
cint_q = 2*maxq/20;
cint_p = 2*maxp/20;
cint_qi = 2*maxqi/20;
cint_ql = 2*maxql/20;
cint_o3 = 2*maxo3/20;

%Make figures of pert trajectories.
if makeplots == 1 

    if plot_u == 1
        
        figure
        set(gcf,'position',[621 30 658 889])

        subplot(3,1,1)
        [C,h] = contour(lon,lat,((u_a_lev))','LineWidth',line_wid_cont);
        hold on; 
        plot(coast_lon,coast_lat,'Color',[grey grey grey])
        caxis([-maxu maxu])
        title('u')
        colormap(cmap)
        colbar = colorbar;
        set(h,'LevelStep',cint_u);

        subplot(3,1,2)
        [C,h] = contour(lon,lat,((u_b_lev))','LineWidth',line_wid_cont);
        hold on; 
        plot(coast_lon,coast_lat,'Color',[grey grey grey])
        caxis([-maxu maxu])
        title('u')
        colormap(cmap)
        colbar = colorbar;
        set(h,'LevelStep',cint_u);

        subplot(3,1,3)
        [C,h] = contour(lon,lat,((u_c_lev))','LineWidth',line_wid_cont);
        hold on; 
        plot(coast_lon,coast_lat,'Color',[grey grey grey])
        caxis([-maxu maxu])
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
        caxis([-maxv maxv])
        title('v')
        colormap(cmap)
        colbar = colorbar;
        set(h,'LevelStep',cint_v);

        subplot(3,1,2)
        [C,h] = contour(lon,lat,((v_b_lev))','LineWidth',line_wid_cont);
        hold on; 
        plot(coast_lon,coast_lat,'Color',[grey grey grey])
        caxis([-maxv maxv])
        title('v')
        colormap(cmap)
        colbar = colorbar;
        set(h,'LevelStep',cint_v);

        subplot(3,1,3)
        [C,h] = contour(lon,lat,((v_c_lev))','LineWidth',line_wid_cont);
        hold on; 
        plot(coast_lon,coast_lat,'Color',[grey grey grey])
        caxis([-maxv maxv])
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
        caxis([-maxt maxt])
        title('t')
        colormap(cmap)
        colbar = colorbar;
        set(h,'LevelStep',cint_t);

        subplot(3,1,2)
        [C,h] = contour(lon,lat,((t_b_lev))','LineWidth',line_wid_cont);
        hold on; 
        plot(coast_lon,coast_lat,'Color',[grey grey grey])
        caxis([-maxt maxt])
        title('t')
        colormap(cmap)
        colbar = colorbar;
        set(h,'LevelStep',cint_t);

        subplot(3,1,3)
        [C,h] = contour(lon,lat,((t_c_lev))','LineWidth',line_wid_cont);
        hold on; 
        plot(coast_lon,coast_lat,'Color',[grey grey grey])
        caxis([-maxt maxt])
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
%         caxis([-max(abs(q_a_lev(:))) max(abs(q_a_lev(:))) ])
        caxis([-maxq maxq])
        title('q')
        colormap(cmap)
        colbar = colorbar;
        set(h,'LevelStep',cint_q);

        subplot(3,1,2)
        [C,h] = contour(lon,lat,((q_b_lev))','LineWidth',line_wid_cont);
        hold on; 
        plot(coast_lon,coast_lat,'Color',[grey grey grey])
%         caxis([-max(abs(q_b_lev(:))) max(abs(q_b_lev(:))) ])
        caxis([-maxq maxq])
        title('q')
        colormap(cmap)
        colbar = colorbar;
        set(h,'LevelStep',cint_q);

        subplot(3,1,3)
        [C,h] = contour(lon,lat,((q_c_lev))','LineWidth',line_wid_cont);
        hold on; 
        plot(coast_lon,coast_lat,'Color',[grey grey grey])
%         caxis([-max(abs(q_a_lev(:))) max(abs(q_a_lev(:))) ])
        caxis([-maxq maxq])
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
        caxis([-maxp maxp])
        title('p')
        colormap(cmap)
        colbar = colorbar;
        set(h,'LevelStep',cint_p);

        subplot(3,1,2)
        [C,h] = contour(lon,lat,((p_b_lev))','LineWidth',line_wid_cont);
        hold on; 
        plot(coast_lon,coast_lat,'Color',[grey grey grey])
        caxis([-maxp maxp])
        title('p')
        colormap(cmap)
        colbar = colorbar;
        set(h,'LevelStep',cint_p);

        subplot(3,1,3)
        [C,h] = contour(lon,lat,((p_c_lev))','LineWidth',line_wid_cont);
        hold on; 
        plot(coast_lon,coast_lat,'Color',[grey grey grey])
        caxis([-maxp maxp])
        title('p')
        colormap(cmap)
        colbar = colorbar;
        set(h,'LevelStep',cint_p);
    
    end
        
    if plot_qi == 1

        figure
        set(gcf,'position',[621 30 658 889])

        if scale == 0 || scale == 1
            maxqi = max(abs(qi_a_lev(:)));
        end
        cint_qi = 2*maxqi/20;
        
        subplot(3,1,1)
        contour(lon,lat,qi_a_lev','LineWidth',line_wid_cont,'LevelStep',cint_qi);
        hold on;
        plot(coast_lon,coast_lat,'Color',[grey grey grey])
        caxis([-maxqi maxqi])
        title('qi nonlinear perturbation')
        colormap(cmap)
        colorbar;

        if scale == 0
            maxqi = max(abs(qi_b_lev(:)));
        end
        cint_qi = 2*maxqi/20;
        
        subplot(3,1,2)
        contour(lon,lat,qi_b_lev','LineWidth',line_wid_cont,'LevelStep',cint_qi);
        hold on;
        plot(coast_lon,coast_lat,'Color',[grey grey grey])
        caxis([-maxqi maxqi])
        title('qi without linear cloud scheme')
        colormap(cmap)
        colorbar;

        if scale == 0
            maxqi = max(abs(qi_c_lev(:)));
        end
        cint_qi = 2*maxqi/20;
        
        subplot(3,1,3)
        contour(lon,lat,qi_c_lev','LineWidth',line_wid_cont,'LevelStep',cint_qi);
        hold on;
        plot(coast_lon,coast_lat,'Color',[grey grey grey])
        caxis([-maxqi maxqi])
        title('qi with linear cloud scheme')
        colormap(cmap)
        colorbar;
    
    end
        
    if plot_ql == 1

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
        title('ql nonlinear perturbation')
        colormap(cmap)
        colorbar;

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
        title('ql without linear cloud scheme')
        colormap(cmap)
        colorbar;

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
        title('ql with linear cloud scheme')
        colormap(cmap)
        colorbar;
        
    end
        
    if plot_o3 == 1

        figure
        set(gcf,'position',[621 30 658 889])
    
        subplot(3,1,1)
        [C,h] = contour(lon,lat,((o3_a_lev))','LineWidth',line_wid_cont);
        hold on; 
        plot(coast_lon,coast_lat,'Color',[grey grey grey])
%         caxis([-maxo3 maxo3])
        title('o3')
%         colormap(cmap)
        colbar = colorbar;
        set(h,'LevelStep',cint_o3);
        
        subplot(3,1,2)
        [C,h] = contour(lon,lat,((o3_b_lev))','LineWidth',line_wid_cont);
        hold on; 
        plot(coast_lon,coast_lat,'Color',[grey grey grey])
%         caxis([-maxo3 maxo3])
        title('o3')
%         colormap(cmap)
        colbar = colorbar;
        set(h,'LevelStep',cint_o3);
            
        subplot(3,1,3)
        [C,h] = contour(lon,lat,((o3_c_lev))','LineWidth',line_wid_cont);
        hold on; 
        plot(coast_lon,coast_lat,'Color',[grey grey grey])
%         caxis([-maxo3 maxo3])
        title('o3')
%         colormap(cmap)
        colbar = colorbar;
        set(h,'LevelStep',cint_o3);
        
    end
    
    
end


