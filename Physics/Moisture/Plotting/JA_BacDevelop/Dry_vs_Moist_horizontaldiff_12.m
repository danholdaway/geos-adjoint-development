close all
clear
clc

load coast
coast_lat = lat; clear lat
coast_lon = long; clear long

%Choose whether you want plots of pert traj (0/1) and which model level to plot 
makeplots = 1;
plot_level = 67;
plot_level_q  = 63;
plot_level_ql = 70;
plot_level_qi = 44;

plot_u = 0;
plot_v = 0;
plot_t = 0;
plot_q = 0;
plot_qi = 1;
plot_ql = 1;

%Choose Region,
% 1 = Global
% 2 = Tropics 30S to 30N, below 100hPa
% 3 = Northern hemisphere
% 4 = Southern hemisphere
region = 1;

%Choose the time being considered: 16th at 12z would be 6_12
% file_end = '6_0600'; % 3 hours
% file_end = '6_0900'; % 6 hours
% file_end = '6_1200'; % 9 hours
file_end = '6_1500'; % 12 hours
% file_end = '6_1800'; % 15 hours
% file_end = '6_2100'; % 18 hours
% file_end = '7_0000'; % 21 hours
% file_end = '7_0300'; % 24 hours
 
%Load Free (background) State.
cd /discover/nobackup/drholdaw/atmp.22292/prog/prog_free/
file = ['x0011dh_a.prog.eta.2013011',file_end,'z.nc4'];

lon = ncread(file,'lon');
lat = ncread(file,'lat');
lev = ncread(file,'lev');

u_free = ncread(file,'u');
v_free = ncread(file,'v');
t_free = ncread(file,'tv');
q_free = ncread(file,'sphu');
qi_free = ncread(file,'qitot');
ql_free = ncread(file,'qltot');

%Load Perturbed (analysis) st
cd /discover/nobackup/drholdaw/atmp.22292/prog/prog_replay10/
file = ['x0011dh_a.prog.eta.2013011',file_end,'z.nc4'];

u_replay = ncread(file,'u');
v_replay = ncread(file,'v');
t_replay = ncread(file,'tv');
q_replay = ncread(file,'sphu');
qi_replay = ncread(file,'qitot');
ql_replay = ncread(file,'qltot');

%Load first TLM state to compare.
cd /discover/nobackup/drholdaw/atmp.22293/sens.20130117.000000%/fvpert_
file = ['x0011dh_a.fvpert.eta.2013011',file_end,'z_BL2_MOI0_GWD0_GQ1_TF900.nc4'];
% file = 'fvpert.eta.nc4'

u_tlmdry = ncread(file,'u');
v_tlmdry = ncread(file,'v');
t_tlmdry = ncread(file,'tv');
q_tlmdry = ncread(file,'sphu');
qi_tlmdry = ncread(file,'qitot');
ql_tlmdry = ncread(file,'qltot');

%Load second TLM state to compare.
cd /discover/nobackup/drholdaw/atmp.22293/sens.20130117.000000%/fvpert_old/
file = ['x0011dh_a.fvpert.eta.2013011',file_end,'z_BL2_MOI2_GWD0_GQ1_TF900_CFTANH20.nc4'];
% file = 'fvpert.eta.nc4'

u_tlmmoi = ncread(file,'u');
v_tlmmoi = ncread(file,'v');
t_tlmmoi = ncread(file,'tv');
q_tlmmoi = ncread(file,'sphu');
qi_tlmmoi = ncread(file,'qitot');
ql_tlmmoi = ncread(file,'qltot');

cd /home/drholdaw/Lin_Moist_Physics/Bac_paper_figs/

%Compute NL perturbation trajectory.
u_nlm = u_replay - u_free;
v_nlm = v_replay - v_free;
t_nlm = t_replay - t_free;
q_nlm = q_replay - q_free;
qi_nlm = qi_replay - qi_free;
ql_nlm = ql_replay - ql_free;

%Pick out just the level to be plotted.
u_a_lev = u_nlm(:,:,plot_level);
v_a_lev = v_nlm(:,:,plot_level);
t_a_lev = t_nlm(:,:,plot_level);
q_a_lev = q_nlm(:,:,plot_level_q);
qi_a_lev = qi_nlm(:,:,plot_level_qi);
ql_a_lev = ql_nlm(:,:,plot_level_ql);

%Uncomment here to plot actual TLM pert trajectories
% u_b_lev = u_tlmdry(:,:,plot_level);
% v_b_lev = v_tlmdry(:,:,plot_level);
% t_b_lev = t_tlmdry(:,:,plot_level);
% q_b_lev = q_tlmdry(:,:,plot_level_q);
% qi_b_lev = qi_tlmdry(:,:,plot_level_qi);
% ql_b_lev = ql_tlmdry(:,:,plot_level_ql);
% 
% u_c_lev = u_tlmmoi(:,:,plot_level);
% v_c_lev = v_tlmmoi(:,:,plot_level);
% t_c_lev = t_tlmmoi(:,:,plot_level);
% q_c_lev = q_tlmmoi(:,:,plot_level_q);
% qi_c_lev = qi_tlmmoi(:,:,plot_level_qi);
% ql_c_lev = ql_tlmmoi(:,:,plot_level_ql);

% % %Uncomment here to plot the difference between the NL and TLM pert trajectories
u_b_lev = u_tlmdry(:,:,plot_level) - u_nlm(:,:,plot_level);
v_b_lev = v_tlmdry(:,:,plot_level) - v_nlm(:,:,plot_level);
t_b_lev = t_tlmdry(:,:,plot_level) - t_nlm(:,:,plot_level);
q_b_lev = q_tlmdry(:,:,plot_level_q) - q_nlm(:,:,plot_level_q);
qi_b_lev = qi_tlmdry(:,:,plot_level_qi) - qi_nlm(:,:,plot_level_qi);
ql_b_lev = ql_tlmdry(:,:,plot_level_ql) - ql_nlm(:,:,plot_level_ql);

u_c_lev = u_tlmmoi(:,:,plot_level) - u_nlm(:,:,plot_level);
v_c_lev = v_tlmmoi(:,:,plot_level) - v_nlm(:,:,plot_level);
t_c_lev = t_tlmmoi(:,:,plot_level) - t_nlm(:,:,plot_level);
q_c_lev = q_tlmmoi(:,:,plot_level_q) - q_nlm(:,:,plot_level_q);
qi_c_lev = qi_tlmmoi(:,:,plot_level_qi) - qi_nlm(:,:,plot_level_qi);
ql_c_lev = ql_tlmmoi(:,:,plot_level_ql) - ql_nlm(:,:,plot_level_ql);


lim = 6.9009e-04;
for i = 1:length(lon)
    for j = 1:length(lat)
        
        if abs(ql_a_lev(i,j)) >= lim
            ql_a_lev(i,j) = ql_a_lev(i,j)/abs(ql_a_lev(i,j)) * lim;
        else
            ql_a_lev(i,j) = ql_a_lev(i,j);
        end
        
        if abs(ql_b_lev(i,j)) >= lim
            ql_b_lev(i,j) = ql_b_lev(i,j)/abs(ql_b_lev(i,j)) * lim;
        else
            ql_b_lev(i,j) = ql_b_lev(i,j);
        end
        
        if abs(ql_c_lev(i,j)) >= lim
            ql_c_lev(i,j) = ql_c_lev(i,j)/abs(ql_c_lev(i,j)) * lim;
        else
            ql_c_lev(i,j) = ql_c_lev(i,j);
        end
        
    end
end

%Compute maximum absolute values for contour limits and intervals.
maxu = max([max(abs(u_a_lev(:))) max(abs(u_b_lev(:)))  max(abs(u_c_lev(:)))]);
minu = -maxu; 
maxv = max([max(abs(v_a_lev(:))) max(abs(v_b_lev(:)))  max(abs(v_c_lev(:)))]);
minv = -maxv;
maxt = max([max(abs(t_a_lev(:))) max(abs(t_b_lev(:)))  max(abs(t_c_lev(:)))]);
mint = -maxt;
maxq = max([max(abs(q_a_lev(:))) max(abs(q_b_lev(:)))  max(abs(q_c_lev(:)))]);
minq = -maxq;
maxqi = max([max(abs(qi_a_lev(:))) max(abs(qi_b_lev(:))) max(abs(qi_c_lev(:)))]);
minqi = -maxqi;
maxql = max([max(abs(ql_a_lev(:))) max(abs(ql_b_lev(:))) max(abs(ql_c_lev(:)))]);
minql = -maxql;

%Set contour intervals for each variable
cint_u = 2*maxu/15;
cint_v = 2*maxv/10;
cint_t = 2*maxt/10;
cint_q = 2*maxq/10;
cint_qi = 2*maxqi/20;
cint_ql = 2*maxql/20;

%Make center of colormap white
cmap = colormap;
cmap(32,:) = [1 1 1];
cmap(33,:) = [1 1 1];
close

grey = 0.75;
fontsize = 11;
line_wid_cont = 1.0;
line_wid_det = 0.6;


    figure('Position',[ 164 70 1058 849])
    
    subplot(3,2,1)
    [C,h] = contour(lon,lat,ql_a_lev','LineWidth',line_wid_cont);
    hold on; 
    plot(coast_lon,coast_lat,'Color',[grey grey grey])

    caxis([minql maxql])
    set(h,'LevelStep',cint_ql);

    colormap(cmap)
    pos = get(gca,'Position');
    colbar = colorbar;
    set(gca,'Position',[pos(1)-0.05 pos(2) pos(3) pos(4)])
    title(colbar,'(kgkg^{-1})','FontSize',fontsize,'FontName','TimesNewRoman')

    axis equal; box on; xlim([-180 180]); ylim([-90 90])
    set(gca,'XTick',[-180 -120 -60 0 60 120 180],'Ytick',[-90 -60 -30 0 30 60 90]);
    set(gca,'XTickLabel',{'180W', '120W', '60W', '0', '60E', '120E', '180E'}...
    ,'YtickLabel',{'90S' '60S' '30S' '0' '30N' '60N' '90N'});

    title('(a) q_l Nonlinear Perturbation Trajectory','FontSize',fontsize,'FontName','TimesNewRoman')
    set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')


    subplot(3,2,3)
    [C,h] = contour(lon,lat,ql_b_lev','LineWidth',line_wid_cont);
    hold on; 
    plot(coast_lon,coast_lat,'Color',[grey grey grey])

    caxis([minql maxql])
    set(h,'LevelStep',cint_ql);

    colormap(cmap)
    pos = get(gca,'Position');
    colbar = colorbar;
    set(gca,'Position',[pos(1)-0.05 pos(2) pos(3) pos(4)])
    title(colbar,'(kgkg^{-1})','FontSize',fontsize,'FontName','TimesNewRoman')

    axis equal; box on; xlim([-180 180]); ylim([-90 90])
    set(gca,'XTick',[-180 -120 -60 0 60 120 180],'Ytick',[-90 -60 -30 0 30 60 90]);
    set(gca,'XTickLabel',{'180W', '120W', '60W', '0', '60E', '120E', '180E'}...
    ,'YtickLabel',{'90S' '60S' '30S' '0' '30N' '60N' '90N'});

    title('(c) q_l Difference (Dry Physics)','FontSize',fontsize,'FontName','TimesNewRoman')
    set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')


    subplot(3,2,5)
    [C,h] = contour(lon,lat,ql_c_lev','LineWidth',line_wid_cont);
    hold on; 
    plot(coast_lon,coast_lat,'Color',[grey grey grey])

    caxis([minql maxql])
    set(h,'LevelStep',cint_ql);

    colormap(cmap)
    pos = get(gca,'Position');
    colbar = colorbar;
    set(gca,'Position',[pos(1)-0.05 pos(2) pos(3) pos(4)])
    title(colbar,'(kgkg^{-1})','FontSize',fontsize,'FontName','TimesNewRoman')

    axis equal; box on; xlim([-180 180]); ylim([-90 90])
    set(gca,'XTick',[-180 -120 -60 0 60 120 180],'Ytick',[-90 -60 -30 0 30 60 90]);
    set(gca,'XTickLabel',{'180W', '120W', '60W', '0', '60E', '120E', '180E'}...
    ,'YtickLabel',{'90S' '60S' '30S' '0' '30N' '60N' '90N'});

    title('(e) q_l Difference (Moist Physics)','FontSize',fontsize,'FontName','TimesNewRoman')
    set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
    
    
    
    
    
    subplot(3,2,2)
    [C,h] = contour(lon,lat,qi_a_lev','LineWidth',line_wid_cont);
    hold on; 
    plot(coast_lon,coast_lat,'Color',[grey grey grey])

    caxis([minqi maxqi])
    set(h,'LevelStep',cint_qi);

    colormap(cmap)
    pos = get(gca,'Position');
    colbar = colorbar;
    set(gca,'Position',[pos(1) pos(2) pos(3) pos(4)])
    title(colbar,'(kgkg^{-1})','FontSize',fontsize,'FontName','TimesNewRoman')

    axis equal; box on; xlim([-180 180]); ylim([-90 90])
    set(gca,'XTick',[-180 -120 -60 0 60 120 180],'Ytick',[-90 -60 -30 0 30 60 90]);
    set(gca,'XTickLabel',{'180W', '120W', '60W', '0', '60E', '120E', '180E'}...
    ,'YtickLabel',{'90S' '60S' '30S' '0' '30N' '60N' '90N'});

    title('(b) q_i Nonlinear Perturbation Trajectory','FontSize',fontsize,'FontName','TimesNewRoman')
    set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')


    subplot(3,2,4)
    [C,h] = contour(lon,lat,qi_b_lev','LineWidth',line_wid_cont);
    hold on; 
    plot(coast_lon,coast_lat,'Color',[grey grey grey])

    caxis([minqi maxqi])
    set(h,'LevelStep',cint_qi);

    colormap(cmap)
    pos = get(gca,'Position');
    colbar = colorbar;
    set(gca,'Position',[pos(1) pos(2) pos(3) pos(4)])
    title(colbar,'(kgkg^{-1})','FontSize',fontsize,'FontName','TimesNewRoman')

    axis equal; box on; xlim([-180 180]); ylim([-90 90])
    set(gca,'XTick',[-180 -120 -60 0 60 120 180],'Ytick',[-90 -60 -30 0 30 60 90]);
    set(gca,'XTickLabel',{'180W', '120W', '60W', '0', '60E', '120E', '180E'}...
    ,'YtickLabel',{'90S' '60S' '30S' '0' '30N' '60N' '90N'});

    title('(d) q_i Difference (Dry Physics)','FontSize',fontsize,'FontName','TimesNewRoman')
    set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')


    subplot(3,2,6)
    [C,h] = contour(lon,lat,qi_c_lev','LineWidth',line_wid_cont);
    hold on; 
    plot(coast_lon,coast_lat,'Color',[grey grey grey])

    caxis([minqi maxqi])
    set(h,'LevelStep',cint_qi);

    colormap(cmap)
    pos = get(gca,'Position');
    colbar = colorbar;
    set(gca,'Position',[pos(1) pos(2) pos(3) pos(4)])
    title(colbar,'(kgkg^{-1})','FontSize',fontsize,'FontName','TimesNewRoman')

    axis equal; box on; xlim([-180 180]); ylim([-90 90])
    set(gca,'XTick',[-180 -120 -60 0 60 120 180],'Ytick',[-90 -60 -30 0 30 60 90]);
    set(gca,'XTickLabel',{'180W', '120W', '60W', '0', '60E', '120E', '180E'}...
    ,'YtickLabel',{'90S' '60S' '30S' '0' '30N' '60N' '90N'});

    title('(f) q_i Difference (Moist Physics)','FontSize',fontsize,'FontName','TimesNewRoman')
    set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
    





