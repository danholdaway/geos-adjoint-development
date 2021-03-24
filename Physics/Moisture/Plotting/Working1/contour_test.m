clear
close all
clc

%Choose whether you want plots of pert traj (0/1) and which model level to plot 
makeplots = 1;
plot_level_qi = 38;
plot_level_ql = 53;

%Choose whether to compute the correlations.
calc_corr = 1;

%Choose Region,
% 1 = Global
% 2 = Tropics 30S to 30N, below 100hPa
% 3 = Northern hemisphere
% 4 = Southern hemisphere
region = 1;

%Choose the time being considered: 16th at 12z would be 6_12
% file_end = '6_0300';
% file_end = '6_0600';
% file_end = '6_1200';
% file_end = '6_1500';
 file_end = '7_0000';

%Load Free (background) State.
cd /discover/nobackup/drholdaw/tmp.22292/nlm_runs/freeClouds/
file = ['x0011dh_a.prog.eta.2013011',file_end(1:4),'z.nc4'];

lon = ncread(file,'lon');
lat = ncread(file,'lat');
lev = ncread(file,'lev');

qi_free = ncread(file,'qitot');
ql_free = ncread(file,'qltot');

%Load Perturbed (analysis) state.
cd /discover/nobackup/drholdaw/tmp.22292/nlm_runs/replayClouds/
file = ['x0011dh_a.prog.eta.2013011',file_end(1:4),'z.nc4'];

qi_replay = ncread(file,'qitot');
ql_replay = ncread(file,'qltot');

%Load first TLM state to compare.
cd /discover/nobackup/drholdaw/tmp.22292/sens.20130117.000000%/fvpert_old/
file = ['x0011dh_a.fvpert.eta.2013011',file_end,'z_dry.nc4'];

qi_tlmdry = ncread(file,'qitot');
ql_tlmdry = ncread(file,'qltot');

%Load second TLM state to compare.
cd /discover/nobackup/drholdaw/tmp.22292/sens.20130117.000000%/fvpert_old/
file = ['x0011dh_a.fvpert.eta.2013011',file_end,'z_moi3.nc4'];

qi_tlmmoi = ncread(file,'qitot');
ql_tlmmoi = ncread(file,'qltot');

cd /home/drholdaw/Lin_Moist_Physics/DAS_plotting/Cloud' Fraction'/

%Compute NL perturbation trajectory.
qi_nlm = qi_replay - qi_free;
ql_nlm = ql_replay - ql_free;

%Pick out just the level to be plotted.
qi_a_lev = qi_nlm(:,:,plot_level_qi);
ql_a_lev = ql_nlm(:,:,plot_level_ql);

%Uncomment here to plot actual TLM pert trajectories
qi_b_lev = qi_tlmdry(:,:,plot_level_qi);
ql_b_lev = ql_tlmdry(:,:,plot_level_ql);

qi_c_lev = qi_tlmmoi(:,:,plot_level_qi);
ql_c_lev = ql_tlmmoi(:,:,plot_level_ql);

% %Uncomment here to plot the difference between the NL and TLM pert trajectories
% qi_b_lev = qi_tlmdry(:,:,plot_level) - qi_nlm(:,:,plot_level);
% ql_b_lev = ql_tlmdry(:,:,plot_level) - ql_nlm(:,:,plot_level);
% 
% qi_c_lev = qi_tlmmoi(:,:,plot_level) - qi_nlm(:,:,plot_level);
% ql_c_lev = ql_tlmmoi(:,:,plot_level) - ql_nlm(:,:,plot_level);


%Compute maximum absolute values for contour limits and intervals.
maxqi = max([max((log(abs(qi_a_lev(:))))) max((log(abs(qi_b_lev(:))))) max((log(abs(qi_c_lev(:)))))])
minqi = -100 %min([min((log(abs(qi_a_lev(:))))) min((log(abs(qi_b_lev(:))))) min((log(abs(qi_c_lev(:)))))])
maxql = max([max((log(abs(ql_a_lev(:))))) max((log(abs(ql_b_lev(:))))) max((log(abs(qi_c_lev(:)))))])
minql = -100 %min([min((log(abs(ql_a_lev(:))))) min((log(abs(ql_b_lev(:))))) min((log(abs(qi_c_lev(:)))))])

%Set contour intervals for each variable
cint_qi = 2*maxqi/10;
cint_ql = 2*maxql/10;

%Make figures of pert trajectories.
if makeplots == 1 

    figure
    set(gcf,'position',[621 30 658 889])

    subplot(3,1,1)
    [C,h] = contourf(lon,lat,log(abs(qi_a_lev))','LineStyle','none');
    caxis([minql maxql])    
    title('qi')
    colorbar
    set(h,'LevelStep',cint_ql);
    
    subplot(3,1,2)
    [C,h] = contourf(lon,lat,log(abs(qi_b_lev))','LineStyle','none');
    caxis([minql maxql])    
    title('qi')
    colorbar
    set(h,'LevelStep',cint_ql);
        
    subplot(3,1,3)
    [C,h] = contourf(lon,lat,log(abs(qi_c_lev))','LineStyle','none');
    caxis([minql maxql])    
    title('qi')
    colorbar
    set(h,'LevelStep',cint_ql);
    
    figure
    set(gcf,'position',[621 30 658 889])

    subplot(3,1,1)
    [C,h] = contourf(lon,lat,log(abs(ql_a_lev))','LineStyle','none');
    caxis([minql maxql])    
    title('ql')
    colorbar
    set(h,'LevelStep',cint_ql);
    
    subplot(3,1,2)
    [C,h] = contourf(lon,lat,log(abs(ql_b_lev))','LineStyle','none');
    caxis([minql maxql])    
    title('ql')
    colorbar
    set(h,'LevelStep',cint_ql);
        
    subplot(3,1,3)
    [C,h] = contourf(lon,lat,log(abs(ql_c_lev))','LineStyle','none');
    caxis([minql maxql])    
    title('ql')
    colorbar
    set(h,'LevelStep',cint_ql);

end
