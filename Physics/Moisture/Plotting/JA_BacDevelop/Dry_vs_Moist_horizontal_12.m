close all
clear
clc

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

%Choose whether to compute the correlations.
calc_corr = 0;

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

cd /home/drholdaw/Lin_Moist_Physics/Inputs/
pref
 
%Load Free (background) State.
cd /discover/nobackup/drholdaw/tmp.22292/prog/prog_free/
file = ['x0011dh_a.prog.eta.2013011',file_end,'z.nc4'];

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

%Load Perturbed (analysis) st
cd /discover/nobackup/drholdaw/tmp.22292/prog/prog_replay10/
file = ['x0011dh_a.prog.eta.2013011',file_end,'z.nc4'];

u_replay = ncread(file,'u');
v_replay = ncread(file,'v');
t_replay = ncread(file,'tv');
q_replay = ncread(file,'sphu');
p_replay = ncread(file,'delp');
qi_replay = ncread(file,'qitot');
ql_replay = ncread(file,'qltot');
o3_replay = ncread(file,'ozone');

%Load first TLM state to compare.
cd /discover/nobackup/drholdaw/tmp.22293/sens.20130117.000000%/fvpert_
file = ['x0011dh_a.fvpert.eta.2013011',file_end,'z_BL2_MOI0_GQ0_TF900.nc4']
% file = 'fvpert.eta.nc4'

u_tlmdry = ncread(file,'u');
v_tlmdry = ncread(file,'v');
t_tlmdry = ncread(file,'tv');
q_tlmdry = ncread(file,'sphu');
p_tlmdry = ncread(file,'delp');
qi_tlmdry = ncread(file,'qitot');
ql_tlmdry = ncread(file,'qltot');
o3_tlmdry = ncread(file,'ozone');

%Load second TLM state to compare.
cd /discover/nobackup/drholdaw/tmp.22293/sens.20130117.000000%/fvpert_old/
file = ['x0011dh_a.fvpert.eta.2013011',file_end,'z_BL2_MOI2_GQ0_TF900_FILT24.nc4']
% file = 'fvpert.eta.nc4'

u_tlmmoi = ncread(file,'u');
v_tlmmoi = ncread(file,'v');
t_tlmmoi = ncread(file,'tv');
q_tlmmoi = ncread(file,'sphu');
p_tlmmoi = ncread(file,'delp');
qi_tlmmoi = ncread(file,'qitot');
ql_tlmmoi = ncread(file,'qltot');
o3_tlmmoi = ncread(file,'ozone');

cd /home/drholdaw/Lin_Moist_Physics/DAS_plotting/Cloud' Fraction'/

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
qi_a_lev = qi_nlm(:,:,plot_level);
ql_a_lev = ql_nlm(:,:,plot_level);
o3_a_lev = o3_nlm(:,:,plot_level);

%Uncomment here to plot actual TLM pert trajectories
u_b_lev = u_tlmdry(:,:,plot_level);
v_b_lev = v_tlmdry(:,:,plot_level);
t_b_lev = t_tlmdry(:,:,plot_level);
q_b_lev = q_tlmdry(:,:,plot_level);
p_b_lev = p_tlmdry(:,:,plot_level);
qi_b_lev = qi_tlmdry(:,:,plot_level);
ql_b_lev = ql_tlmdry(:,:,plot_level);
o3_b_lev = o3_tlmdry(:,:,plot_level);

u_c_lev = u_tlmmoi(:,:,plot_level);
v_c_lev = v_tlmmoi(:,:,plot_level);
t_c_lev = t_tlmmoi(:,:,plot_level);
q_c_lev = q_tlmmoi(:,:,plot_level);
p_c_lev = p_tlmmoi(:,:,plot_level);
qi_c_lev = qi_tlmmoi(:,:,plot_level);
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
maxu = max([max(abs(u_a_lev(:))) max(abs(u_b_lev(:)))  max(abs(u_c_lev(:)))]);
minu = -maxu; 
maxv = max([max(abs(v_a_lev(:))) max(abs(v_b_lev(:)))  max(abs(v_c_lev(:)))]);
minv = -maxv;
maxt = max([max(abs(t_a_lev(:))) max(abs(t_b_lev(:)))  max(abs(t_c_lev(:)))]);
mint = -maxt;
maxq = max([max(abs(q_a_lev(:))) max(abs(q_b_lev(:)))  max(abs(q_c_lev(:)))]);
minq = -maxq;
maxp = max([max(abs(p_a_lev(:))) max(abs(p_b_lev(:)))  max(abs(p_c_lev(:)))]);
minp = -maxp; 
maxqi = max([max(abs(qi_a_lev(:))) max(abs(qi_b_lev(:))) max(abs(qi_c_lev(:)))]);
minqi = -maxqi;
maxql = max([max(abs(ql_a_lev(:))) max(abs(ql_b_lev(:))) max(abs(qi_c_lev(:)))]);
minql = -maxql;
maxo3 = max([max(abs(o3_a_lev(:))) max(abs(o3_b_lev(:))) max(abs(qi_c_lev(:)))]);
mino3 = -maxo3;

%Set contour intervals for each variable
cint_u = 2*maxu/10;
cint_v = 2*maxv/10;
cint_t = 2*maxt/10;
cint_q = 2*maxq/10;
cint_p = 2*maxp/10;
cint_qi = 2*maxqi/10;
cint_ql = 2*maxql/10;
cint_o3 = 2*maxo3/10;

%Make figures of pert trajectories.
if makeplots == 1 

    if plot_u == 1
        
        figure
        set(gcf,'position',[621 30 658 889])

        subplot(3,1,1)
        [C,h] = contourf(lon,lat,((u_a_lev))','LineStyle','none');
        caxis([minu maxu])
        title('u')
        colorbar
        set(h,'LevelStep',cint_u);

        subplot(3,1,2)
        [C,h] = contourf(lon,lat,((u_b_lev))','LineStyle','none');
        caxis([minu maxu])
        title('u')
        colorbar
        set(h,'LevelStep',cint_u);

        subplot(3,1,3)
        [C,h] = contourf(lon,lat,((u_c_lev))','LineStyle','none');
        caxis([minu maxu])
        title('v')
        colorbar
        set(h,'LevelStep',cint_u);
    
    end
        
    if plot_v == 1
    
        figure
        set(gcf,'position',[621 30 658 889])

        subplot(3,1,1)
        [C,h] = contourf(lon,lat,((v_a_lev))','LineStyle','none');
        caxis([minv maxv])
        title('v')
        colorbar
        set(h,'LevelStep',cint_v);

        subplot(3,1,2)
        [C,h] = contourf(lon,lat,((v_b_lev))','LineStyle','none');
        caxis([minv maxv])
        title('v')
        colorbar
        set(h,'LevelStep',cint_v);

        subplot(3,1,3)
        [C,h] = contourf(lon,lat,((v_c_lev))','LineStyle','none');
        caxis([minv maxv])
        title('v')
        colorbar
        set(h,'LevelStep',cint_v);
    
    end
        
    if plot_t == 1

        figure
        set(gcf,'position',[621 30 658 889])

        subplot(3,1,1)
        [C,h] = contourf(lon,lat,((t_a_lev))','LineStyle','none');
        caxis([mint maxt])
        title('t')
        colorbar
        set(h,'LevelStep',cint_t);

        subplot(3,1,2)
        [C,h] = contourf(lon,lat,((t_b_lev))','LineStyle','none');
        caxis([mint maxt])
        title('t')
        colorbar
        set(h,'LevelStep',cint_t);

        subplot(3,1,3)
        [C,h] = contourf(lon,lat,((t_c_lev))','LineStyle','none');
        caxis([mint maxt])
        title('t')
        colorbar
        set(h,'LevelStep',cint_t);
    
    end
        
    if plot_q == 1

        figure
        set(gcf,'position',[621 30 658 889])

        subplot(3,1,1)
        [C,h] = contourf(lon,lat,((q_a_lev))','LineStyle','none');
        caxis([minq maxq])
        title('q')
        colorbar
        set(h,'LevelStep',cint_q);

        subplot(3,1,2)
        [C,h] = contourf(lon,lat,((q_b_lev))','LineStyle','none');
        caxis([minq maxq])
        title('q')
        colorbar
        set(h,'LevelStep',cint_q);

        subplot(3,1,3)
        [C,h] = contourf(lon,lat,((q_c_lev))','LineStyle','none');
        caxis([minq maxq])
        title('q')
        colorbar
        set(h,'LevelStep',cint_q);
    
    end
        
    if plot_p == 1

        figure
        set(gcf,'position',[621 30 658 889])

        subplot(3,1,1)
        [C,h] = contourf(lon,lat,((p_a_lev))','LineStyle','none');
        caxis([minp maxp])
        title('p')
        colorbar
        set(h,'LevelStep',cint_p);

        subplot(3,1,2)
        [C,h] = contourf(lon,lat,((p_b_lev))','LineStyle','none');
        caxis([minp maxp])
        title('p')
        colorbar
        set(h,'LevelStep',cint_p);

        subplot(3,1,3)
        [C,h] = contourf(lon,lat,((p_c_lev))','LineStyle','none');
        caxis([minp maxp])
        title('p')
        colorbar
        set(h,'LevelStep',cint_p);
    
    end
        
    if plot_qi == 1

        figure
        set(gcf,'position',[621 30 658 889])

        subplot(3,1,1)
        [C,h] = contourf(lon,lat,((qi_a_lev))','LineStyle','none');
        caxis([minqi maxqi])
        title('qi')
        colorbar
        set(h,'LevelStep',cint_qi);

        subplot(3,1,2)
        [C,h] = contourf(lon,lat,((qi_b_lev))','LineStyle','none');
        caxis([minqi maxqi])
        title('qi')
        colorbar
        set(h,'LevelStep',cint_qi);

        subplot(3,1,3)
        [C,h] = contourf(lon,lat,((qi_c_lev))','LineStyle','none');
        caxis([minqi maxqi])
        title('qi')
        colorbar
        set(h,'LevelStep',cint_qi);
    
    end
        
    if plot_ql == 1

        figure
        set(gcf,'position',[621 30 658 889])

        subplot(3,1,1)
        [C,h] = contourf(lon,lat,((ql_a_lev))','LineStyle','none');
        caxis([minql maxql])    
        title('ql')
        colorbar
        set(h,'LevelStep',cint_ql/0.5);

        subplot(3,1,2)
        [C,h] = contourf(lon,lat,((ql_b_lev))','LineStyle','none');
        caxis([minql maxql])    
        title('ql')
        colorbar
        set(h,'LevelStep',cint_ql/0.5);

        subplot(3,1,3)
        [C,h] = contourf(lon,lat,((ql_c_lev))','LineStyle','none');
        caxis([minql maxql])    
        title('ql')
        colorbar
        set(h,'LevelStep',cint_ql/0.5);
    
    end
        
    if plot_o3 == 1

        figure
        set(gcf,'position',[621 30 658 889])
    
        subplot(3,1,1)
        [C,h] = contourf(lon,lat,((o3_a_lev))','LineStyle','none');
        caxis([mino3 maxo3])
        title('o3')
        colorbar
        set(h,'LevelStep',cint_o3);
        
        subplot(3,1,2)
        [C,h] = contourf(lon,lat,((o3_b_lev))','LineStyle','none');
        caxis([mino3 maxo3])
        title('o3')
        colorbar
        set(h,'LevelStep',cint_o3);
            
        subplot(3,1,3)
        [C,h] = contourf(lon,lat,((o3_c_lev))','LineStyle','none');
        caxis([mino3 maxo3])
        title('o3')
        colorbar
        set(h,'LevelStep',cint_o3);
        
    end
    
    
end

