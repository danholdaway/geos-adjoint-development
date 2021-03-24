close all
clear
clc

%Choose model level to plot.
plot_level = 63;

fontsize = 13;

line_wid_cont = 1.0;
line_wid_det = 0.6;

linewid = 1.5;

grey = 0.75;

%Choose whether to compute the correlations.
calc_corr = 1;

%Choose Region,
region = 1;

%Choose whether you want plots of pert traj (0/1) and which model level to plot 
makeplots = 0;
plot_level = 50;

plot_u = 0;
plot_v = 0;
plot_t = 1;
plot_q = 0;
plot_p = 0;
plot_qi = 0;
plot_ql = 0;
plot_o3 = 0;

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
% file_end = '6_0900';
% file_end = '6_1200';
% file_end = '6_1500';
% file_end = '6_1800';
% file_end = '6_2100';
file_end = '7_0000';

cd /home/drholdaw/Lin_Moist_Physics/Inputs/
pref
 
%Load Free (background) State.
cd /discover/nobackup/drholdaw/tmp.22292/prog/prog_free/
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
cd /discover/nobackup/drholdaw/tmp.22292/prog/prog_replay10/
file = ['x0011dh_a.prog.eta.2013011',file_end(1:4),'z.nc4'];

u_replay = ncread(file,'u');
v_replay = ncread(file,'v');
t_replay = ncread(file,'tv');
q_replay = ncread(file,'sphu');
p_replay = ncread(file,'delp');
qi_replay = ncread(file,'qitot');
ql_replay = ncread(file,'qltot');
o3_replay = ncread(file,'ozone');

%Load first TLM state to compare.
cd /discover/nobackup/drholdaw/tmp.22292/sens.20130117.000000%/fvpert_old/
file = ['x0011dh_a.fvpert.eta.2013011',file_end,'z_newBL_03rst.nc4']

u_tlmdry = ncread(file,'u');
v_tlmdry = ncread(file,'v');
t_tlmdry = ncread(file,'tv');
q_tlmdry = ncread(file,'sphu');
p_tlmdry = ncread(file,'delp');
qi_tlmdry = ncread(file,'qitot');
ql_tlmdry = ncread(file,'qltot');
o3_tlmdry = ncread(file,'ozone');

%Load second TLM state to compare.
cd /discover/nobackup/drholdaw/tmp.22292/sens.20130117.000000%/fvpert_old/
file = ['x0011dh_a.fvpert.eta.2013011',file_end,'z_newMOInewBL_03rst.nc4']

u_tlmmoi = ncread(file,'u');
v_tlmmoi = ncread(file,'v');
t_tlmmoi = ncread(file,'tv');
q_tlmmoi = ncread(file,'sphu');
p_tlmmoi = ncread(file,'delp');
qi_tlmmoi = ncread(file,'qitot');
ql_tlmmoi = ncread(file,'qltot');
o3_tlmmoi = ncread(file,'ozone');

cd /home/drholdaw/Lin_Moist_Physics/DAS_plotting/GSI' Talk'/

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
        set(h,'LevelStep',cint_ql);

        subplot(3,1,2)
        [C,h] = contourf(lon,lat,((ql_b_lev))','LineStyle','none');
        caxis([minql maxql])    
        title('ql')
        colorbar
        set(h,'LevelStep',cint_ql);

        subplot(3,1,3)
        [C,h] = contourf(lon,lat,((ql_c_lev))','LineStyle','none');
        caxis([minql maxql])    
        title('ql')
        colorbar
        set(h,'LevelStep',cint_ql);
    
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

%Compute correlations.
if calc_corr == 1

    im = length(lon);
    jm = length(lat);
    lm = length(lev);
    
    %Weighting 
    r_earth = 6378.1;
    lat_pole = lat+lat(1)*-1;
    A = zeros(length(lon),length(lat));
    A(1,:) = (2*pi^2*r_earth^2/(length(lon)*length(lat))) * sind(lat_pole);
    for i = 2:length(lon)
        A(i,:) = A(1,:); %Make into rank-2 array for simplicity
    end

    dA = 0.0;

    if region == 1
        lon_min = 1;
        lon_max = im;
        lat_min = 1;
        lat_max = jm;
        lev_min = 30;
        lev_max = lm;
    elseif region == 2
        lon_min = 1;
        lon_max = im;
        lat_min = 121;
        lat_max = 241;
        lev_min = 30;
        lev_max = lm;
    elseif region == 3
        lon_min = 1;
        lon_max = im;
        lat_min = ceil(jm/2);
        lat_max = jm;
        lev_min = 30;
        lev_max = lm;
    elseif region == 4
        lon_min = 1;
        lon_max = im;
        lat_min = 1;
        lat_max = floor(jm/2);
        lev_min = 30;
        lev_max = lm;    
    end
    
    %Variance, covariance and correlation
    %Second index is number of variables.
    var_nlm = zeros(lm,8);
    var_dry = zeros(lm,8);
    var_moi = zeros(lm,8);
    
    cov_dry = zeros(lm,8);
    cov_moi = zeros(lm,8);
    
    corr_dry = zeros(lm,8);
    corr_moi = zeros(lm,8);
    
    %Loop over model all levels.
    for k = 1:lm

        dA = 0.0;

        if k >= lev_min
        
            %Loop over domain of interest
            for i = lon_min:lon_max
                for j = lat_min:lat_max

                    dA = dA + A(i,j);

                    cov_dry(k,1) = cov_dry(k,1) + u_nlm(i,j,k)*u_tlmdry(i,j,k)*A(i,j);
                    cov_moi(k,1) = cov_moi(k,1) + u_nlm(i,j,k)*u_tlmmoi(i,j,k)*A(i,j);

                    var_nlm(k,1) = var_nlm(k,1) + u_nlm(i,j,k)*u_nlm(i,j,k)*A(i,j);
                    var_dry(k,1) = var_dry(k,1) + u_tlmdry(i,j,k)*u_tlmdry(i,j,k)*A(i,j);
                    var_moi(k,1) = var_moi(k,1) + u_tlmmoi(i,j,k)*u_tlmmoi(i,j,k)*A(i,j);

                    cov_dry(k,2) = cov_dry(k,2) + v_nlm(i,j,k)*v_tlmdry(i,j,k)*A(i,j);
                    cov_moi(k,2) = cov_moi(k,2) + v_nlm(i,j,k)*v_tlmmoi(i,j,k)*A(i,j);

                    var_nlm(k,2) = var_nlm(k,2) + v_nlm(i,j,k)*v_nlm(i,j,k)*A(i,j);
                    var_dry(k,2) = var_dry(k,2) + v_tlmdry(i,j,k)*v_tlmdry(i,j,k)*A(i,j);
                    var_moi(k,2) = var_moi(k,2) + v_tlmmoi(i,j,k)*v_tlmmoi(i,j,k)*A(i,j);

                    cov_dry(k,3) = cov_dry(k,3) + t_nlm(i,j,k)*t_tlmdry(i,j,k)*A(i,j);
                    cov_moi(k,3) = cov_moi(k,3) + t_nlm(i,j,k)*t_tlmmoi(i,j,k)*A(i,j);

                    var_nlm(k,3) = var_nlm(k,3) + t_nlm(i,j,k)*t_nlm(i,j,k)*A(i,j);
                    var_dry(k,3) = var_dry(k,3) + t_tlmdry(i,j,k)*t_tlmdry(i,j,k)*A(i,j);
                    var_moi(k,3) = var_moi(k,3) + t_tlmmoi(i,j,k)*t_tlmmoi(i,j,k)*A(i,j);

                    cov_dry(k,4) = cov_dry(k,4) + q_nlm(i,j,k)*q_tlmdry(i,j,k)*A(i,j);
                    cov_moi(k,4) = cov_moi(k,4) + q_nlm(i,j,k)*q_tlmmoi(i,j,k)*A(i,j);

                    var_nlm(k,4) = var_nlm(k,4) + q_nlm(i,j,k)*q_nlm(i,j,k)*A(i,j);
                    var_dry(k,4) = var_dry(k,4) + q_tlmdry(i,j,k)*q_tlmdry(i,j,k)*A(i,j);
                    var_moi(k,4) = var_moi(k,4) + q_tlmmoi(i,j,k)*q_tlmmoi(i,j,k)*A(i,j);

                    cov_dry(k,5) = cov_dry(k,5) + p_nlm(i,j,k)*p_tlmdry(i,j,k)*A(i,j);
                    cov_moi(k,5) = cov_moi(k,5) + p_nlm(i,j,k)*p_tlmmoi(i,j,k)*A(i,j);

                    var_nlm(k,5) = var_nlm(k,5) + p_nlm(i,j,k)*p_nlm(i,j,k)*A(i,j);
                    var_dry(k,5) = var_dry(k,5) + p_tlmdry(i,j,k)*p_tlmdry(i,j,k)*A(i,j);
                    var_moi(k,5) = var_moi(k,5) + p_tlmmoi(i,j,k)*p_tlmmoi(i,j,k)*A(i,j);

                    cov_dry(k,6) = cov_dry(k,6) + qi_nlm(i,j,k)*qi_tlmdry(i,j,k)*A(i,j);
                    cov_moi(k,6) = cov_moi(k,6) + qi_nlm(i,j,k)*qi_tlmmoi(i,j,k)*A(i,j);

                    var_nlm(k,6) = var_nlm(k,6) + qi_nlm(i,j,k)*qi_nlm(i,j,k)*A(i,j);
                    var_dry(k,6) = var_dry(k,6) + qi_tlmdry(i,j,k)*qi_tlmdry(i,j,k)*A(i,j);
                    var_moi(k,6) = var_moi(k,6) + qi_tlmmoi(i,j,k)*qi_tlmmoi(i,j,k)*A(i,j);

                    cov_dry(k,7) = cov_dry(k,7) + ql_nlm(i,j,k)*ql_tlmdry(i,j,k)*A(i,j);
                    cov_moi(k,7) = cov_moi(k,7) + ql_nlm(i,j,k)*ql_tlmmoi(i,j,k)*A(i,j);

                    var_nlm(k,7) = var_nlm(k,7) + ql_nlm(i,j,k)*ql_nlm(i,j,k)*A(i,j);
                    var_dry(k,7) = var_dry(k,7) + ql_tlmdry(i,j,k)*ql_tlmdry(i,j,k)*A(i,j);
                    var_moi(k,7) = var_moi(k,7) + ql_tlmmoi(i,j,k)*ql_tlmmoi(i,j,k)*A(i,j);

                    cov_dry(k,8) = cov_dry(k,8) + o3_nlm(i,j,k)*o3_tlmdry(i,j,k)*A(i,j);
                    cov_moi(k,8) = cov_moi(k,8) + o3_nlm(i,j,k)*o3_tlmmoi(i,j,k)*A(i,j);

                    var_nlm(k,8) = var_nlm(k,8) + o3_nlm(i,j,k)*o3_nlm(i,j,k)*A(i,j);
                    var_dry(k,8) = var_dry(k,8) + o3_tlmdry(i,j,k)*o3_tlmdry(i,j,k)*A(i,j);
                    var_moi(k,8) = var_moi(k,8) + o3_tlmmoi(i,j,k)*o3_tlmmoi(i,j,k)*A(i,j);

                end
            end

            cov_dry(k,1) = cov_dry(k,1)/dA;
            cov_moi(k,1) = cov_moi(k,1)/dA;
            var_nlm(k,1) = var_nlm(k,1)/dA;
            var_dry(k,1) = var_dry(k,1)/dA;
            var_moi(k,1) = var_moi(k,1)/dA;

            cov_dry(k,2) = cov_dry(k,2)/dA;
            cov_moi(k,2) = cov_moi(k,2)/dA;
            var_nlm(k,2) = var_nlm(k,2)/dA;
            var_dry(k,2) = var_dry(k,2)/dA;
            var_moi(k,2) = var_moi(k,2)/dA;

            cov_dry(k,3) = cov_dry(k,3)/dA;
            cov_moi(k,3) = cov_moi(k,3)/dA;
            var_nlm(k,3) = var_nlm(k,3)/dA;
            var_dry(k,3) = var_dry(k,3)/dA;
            var_moi(k,3) = var_moi(k,3)/dA;

            cov_dry(k,4) = cov_dry(k,4)/dA;
            cov_moi(k,4) = cov_moi(k,4)/dA;
            var_nlm(k,4) = var_nlm(k,4)/dA;
            var_dry(k,4) = var_dry(k,4)/dA;
            var_moi(k,4) = var_moi(k,4)/dA;

            cov_dry(k,5) = cov_dry(k,5)/dA;
            cov_moi(k,5) = cov_moi(k,5)/dA;
            var_nlm(k,5) = var_nlm(k,5)/dA;
            var_dry(k,5) = var_dry(k,5)/dA;
            var_moi(k,5) = var_moi(k,5)/dA;

            cov_dry(k,6) = cov_dry(k,6)/dA;
            cov_moi(k,6) = cov_moi(k,6)/dA;
            var_nlm(k,6) = var_nlm(k,6)/dA;
            var_dry(k,6) = var_dry(k,6)/dA;
            var_moi(k,6) = var_moi(k,6)/dA;

            cov_dry(k,7) = cov_dry(k,7)/dA;
            cov_moi(k,7) = cov_moi(k,7)/dA;
            var_nlm(k,7) = var_nlm(k,7)/dA;
            var_dry(k,7) = var_dry(k,7)/dA;
            var_moi(k,7) = var_moi(k,7)/dA;

            cov_dry(k,8) = cov_dry(k,8)/dA;
            cov_moi(k,8) = cov_moi(k,8)/dA;
            var_nlm(k,8) = var_nlm(k,8)/dA;
            var_dry(k,8) = var_dry(k,8)/dA;
            var_moi(k,8) = var_moi(k,8)/dA;

            corr_dry(k,1) = cov_dry(k,1)/sqrt(var_nlm(k,1)*var_dry(k,1));
            corr_moi(k,1) = cov_moi(k,1)/sqrt(var_nlm(k,1)*var_moi(k,1));

            corr_dry(k,2) = cov_dry(k,2)/sqrt(var_nlm(k,2)*var_dry(k,2));
            corr_moi(k,2) = cov_moi(k,2)/sqrt(var_nlm(k,2)*var_moi(k,2));

            corr_dry(k,3) = cov_dry(k,3)/sqrt(var_nlm(k,3)*var_dry(k,3));
            corr_moi(k,3) = cov_moi(k,3)/sqrt(var_nlm(k,3)*var_moi(k,3));

            corr_dry(k,4) = cov_dry(k,4)/sqrt(var_nlm(k,4)*var_dry(k,4));
            corr_moi(k,4) = cov_moi(k,4)/sqrt(var_nlm(k,4)*var_moi(k,4));

            corr_dry(k,5) = cov_dry(k,5)/sqrt(var_nlm(k,5)*var_dry(k,5));
            corr_moi(k,5) = cov_moi(k,5)/sqrt(var_nlm(k,5)*var_moi(k,5));

            corr_dry(k,6) = cov_dry(k,6)/sqrt(var_nlm(k,6)*var_dry(k,6));
            corr_moi(k,6) = cov_moi(k,6)/sqrt(var_nlm(k,6)*var_moi(k,6));

            corr_dry(k,7) = cov_dry(k,7)/sqrt(var_nlm(k,7)*var_dry(k,7));
            corr_moi(k,7) = cov_moi(k,7)/sqrt(var_nlm(k,7)*var_moi(k,7));

            corr_dry(k,8) = cov_dry(k,8)/sqrt(var_nlm(k,8)*var_dry(k,8));
            corr_moi(k,8) = cov_moi(k,8)/sqrt(var_nlm(k,8)*var_moi(k,8));

        end
            
    end

    format short

    %Show the correlations for Dry, moist and then difference.
    fprintf('u wind \n')
    fprintf('     Level       Dry     Moist      Gain \n')
    disp([lev(lev_min:lev_max) corr_dry(lev_min:lev_max,1) corr_moi(lev_min:lev_max,1) corr_moi(lev_min:lev_max,1)-corr_dry(lev_min:lev_max,1)])
    fprintf('v wind \n')
    fprintf('     Level       Dry     Moist      Gain \n')
    disp([lev(lev_min:lev_max) corr_dry(lev_min:lev_max,2) corr_moi(lev_min:lev_max,2) corr_moi(lev_min:lev_max,2)-corr_dry(lev_min:lev_max,2)])
    fprintf('Temperature \n')
    fprintf('     Level       Dry     Moist      Gain \n')
    disp([lev(lev_min:lev_max) corr_dry(lev_min:lev_max,3) corr_moi(lev_min:lev_max,3) corr_moi(lev_min:lev_max,3)-corr_dry(lev_min:lev_max,3)])
    fprintf('Specific Humidity \n')
    fprintf('     Level       Dry     Moist      Gain \n')
    disp([lev(lev_min:lev_max) corr_dry(lev_min:lev_max,4) corr_moi(lev_min:lev_max,4) corr_moi(lev_min:lev_max,4)-corr_dry(lev_min:lev_max,4)])
    fprintf('Pressure Thickness \n')
    fprintf('     Level       Dry     Moist      Gain \n')
    disp([lev(lev_min:lev_max) corr_dry(lev_min:lev_max,5) corr_moi(lev_min:lev_max,5) corr_moi(lev_min:lev_max,5)-corr_dry(lev_min:lev_max,5)])
    fprintf('qi Cloud \n')
    fprintf('     Level       Dry     Moist      Gain \n')
    disp([lev(lev_min:lev_max) corr_dry(lev_min:lev_max,6) corr_moi(lev_min:lev_max,6) corr_moi(lev_min:lev_max,6)-corr_dry(lev_min:lev_max,6)])
    fprintf('ql Cloud \n')
    fprintf('     Level       Dry     Moist      Gain \n')
    disp([lev(lev_min:lev_max) corr_dry(lev_min:lev_max,7) corr_moi(lev_min:lev_max,7) corr_moi(lev_min:lev_max,7)-corr_dry(lev_min:lev_max,7)])
    fprintf('Ozone \n')
    fprintf('     Level       Dry     Moist      Gain \n')
    disp([lev(lev_min:lev_max) corr_dry(lev_min:lev_max,8) corr_moi(lev_min:lev_max,8) corr_moi(lev_min:lev_max,8)-corr_dry(lev_min:lev_max,8)])

    % disp([(1:72)' lev corr_dry(:,1)-corr_moi(:,1) corr_dry(:,2)-corr_moi(:,2) corr_dry(:,3)-corr_moi(:,3) corr_dry(:,4)-corr_moi(:,4)])

    %High up correlations for clouds and pressure can include divide by zero so ignore here
    for i = lm:-1:lev_min
        if isnan(corr_dry(i,5)) == 0
           pnan_dry = i;
        end
        if isnan(corr_dry(i,6)) == 0
           qinan_dry = i;
        end
        if isnan(corr_dry(i,7)) == 0
           qlnan_dry = i;
        end
        if isnan(corr_moi(i,5)) == 0
           pnan_moi = i;
        end
        if isnan(corr_moi(i,6)) == 0
           qinan_moi = i;
        end
        if isnan(corr_moi(i,7)) == 0
           qlnan_moi = i;
        end
    end

    pnan = max([pnan_dry pnan_moi]);
    qinan = max([qinan_dry qinan_moi]);
    qlnan = max([qlnan_dry qlnan_moi]);

    %Show mean correlation for the levels of interest.
    mean_u_corr = [ mean(corr_dry(lev_min:lev_max,1)) mean(corr_moi(lev_min:lev_max,1)) ];
    fprintf('Mean u  corr. Dry: %f  Moist: %f  Gain: %f \n\n', mean_u_corr(1), mean_u_corr(2), mean_u_corr(2)-mean_u_corr(1));
    mean_v_corr = [ mean(corr_dry(lev_min:lev_max,2)) mean(corr_moi(lev_min:lev_max,2)) ];
    fprintf('Mean v  corr. Dry: %f  Moist: %f  Gain: %f \n\n', mean_v_corr(1), mean_v_corr(2), mean_v_corr(2)-mean_v_corr(1));
    mean_t_corr = [ mean(corr_dry(lev_min:lev_max,3)) mean(corr_moi(lev_min:lev_max,3)) ];
    fprintf('Mean t  corr. Dry: %f  Moist: %f  Gain: %f \n\n', mean_t_corr(1), mean_t_corr(2), mean_t_corr(2)-mean_t_corr(1));
    mean_q_corr = [ mean(corr_dry(lev_min:lev_max,4)) mean(corr_moi(lev_min:lev_max,4)) ];
    fprintf('Mean q  corr. Dry: %f  Moist: %f  Gain: %f <-- \n\n', mean_q_corr(1), mean_q_corr(2), mean_q_corr(2)-mean_q_corr(1));
    mean_p_corr = [ mean(corr_dry(pnan:lm,5)) mean(corr_moi(pnan:lm,5)) ];
    fprintf('Mean p  corr. Dry: %f  Moist: %f  Gain: %f \n\n', mean_p_corr(1), mean_p_corr(2), mean_p_corr(2)-mean_p_corr(1));
    mean_qi_corr = [ mean(corr_dry(qinan:lm,6)) mean(corr_moi(qinan:lm,6)) ];
    fprintf('Mean qi corr. Dry: %f  Moist: %f  Gain: %f <-- \n\n', mean_qi_corr(1), mean_qi_corr(2), mean_qi_corr(2)-mean_qi_corr(1));
    mean_ql_corr = [ mean(corr_dry(qlnan:lm,7)) mean(corr_moi(qlnan:lm,7)) ];
    fprintf('Mean ql corr. Dry: %f  Moist: %f  Gain: %f <-- \n\n', mean_ql_corr(1), mean_ql_corr(2), mean_ql_corr(2)-mean_ql_corr(1));
    mean_o3_corr = [ mean(corr_dry(:,8)) mean(corr_moi(:,8)) ];
    fprintf('Mean o3 corr. Dry: %f  Moist: %f  Gain: %f \n\n', mean_o3_corr(1), mean_o3_corr(2), mean_o3_corr(2)-mean_o3_corr(1));
    
end

figure
set(gcf,'position',[3 343 1276 576])

subplot(1,8,1)
plot(corr_dry(lev_min:lev_max,1),p_ref(lev_min+1:lev_max+1),'LineWidth',linewid)
hold on
plot(corr_moi(lev_min:lev_max,1),p_ref(lev_min+1:lev_max+1),'r','LineWidth',linewid)
set(gca,'YDir','reverse')
xlim([0 1])
ylim([p_ref(lev_min+1) p_ref(lev_max+1)])
title('u','FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
ylabel('p (hPa)','FontSize',fontsize,'FontName','TimesNewRoman')

subplot(1,8,2)
plot(corr_dry(lev_min:lev_max,2),p_ref(lev_min+1:lev_max+1),'LineWidth',linewid)
hold on
plot(corr_moi(lev_min:lev_max,2),p_ref(lev_min+1:lev_max+1),'r','LineWidth',linewid)
set(gca,'YDir','reverse')
xlim([0 1])
ylim([p_ref(lev_min+1) p_ref(lev_max+1)])
title('v','FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'YTickLabel',[])

subplot(1,8,3)
plot(corr_dry(lev_min:lev_max,3),p_ref(lev_min+1:lev_max+1),'LineWidth',linewid)
hold on
plot(corr_moi(lev_min:lev_max,3),p_ref(lev_min+1:lev_max+1),'r','LineWidth',linewid)
hold on
set(gca,'YDir','reverse')
xlim([0 1])
ylim([p_ref(lev_min+1) p_ref(lev_max+1)])
title('T_v','FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'YTickLabel',[])

subplot(1,8,4)
plot(corr_dry(lev_min:lev_max,4),p_ref(lev_min+1:lev_max+1),'LineWidth',linewid)
hold on
plot(corr_moi(lev_min:lev_max,4),p_ref(lev_min+1:lev_max+1),'r','LineWidth',linewid)
hold on
set(gca,'YDir','reverse')
xlim([0 1])
ylim([p_ref(lev_min+1) p_ref(lev_max+1)])
title('q','FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'YTickLabel',[])

subplot(1,8,5)
plot(corr_dry(lev_min:lev_max,5),p_ref(lev_min+1:lev_max+1),'LineWidth',linewid)
hold on
plot(corr_moi(lev_min:lev_max,5),p_ref(lev_min+1:lev_max+1),'r','LineWidth',linewid)
set(gca,'YDir','reverse')
xlim([0 1])
ylim([p_ref(lev_min+1) p_ref(lev_max+1)])
title('\delta p','FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'YTickLabel',[])

subplot(1,8,6)
plot(corr_dry(lev_min:lev_max,6),p_ref(lev_min+1:lev_max+1),'LineWidth',linewid)
hold on
plot(corr_moi(lev_min:lev_max,6),p_ref(lev_min+1:lev_max+1),'r','LineWidth',linewid)
set(gca,'YDir','reverse')
xlim([0 1])
ylim([p_ref(lev_min+1) p_ref(lev_max+1)])
title('qi','FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'YTickLabel',[])

subplot(1,8,7)
plot(corr_dry(lev_min:lev_max,7),p_ref(lev_min+1:lev_max+1),'LineWidth',linewid)
hold on
plot(corr_moi(lev_min:lev_max,7),p_ref(lev_min+1:lev_max+1),'r','LineWidth',linewid)
set(gca,'YDir','reverse')
xlim([0 1])
ylim([p_ref(lev_min+1) p_ref(lev_max+1)])
title('ql','FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'YTickLabel',[])

subplot(1,8,8)
plot(corr_dry(lev_min:lev_max,8),p_ref(lev_min+1:lev_max+1),'LineWidth',linewid)
hold on
plot(corr_moi(lev_min:lev_max,8),p_ref(lev_min+1:lev_max+1),'r','LineWidth',linewid)
set(gca,'YDir','reverse')
xlim([0 1])
ylim([p_ref(lev_min+1) p_ref(lev_max+1)])
title('o_3','FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'YAxisLocation','right')