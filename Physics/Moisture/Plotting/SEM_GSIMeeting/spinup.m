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

cd /home/drholdaw/Lin_Moist_Physics/Inputs/
pref
 
%Load Free (background) State.
cd /discover/nobackup/drholdaw/tmp.22292a/prog/prog_free/
file = ['x0011dh_a.prog.eta.20130116_00z.nc4'];

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
cd /discover/nobackup/drholdaw/tmp.22292a/prog/prog_replay10/
file = ['x0011dh_a.prog.eta.20130116_00z.nc4'];

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
file = ['fvpert.eta.nc4']

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
file = ['fvpert.eta.nc4']

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
    
    cov_dry = zeros(lm,8);
    
    corr_dry = zeros(lm,8);
    
    %Loop over model all levels.
    for k = 1:lm

        dA = 0.0;

        if k >= lev_min
        
            %Loop over domain of interest
            for i = lon_min:lon_max
                for j = lat_min:lat_max

                    dA = dA + A(i,j);

                    cov_dry(k,1) = cov_dry(k,1) + u_nlm(i,j,k)*u_tlmdry(i,j,k)*A(i,j);

                    var_nlm(k,1) = var_nlm(k,1) + u_nlm(i,j,k)*u_nlm(i,j,k)*A(i,j);
                    var_dry(k,1) = var_dry(k,1) + u_tlmdry(i,j,k)*u_tlmdry(i,j,k)*A(i,j);

                    cov_dry(k,2) = cov_dry(k,2) + v_nlm(i,j,k)*v_tlmdry(i,j,k)*A(i,j);

                    var_nlm(k,2) = var_nlm(k,2) + v_nlm(i,j,k)*v_nlm(i,j,k)*A(i,j);
                    var_dry(k,2) = var_dry(k,2) + v_tlmdry(i,j,k)*v_tlmdry(i,j,k)*A(i,j);

                    cov_dry(k,3) = cov_dry(k,3) + t_nlm(i,j,k)*t_tlmdry(i,j,k)*A(i,j);

                    var_nlm(k,3) = var_nlm(k,3) + t_nlm(i,j,k)*t_nlm(i,j,k)*A(i,j);
                    var_dry(k,3) = var_dry(k,3) + t_tlmdry(i,j,k)*t_tlmdry(i,j,k)*A(i,j);

                    cov_dry(k,4) = cov_dry(k,4) + q_nlm(i,j,k)*q_tlmdry(i,j,k)*A(i,j);

                    var_nlm(k,4) = var_nlm(k,4) + q_nlm(i,j,k)*q_nlm(i,j,k)*A(i,j);
                    var_dry(k,4) = var_dry(k,4) + q_tlmdry(i,j,k)*q_tlmdry(i,j,k)*A(i,j);

                    cov_dry(k,5) = cov_dry(k,5) + p_nlm(i,j,k)*p_tlmdry(i,j,k)*A(i,j);

                    var_nlm(k,5) = var_nlm(k,5) + p_nlm(i,j,k)*p_nlm(i,j,k)*A(i,j);
                    var_dry(k,5) = var_dry(k,5) + p_tlmdry(i,j,k)*p_tlmdry(i,j,k)*A(i,j);

                    cov_dry(k,6) = cov_dry(k,6) + qi_nlm(i,j,k)*qi_tlmdry(i,j,k)*A(i,j);

                    var_nlm(k,6) = var_nlm(k,6) + qi_nlm(i,j,k)*qi_nlm(i,j,k)*A(i,j);
                    var_dry(k,6) = var_dry(k,6) + qi_tlmdry(i,j,k)*qi_tlmdry(i,j,k)*A(i,j);

                    cov_dry(k,7) = cov_dry(k,7) + ql_nlm(i,j,k)*ql_tlmdry(i,j,k)*A(i,j);

                    var_nlm(k,7) = var_nlm(k,7) + ql_nlm(i,j,k)*ql_nlm(i,j,k)*A(i,j);
                    var_dry(k,7) = var_dry(k,7) + ql_tlmdry(i,j,k)*ql_tlmdry(i,j,k)*A(i,j);

                    cov_dry(k,8) = cov_dry(k,8) + o3_nlm(i,j,k)*o3_tlmdry(i,j,k)*A(i,j);

                    var_nlm(k,8) = var_nlm(k,8) + o3_nlm(i,j,k)*o3_nlm(i,j,k)*A(i,j);
                    var_dry(k,8) = var_dry(k,8) + o3_tlmdry(i,j,k)*o3_tlmdry(i,j,k)*A(i,j);

                end
            end

            cov_dry(k,1) = cov_dry(k,1)/dA;
            var_nlm(k,1) = var_nlm(k,1)/dA;
            var_dry(k,1) = var_dry(k,1)/dA;

            cov_dry(k,2) = cov_dry(k,2)/dA;
            var_nlm(k,2) = var_nlm(k,2)/dA;
            var_dry(k,2) = var_dry(k,2)/dA;

            cov_dry(k,3) = cov_dry(k,3)/dA;
            var_nlm(k,3) = var_nlm(k,3)/dA;
            var_dry(k,3) = var_dry(k,3)/dA;

            cov_dry(k,4) = cov_dry(k,4)/dA;
            var_nlm(k,4) = var_nlm(k,4)/dA;
            var_dry(k,4) = var_dry(k,4)/dA;

            cov_dry(k,5) = cov_dry(k,5)/dA;
            var_nlm(k,5) = var_nlm(k,5)/dA;
            var_dry(k,5) = var_dry(k,5)/dA;

            cov_dry(k,6) = cov_dry(k,6)/dA;
            var_nlm(k,6) = var_nlm(k,6)/dA;
            var_dry(k,6) = var_dry(k,6)/dA;

            cov_dry(k,7) = cov_dry(k,7)/dA;
            var_nlm(k,7) = var_nlm(k,7)/dA;
            var_dry(k,7) = var_dry(k,7)/dA;

            cov_dry(k,8) = cov_dry(k,8)/dA;
            var_nlm(k,8) = var_nlm(k,8)/dA;
            var_dry(k,8) = var_dry(k,8)/dA;

            corr_dry(k,1) = cov_dry(k,1)/sqrt(var_nlm(k,1)*var_dry(k,1));

            corr_dry(k,2) = cov_dry(k,2)/sqrt(var_nlm(k,2)*var_dry(k,2));

            corr_dry(k,3) = cov_dry(k,3)/sqrt(var_nlm(k,3)*var_dry(k,3));

            corr_dry(k,4) = cov_dry(k,4)/sqrt(var_nlm(k,4)*var_dry(k,4));

            corr_dry(k,5) = cov_dry(k,5)/sqrt(var_nlm(k,5)*var_dry(k,5));

            corr_dry(k,6) = cov_dry(k,6)/sqrt(var_nlm(k,6)*var_dry(k,6));

            corr_dry(k,7) = cov_dry(k,7)/sqrt(var_nlm(k,7)*var_dry(k,7));

            corr_dry(k,8) = cov_dry(k,8)/sqrt(var_nlm(k,8)*var_dry(k,8));

        end
            
    end


    
end

figure
set(gcf,'position',[3 343 1276 576])

subplot(1,8,1)
plot(corr_dry(lev_min:lev_max,1),p_ref(lev_min+1:lev_max+1),'LineWidth',linewid)
set(gca,'YDir','reverse')
xlim([0 1])
ylim([p_ref(lev_min+1) p_ref(lev_max+1)])
title('u','FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
ylabel('p (hPa)','FontSize',fontsize,'FontName','TimesNewRoman')

subplot(1,8,2)
plot(corr_dry(lev_min:lev_max,2),p_ref(lev_min+1:lev_max+1),'LineWidth',linewid)
set(gca,'YDir','reverse')
xlim([0 1])
ylim([p_ref(lev_min+1) p_ref(lev_max+1)])
title('v','FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'YTickLabel',[])

subplot(1,8,3)
plot(corr_dry(lev_min:lev_max,3),p_ref(lev_min+1:lev_max+1),'LineWidth',linewid)
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
set(gca,'YDir','reverse')
xlim([0 1])
ylim([p_ref(lev_min+1) p_ref(lev_max+1)])
title('q','FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'YTickLabel',[])

subplot(1,8,5)
plot(corr_dry(lev_min:lev_max,5),p_ref(lev_min+1:lev_max+1),'LineWidth',linewid)
set(gca,'YDir','reverse')
xlim([0 1])
ylim([p_ref(lev_min+1) p_ref(lev_max+1)])
title('\delta p','FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'YTickLabel',[])

subplot(1,8,6)
plot(corr_dry(lev_min:lev_max,6),p_ref(lev_min+1:lev_max+1),'LineWidth',linewid)
set(gca,'YDir','reverse')
xlim([-0.2 1])
ylim([p_ref(lev_min+1) p_ref(lev_max+1)])
title('qi','FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'YTickLabel',[])

subplot(1,8,7)
plot(corr_dry(lev_min:lev_max,7),p_ref(lev_min+1:lev_max+1),'LineWidth',linewid)
set(gca,'YDir','reverse')
xlim([-0.2 1])
ylim([p_ref(lev_min+1) p_ref(lev_max+1)])
title('ql','FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'YTickLabel',[])

subplot(1,8,8)
plot(corr_dry(lev_min:lev_max,8),p_ref(lev_min+1:lev_max+1),'LineWidth',linewid)
set(gca,'YDir','reverse')
xlim([0 1])
ylim([p_ref(lev_min+1) p_ref(lev_max+1)])
title('o_3','FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'YAxisLocation','right')

