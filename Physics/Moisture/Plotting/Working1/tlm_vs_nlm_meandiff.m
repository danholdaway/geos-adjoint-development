close all
clear
clc
mydir = pwd;


%Choose whether to compute the correlations.
calc_rms = 1;

%Choose Region,
% 1 = Global
% 2 = Tropics 23S to 23N, below 100hPa
% 3 = Northern hemisphere
% 4 = Southern hemisphere
region = 1;

if region == 1 
    fprintf('Computing RMSEs globally \n')
elseif region == 2
    fprintf('Computing RMSEs for the tropics \n')  
elseif region == 3
    fprintf('Computing RMSEs for the northern hemisphere \n')
elseif region == 4
    fprintf('Computing RMSEs for the southern hemisphere \n')
else
    fprintf('Not a valid region selection \n')
end

%Choose the time being considered: 16th at 12z would be 6_12
% file_end = '1_0100';
% file_end = '1_0300';
% file_end = '1_0600';
% file_end = '1_0900';
% file_end = '1_1200';
file_end = '1_1500';
% file_end = '1_1800';
% file_end = '1_2100';
% file_end = '2_0000';
% file_end = '2_0300';

% file_end = '2_0300'; % 24 hours

im = 576;
jm = 361;
lm = 72;

LevMin = 30;
LevMax = lm;

if region == 1
    LonMin = 1;
    LonMax = im;
    LatMin = 1;
    LatMax = jm;
elseif region == 2
    LonMin = 1;
    LonMax = im;
    LatMin = 135;
    LatMax = 227;
elseif region == 3
    LonMin = 1;
    LonMax = im;
    LatMin = ceil(jm/2);
    LatMax = jm;
elseif region == 4
    LonMin = 1;
    LonMax = im;
    LatMin = 1;
    LatMax = floor(jm/2);
end

cd /home/drholdaw/LinearisedPhysics/Inputs/
pref

%Load Free (background) State.
cd /discover/nobackup/drholdaw/btmp.25956/prog/prog_free_onlyremaptracer/
file = ['v000_C180.prog.eta.2014020',file_end,'z.nc4'];

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
cd /discover/nobackup/drholdaw/btmp.25956/prog/prog_replay_onlyremaptracer/
file = ['v000_C180.prog.eta.2014020',file_end,'z.nc4'];

u_replay = ncread(file,'u');
v_replay = ncread(file,'v');
t_replay = ncread(file,'tv');
q_replay = ncread(file,'sphu');
p_replay = ncread(file,'delp');
qi_replay = ncread(file,'qitot');
ql_replay = ncread(file,'qltot');
o3_replay = ncread(file,'ozone');

%Load first TLM state to compare.
cd /discover/nobackup/drholdaw/btmp.25956/sens.20140202.000000/
% file = ['v000_C180.fvpert.eta.2014020',file_end,'z_03z_P200000.nc4'];
file = 'fvpertX.eta.nc4';

u_tlma = ncread(file,'U');
v_tlma = ncread(file,'V');
t_tlma = ncread(file,'TV');
q_tlma = ncread(file,'QV');
p_tlma = ncread(file,'DP');
qi_tlma = ncread(file,'QI');
ql_tlma = ncread(file,'QL');
o3_tlma = ncread(file,'O3');

%Load second TLM state to compare.
cd /discover/nobackup/drholdaw/btmp.25956/sens.20140202.000000/
file = ['v000_C180.fvpert.eta.2014020',file_end,'z_03z_ONLYREMAPTRACER.nc4'];
% file = 'fvpertX.eta.nc4';

u_tlmb = ncread(file,'U');
v_tlmb = ncread(file,'V');
t_tlmb = ncread(file,'TV');
q_tlmb = ncread(file,'QV');
p_tlmb = ncread(file,'DP');
qi_tlmb = ncread(file,'QI');
ql_tlmb = ncread(file,'QL');
o3_tlmb = ncread(file,'O3');


%Load second TLM state to compare.
cd /discover/nobackup/drholdaw/btmp.25956/sens.20140202.000000/
file = ['v000_C180.fvpert.eta.2014020',file_end,'z_03z_ONLYREMAPTRACER.nc4'];
% file = 'fvpertX.eta.nc4';

u_tlmc = ncread(file,'U');
v_tlmc = ncread(file,'V');
t_tlmc = ncread(file,'TV');
q_tlmc = ncread(file,'QV');
p_tlmc = ncread(file,'DP');
qi_tlmc = ncread(file,'QI');
ql_tlmc = ncread(file,'QL');
o3_tlmc = ncread(file,'O3');

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


%Compute correlations.
if calc_rms == 1

    im = length(lon);
    jm = length(lat);
    lm = length(lev);
    
    %Weighting 
    r_earth = 6378.1;
    lat_pole = lat+lat(1)*-1;
    A = ones(im,jm,lm);
%     for k = 1:lm
%         A(1,:,k) = (2*pi^2*r_earth^2/(length(lon)*length(lat))) * sind(lat_pole);
%         for i = 2:length(lon)
%             A(i,:,k) = A(1,:,k); %Make into rank-2 array for simplicity
%         end
%     end
        
    A = 1.0;
    
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
        lat_min = 135;
        lat_max = 227;
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
    
    
    diff2_u_a = abs(u_nlm - u_tlma);
    diff2_u_b = abs(u_nlm - u_tlmb);
    diff2_u_c = abs(u_nlm - u_tlmc);
    diff2_v_a = abs(v_nlm - v_tlma);
    diff2_v_b = abs(v_nlm - v_tlmb);
    diff2_v_c = abs(v_nlm - v_tlmc);
    diff2_t_a = abs(t_nlm - t_tlma);
    diff2_t_b = abs(t_nlm - t_tlmb);
    diff2_t_c = abs(t_nlm - t_tlmc);
    diff2_q_a = abs(q_nlm - q_tlma);
    diff2_q_b = abs(q_nlm - q_tlmb);
    diff2_q_c = abs(q_nlm - q_tlmc);
    diff2_p_a = abs(p_nlm - p_tlma);
    diff2_p_b = abs(p_nlm - p_tlmb);
    diff2_p_c = abs(p_nlm - p_tlmc);
    diff2_qi_a = abs(qi_nlm - qi_tlma);
    diff2_qi_b = abs(qi_nlm - qi_tlmb);
    diff2_qi_c = abs(qi_nlm - qi_tlmc);
    diff2_ql_a = abs(ql_nlm - ql_tlma);
    diff2_ql_b = abs(ql_nlm - ql_tlmb);
    diff2_ql_c = abs(ql_nlm - ql_tlmc);
    diff2_o3_a = abs(o3_nlm - o3_tlma);
    diff2_o3_b = abs(o3_nlm - o3_tlmb);
    diff2_o3_c = abs(o3_nlm - o3_tlmc);
    
    sumA = sum(sum(A(:,:,1)));
    
    rms_u_a = zeros(lm,1);
    rms_u_b = zeros(lm,1);
    rms_u_c = zeros(lm,1);
    rms_v_a = zeros(lm,1);
    rms_v_b = zeros(lm,1);
    rms_v_c = zeros(lm,1);
    rms_t_a = zeros(lm,1);
    rms_t_b = zeros(lm,1);
    rms_t_c = zeros(lm,1);
    rms_q_a = zeros(lm,1);
    rms_q_b = zeros(lm,1);
    rms_q_c = zeros(lm,1);
    rms_p_a = zeros(lm,1);
    rms_p_b = zeros(lm,1);
    rms_p_c = zeros(lm,1);
    rms_qi_a = zeros(lm,1);
    rms_qi_b = zeros(lm,1);
    rms_qi_c = zeros(lm,1);
    rms_ql_a = zeros(lm,1);
    rms_ql_b = zeros(lm,1);
    rms_ql_c = zeros(lm,1);
    rms_o3_a = zeros(lm,1);
    rms_o3_b = zeros(lm,1);
    rms_o3_c = zeros(lm,1);
    
    n = (lon_max-lon_min+1)*(lat_max-lat_min+1);
    
    %Loop over model all levels.
    for k = 1:lm

        rms_u_a(k) = ( sum(sum(diff2_u_a(lon_min:lon_max,lat_min:lat_max,k))) / n );
        rms_u_b(k) = ( sum(sum(diff2_u_b(lon_min:lon_max,lat_min:lat_max,k))) / n );
        rms_u_c(k) = ( sum(sum(diff2_u_c(lon_min:lon_max,lat_min:lat_max,k))) / n );

        rms_v_a(k) = ( sum(sum(diff2_v_a(lon_min:lon_max,lat_min:lat_max,k))) / n );
        rms_v_b(k) = ( sum(sum(diff2_v_b(lon_min:lon_max,lat_min:lat_max,k))) / n );
        rms_v_c(k) = ( sum(sum(diff2_v_c(lon_min:lon_max,lat_min:lat_max,k))) / n );
        
        rms_t_a(k) = ( sum(sum(diff2_t_a(lon_min:lon_max,lat_min:lat_max,k))) / n );
        rms_t_b(k) = ( sum(sum(diff2_t_b(lon_min:lon_max,lat_min:lat_max,k))) / n );
        rms_t_c(k) = ( sum(sum(diff2_t_c(lon_min:lon_max,lat_min:lat_max,k))) / n );
        
        rms_q_a(k) = ( sum(sum(diff2_q_a(lon_min:lon_max,lat_min:lat_max,k))) / n );
        rms_q_b(k) = ( sum(sum(diff2_q_b(lon_min:lon_max,lat_min:lat_max,k))) / n );
        rms_q_c(k) = ( sum(sum(diff2_q_c(lon_min:lon_max,lat_min:lat_max,k))) / n );
        
        rms_p_a(k) = ( sum(sum(diff2_p_a(lon_min:lon_max,lat_min:lat_max,k))) / n );
        rms_p_b(k) = ( sum(sum(diff2_p_b(lon_min:lon_max,lat_min:lat_max,k))) / n );
        rms_p_c(k) = ( sum(sum(diff2_p_c(lon_min:lon_max,lat_min:lat_max,k))) / n );

        rms_qi_a(k) = ( sum(sum(diff2_qi_a(lon_min:lon_max,lat_min:lat_max,k))) / n );
        rms_qi_b(k) = ( sum(sum(diff2_qi_b(lon_min:lon_max,lat_min:lat_max,k))) / n );
        rms_qi_c(k) = ( sum(sum(diff2_qi_c(lon_min:lon_max,lat_min:lat_max,k))) / n );
        
        rms_ql_a(k) = ( sum(sum(diff2_ql_a(lon_min:lon_max,lat_min:lat_max,k))) / n );
        rms_ql_b(k) = ( sum(sum(diff2_ql_b(lon_min:lon_max,lat_min:lat_max,k))) / n );
        rms_ql_c(k) = ( sum(sum(diff2_ql_c(lon_min:lon_max,lat_min:lat_max,k))) / n );
        
        rms_o3_a(k) = ( sum(sum(diff2_o3_a(lon_min:lon_max,lat_min:lat_max,k))) / n );
        rms_o3_b(k) = ( sum(sum(diff2_o3_b(lon_min:lon_max,lat_min:lat_max,k))) / n );
        rms_o3_c(k) = ( sum(sum(diff2_o3_c(lon_min:lon_max,lat_min:lat_max,k))) / n );
    end

    format short

    %Show the correlations for Dry, bst and then difference.
    fprintf('u wind \n')
    fprintf('     Level       Dry     Moist      Gain \n')
    disp([lev(lev_min:lev_max) rms_u_a(lev_min:lev_max) rms_u_b(lev_min:lev_max) rms_u_a(lev_min:lev_max)-rms_u_b(lev_min:lev_max)])
    fprintf('v wind \n')
    fprintf('     Level       Dry     Moist      Gain \n')
    disp([lev(lev_min:lev_max) rms_v_a(lev_min:lev_max) rms_v_b(lev_min:lev_max) rms_v_b(lev_min:lev_max)-rms_v_a(lev_min:lev_max)])
    fprintf('Temperature \n')
    fprintf('     Level       Dry     Moist      Gain \n')
    disp([lev(lev_min:lev_max) rms_t_a(lev_min:lev_max) rms_t_b(lev_min:lev_max) rms_t_b(lev_min:lev_max)-rms_t_a(lev_min:lev_max)])
    fprintf('Specific Humidity \n')
    fprintf('     Level       Dry     Moist      Gain \n')
    disp([lev(lev_min:lev_max) rms_q_a(lev_min:lev_max) rms_q_b(lev_min:lev_max) rms_q_b(lev_min:lev_max)-rms_q_a(lev_min:lev_max)])
    fprintf('Pressure Thickness \n')
    fprintf('     Level       Dry     Moist      Gain \n')
    disp([lev(lev_min:lev_max) rms_p_a(lev_min:lev_max) rms_p_b(lev_min:lev_max) rms_p_b(lev_min:lev_max)-rms_p_a(lev_min:lev_max)])
    fprintf('qi Cloud \n')
    fprintf('     Level       Dry     Moist      Gain \n')
    disp([lev(lev_min:lev_max) rms_qi_a(lev_min:lev_max) rms_qi_b(lev_min:lev_max) rms_qi_b(lev_min:lev_max)-rms_qi_a(lev_min:lev_max)])
    fprintf('ql Cloud \n')
    fprintf('     Level       Dry     Moist      Gain \n')
    disp([lev(lev_min:lev_max) rms_ql_a(lev_min:lev_max) rms_ql_b(lev_min:lev_max) rms_ql_b(lev_min:lev_max)-rms_ql_a(lev_min:lev_max)])
    fprintf('Ozone \n')
    fprintf('     Level       Dry     Moist      Gain \n')
    disp([lev(lev_min:lev_max) rms_o3_a(lev_min:lev_max) rms_o3_b(lev_min:lev_max) rms_o3_b(lev_min:lev_max)-rms_o3_a(lev_min:lev_max)])
    
    
% 
%     % disp([(1:72)' lev corr_a(:,1)-corr_b(:,1) corr_a(:,2)-corr_b(:,2) corr_a(:,3)-corr_b(:,3) corr_a(:,4)-corr_b(:,4)])
% 
%     %High up correlations for clouds and pressure can include divide by zero so ignore here
%     for i = lm:-1:lev_min
%         if isnan(corr_a(i,5)) == 0
%            pnan_a = i;
%         end
%         if isnan(corr_a(i,6)) == 0
%            qinan_a = i;
%         end
%         if isnan(corr_a(i,7)) == 0
%            qlnan_a = i;
%         end
%         if isnan(corr_b(i,5)) == 0
%            pnan_b = i;
%         end
%         if isnan(corr_b(i,6)) == 0
%            qinan_b = i;
%         end
%         if isnan(corr_b(i,7)) == 0
%            qlnan_b = i;
%         end
%     end
% 
%     pnan = max([pnan_a pnan_b]);
%     qinan = max([qinan_a qinan_b]);
%     qlnan = max([qlnan_a qlnan_b]);
% 
%     %Show mean correlation for the levels of interest.
%     mean_u_corr = [ mean(corr_a(lev_min:lev_max,1)) mean(corr_b(lev_min:lev_max,1)) ];
%     fprintf('Mean u  corr. Dry: %f  Moist: %f  Gain: %f \n\n', mean_u_corr(1), mean_u_corr(2), mean_u_corr(2)-mean_u_corr(1));
%     mean_v_corr = [ mean(corr_a(lev_min:lev_max,2)) mean(corr_b(lev_min:lev_max,2)) ];
%     fprintf('Mean v  corr. Dry: %f  Moist: %f  Gain: %f \n\n', mean_v_corr(1), mean_v_corr(2), mean_v_corr(2)-mean_v_corr(1));
%     mean_t_corr = [ mean(corr_a(lev_min:lev_max,3)) mean(corr_b(lev_min:lev_max,3)) ];
%     fprintf('Mean t  corr. Dry: %f  Moist: %f  Gain: %f \n\n', mean_t_corr(1), mean_t_corr(2), mean_t_corr(2)-mean_t_corr(1));
%     mean_q_corr = [ mean(corr_a(lev_min:lev_max,4)) mean(corr_b(lev_min:lev_max,4)) ];
%     fprintf('Mean q  corr. Dry: %f  Moist: %f  Gain: %f <-- \n\n', mean_q_corr(1), mean_q_corr(2), mean_q_corr(2)-mean_q_corr(1));
%     mean_p_corr = [ mean(corr_a(pnan:lm,5)) mean(corr_b(pnan:lm,5)) ];
%     fprintf('Mean p  corr. Dry: %f  Moist: %f  Gain: %f \n\n', mean_p_corr(1), mean_p_corr(2), mean_p_corr(2)-mean_p_corr(1));
%     mean_qi_corr = [ mean(corr_a(qinan:lm,6)) mean(corr_b(qinan:lm,6)) ];
%     fprintf('Mean qi corr. Dry: %f  Moist: %f  Gain: %f <-- \n\n', mean_qi_corr(1), mean_qi_corr(2), mean_qi_corr(2)-mean_qi_corr(1));
%     mean_ql_corr = [ mean(corr_a(qlnan:lm,7)) mean(corr_b(qlnan:lm,7)) ];
%     fprintf('Mean ql corr. Dry: %f  Moist: %f  Gain: %f <-- \n\n', mean_ql_corr(1), mean_ql_corr(2), mean_ql_corr(2)-mean_ql_corr(1));
%     mean_o3_corr = [ mean(corr_a(:,8)) mean(corr_b(:,8)) ];
%     fprintf('Mean o3 corr. Dry: %f  Moist: %f  Gain: %f \n\n', mean_o3_corr(1), mean_o3_corr(2), mean_o3_corr(2)-mean_o3_corr(1));
%     
% 
    figure
    set(gcf,'position',[3 343 1276 576])

    subplot(1,8,1)
    plot(rms_u_a(lev_min:lev_max),p_ref(lev_min+1:lev_max+1))
    hold on
    plot(rms_u_b(lev_min:lev_max),p_ref(lev_min+1:lev_max+1),'r')
    plot(rms_u_c(lev_min:lev_max),p_ref(lev_min+1:lev_max+1),'g')
    set(gca,'YDir','reverse')
    ylim([p_ref(lev_min+1) p_ref(lev_max+1)])
    
    subplot(1,8,2)
    plot(rms_v_a(lev_min:lev_max),p_ref(lev_min+1:lev_max+1))
    hold on
    plot(rms_v_b(lev_min:lev_max),p_ref(lev_min+1:lev_max+1),'r')
    plot(rms_v_c(lev_min:lev_max),p_ref(lev_min+1:lev_max+1),'g')
    set(gca,'YDir','reverse')
    ylim([p_ref(lev_min+1) p_ref(lev_max+1)])

    subplot(1,8,3)
    plot(rms_t_a(lev_min:lev_max),p_ref(lev_min+1:lev_max+1))
    hold on
    plot(rms_t_b(lev_min:lev_max),p_ref(lev_min+1:lev_max+1),'r')
    plot(rms_t_c(lev_min:lev_max),p_ref(lev_min+1:lev_max+1),'g')
    set(gca,'YDir','reverse')
    ylim([p_ref(lev_min+1) p_ref(lev_max+1)])

    subplot(1,8,4)
    plot(rms_q_a(lev_min:lev_max),p_ref(lev_min+1:lev_max+1))
    hold on
    plot(rms_q_b(lev_min:lev_max),p_ref(lev_min+1:lev_max+1),'r')
    plot(rms_q_c(lev_min:lev_max),p_ref(lev_min+1:lev_max+1),'g')
    set(gca,'YDir','reverse')
    ylim([p_ref(lev_min+1) p_ref(lev_max+1)])

    subplot(1,8,5)
    plot(rms_p_a(lev_min:lev_max),p_ref(lev_min+1:lev_max+1))
    hold on
    plot(rms_p_b(lev_min:lev_max),p_ref(lev_min+1:lev_max+1),'r')
    plot(rms_p_c(lev_min:lev_max),p_ref(lev_min+1:lev_max+1),'g')
    set(gca,'YDir','reverse')
    ylim([p_ref(lev_min+1) p_ref(lev_max+1)])

    subplot(1,8,6)
    plot(rms_qi_a(lev_min:lev_max),p_ref(lev_min+1:lev_max+1))
    hold on
    plot(rms_qi_b(lev_min:lev_max),p_ref(lev_min+1:lev_max+1),'r')
    plot(rms_qi_c(lev_min:lev_max),p_ref(lev_min+1:lev_max+1),'g')
    set(gca,'YDir','reverse')
    ylim([p_ref(lev_min+1) p_ref(lev_max+1)])

    subplot(1,8,7)
    plot(rms_ql_a(lev_min:lev_max),p_ref(lev_min+1:lev_max+1))
    hold on
    plot(rms_ql_b(lev_min:lev_max),p_ref(lev_min+1:lev_max+1),'r')
    plot(rms_ql_c(lev_min:lev_max),p_ref(lev_min+1:lev_max+1),'g')
    set(gca,'YDir','reverse')
    ylim([p_ref(lev_min+1) p_ref(lev_max+1)])

    subplot(1,8,8)
    plot(rms_o3_a(lev_min:lev_max),p_ref(lev_min+1:lev_max+1))
    hold on
    plot(rms_o3_b(lev_min:lev_max),p_ref(lev_min+1:lev_max+1),'r')
    plot(rms_o3_c(lev_min:lev_max),p_ref(lev_min+1:lev_max+1),'g')
    set(gca,'YDir','reverse')
    ylim([p_ref(lev_min+1) p_ref(lev_max+1)])



    figure
    set(gcf,'position',[3 343 1276 576])

    subplot(1,8,1)
    plot(rms_u_a(lev_min:lev_max),lev_min:lev_max)
    hold on
    plot(rms_u_b(lev_min:lev_max),lev_min:lev_max,'r')
    plot(rms_u_c(lev_min:lev_max),lev_min:lev_max,'g')
    set(gca,'YDir','reverse')
    ylim([lev_min lev_max])

    subplot(1,8,2)
    plot(rms_v_a(lev_min:lev_max),lev_min:lev_max)
    hold on
    plot(rms_v_b(lev_min:lev_max),lev_min:lev_max,'r')
    plot(rms_v_c(lev_min:lev_max),lev_min:lev_max,'g')
    set(gca,'YDir','reverse')
    ylim([lev_min lev_max])

    subplot(1,8,3)
    plot(rms_t_a(lev_min:lev_max),lev_min:lev_max)
    hold on
    plot(rms_t_b(lev_min:lev_max),lev_min:lev_max,'r')
    plot(rms_t_c(lev_min:lev_max),lev_min:lev_max,'g')
    set(gca,'YDir','reverse')
    ylim([lev_min lev_max])

    subplot(1,8,4)
    plot(rms_q_a(lev_min:lev_max),lev_min:lev_max)
    hold on
    plot(rms_q_b(lev_min:lev_max),lev_min:lev_max,'r')
    plot(rms_q_c(lev_min:lev_max),lev_min:lev_max,'g')
    set(gca,'YDir','reverse')
    ylim([lev_min lev_max])

    subplot(1,8,5)
    plot(rms_p_a(lev_min:lev_max),lev_min:lev_max)
    hold on
    plot(rms_p_b(lev_min:lev_max),lev_min:lev_max,'r')
    plot(rms_p_c(lev_min:lev_max),lev_min:lev_max,'g')
    set(gca,'YDir','reverse')
    ylim([lev_min lev_max])

    subplot(1,8,6)
    plot(rms_qi_a(lev_min:lev_max),lev_min:lev_max)
    hold on
    plot(rms_qi_b(lev_min:lev_max),lev_min:lev_max,'r')
    plot(rms_qi_c(lev_min:lev_max),lev_min:lev_max,'g')
    set(gca,'YDir','reverse')
    ylim([lev_min lev_max])

    subplot(1,8,7)
    plot(rms_ql_a(lev_min:lev_max),lev_min:lev_max)
    hold on
    plot(rms_ql_b(lev_min:lev_max),lev_min:lev_max,'r')
    plot(rms_ql_c(lev_min:lev_max),lev_min:lev_max,'g')
    set(gca,'YDir','reverse')
    ylim([lev_min lev_max])

    subplot(1,8,8)
    plot(rms_o3_a(lev_min:lev_max),lev_min:lev_max)
    hold on
    plot(rms_o3_b(lev_min:lev_max),lev_min:lev_max,'r')
    plot(rms_o3_c(lev_min:lev_max),lev_min:lev_max,'g')
    set(gca,'YDir','reverse')
    ylim([lev_min lev_max])


end



