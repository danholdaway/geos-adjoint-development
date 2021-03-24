close all
clear
clc

%Choose whether to compute the correlations.
calc_rms = 1;

%Choose Region,
% 1 = Global
% 2 = Tropics 30S to 30N, below 100hPa
% 3 = Northern hemisphere
% 4 = Southern hemisphere
region = 1;

%Choose the time being considered: 16th at 12z would be 6_12
% file_end = '6_0400'; % 1 hours
% file_end = '6_0600'; % 3 hours
% file_end = '6_0900'; % 6 hours
file_end = '6_1200'; % 9 hours
% file_end = '6_1500'; % 12 hours
% file_end = '6_1800'; % 15 hours
% file_end = '6_2100'; % 18 hours
% file_end = '7_0000'; % 21 hours
% file_end = '7_0300'; % 24 hours

cd /home/drholdaw/Lin_Moist_Physics/Inputs/
pref
 
%Load Free (background) State.
cd /discover/nobackup/drholdaw/atmp.22292/prog/prog_free/
file = ['x0011dh_a.prog.eta.2013011',file_end(1:6),'z.nc4'];

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
cd /discover/nobackup/drholdaw/atmp.22292/prog/prog_replay10/
file = ['x0011dh_a.prog.eta.2013011',file_end(1:6),'z.nc4'];

u_replay = ncread(file,'u');
v_replay = ncread(file,'v');
t_replay = ncread(file,'tv');
q_replay = ncread(file,'sphu');
p_replay = ncread(file,'delp');
qi_replay = ncread(file,'qitot');
ql_replay = ncread(file,'qltot');
o3_replay = ncread(file,'ozone');

%Load first TLM state to compare.
cd /discover/nobackup/drholdaw/atmp.22293/sens.20130117.000000%/fvpert_
file = ['x0011dh_a.fvpert.eta.2013011',file_end,'z_BL1_MOI0_GQ0_TF900.nc4'];
% file = 'fvpert.eta.nc4'

u_tlma = ncread(file,'u');
v_tlma = ncread(file,'v');
t_tlma = ncread(file,'tv');
q_tlma = ncread(file,'sphu');
p_tlma = ncread(file,'delp');
qi_tlma = ncread(file,'qitot');
ql_tlma = ncread(file,'qltot');
o3_tlma = ncread(file,'ozone');

%Load second TLM state to compare.
cd /discover/nobackup/drholdaw/atmp.22293/sens.20130117.000000%/fvpert_old/
file = ['x0011dh_a.fvpert.eta.2013011',file_end,'z_BL2_MOI0_GQ0_TF900.nc4'];
% file = 'fvpert.eta.nc4'

u_tlmb = ncread(file,'u');
v_tlmb = ncread(file,'v');
t_tlmb = ncread(file,'tv');
q_tlmb = ncread(file,'sphu');
p_tlmb = ncread(file,'delp');
qi_tlmb = ncread(file,'qitot');
ql_tlmb = ncread(file,'qltot');
o3_tlmb = ncread(file,'ozone');

%Load second TLM state to compare.
cd /discover/nobackup/drholdaw/atmp.22293/sens.20130117.000000%/fvpert_old/
file = ['x0011dh_a.fvpert.eta.2013011',file_end,'z_BL2_MOI2_GQ0_TF900_FILT24.nc4'];
% file = 'fvpert.eta.nc4'

u_tlmc = ncread(file,'u');
v_tlmc = ncread(file,'v');
t_tlmc = ncread(file,'tv');
q_tlmc = ncread(file,'sphu');
p_tlmc = ncread(file,'delp');
qi_tlmc = ncread(file,'qitot');
ql_tlmc = ncread(file,'qltot');
o3_tlmc = ncread(file,'ozone');

cd /home/drholdaw/Lin_Moist_Physics/Bac_paper_figs/

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
    A = zeros(im,jm,lm);
    for k = 1:lm
        A(1,:,k) = (2*pi^2*r_earth^2/(length(lon)*length(lat))) * sind(lat_pole);
        for i = 2:length(lon)
            A(i,:,k) = A(1,:,k); %Make into rank-2 array for simplicity
        end
    end
        
    %A = 1.0;
    
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
    
    
    diff2_u_a = A.*(u_nlm - u_tlma).^2;
    diff2_u_b = A.*(u_nlm - u_tlmb).^2;
    diff2_u_c = A.*(u_nlm - u_tlmc).^2;
    diff2_v_a = A.*(v_nlm - v_tlma).^2;
    diff2_v_b = A.*(v_nlm - v_tlmb).^2;
    diff2_v_c = A.*(v_nlm - v_tlmc).^2;
    diff2_t_a = A.*(t_nlm - t_tlma).^2;
    diff2_t_b = A.*(t_nlm - t_tlmb).^2;
    diff2_t_c = A.*(t_nlm - t_tlmc).^2;
    diff2_q_a = A.*(q_nlm - q_tlma).^2;
    diff2_q_b = A.*(q_nlm - q_tlmb).^2;
    diff2_q_c = A.*(q_nlm - q_tlmc).^2;
    diff2_p_a = A.*(p_nlm - p_tlma).^2;
    diff2_p_b = A.*(p_nlm - p_tlmb).^2;
    diff2_p_c = A.*(p_nlm - p_tlmc).^2;
    diff2_qi_a = A.*(qi_nlm - qi_tlma).^2;
    diff2_qi_b = A.*(qi_nlm - qi_tlmb).^2;
    diff2_qi_c = A.*(qi_nlm - qi_tlmc).^2;
    diff2_ql_a = A.*(ql_nlm - ql_tlma).^2;
    diff2_ql_b = A.*(ql_nlm - ql_tlmb).^2;
    diff2_ql_c = A.*(ql_nlm - ql_tlmc).^2;
    diff2_o3_a = A.*(o3_nlm - o3_tlma).^2;
    diff2_o3_b = A.*(o3_nlm - o3_tlmb).^2;
    diff2_o3_c = A.*(o3_nlm - o3_tlmc).^2;
    
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

        rms_u_a(k) = sqrt(  sum(sum(diff2_u_a(lon_min:lon_max,lat_min:lat_max,k))) / (sumA*n) );
        rms_u_b(k) = sqrt(  sum(sum(diff2_u_b(lon_min:lon_max,lat_min:lat_max,k))) / (sumA*n) );
        rms_u_c(k) = sqrt(  sum(sum(diff2_u_c(lon_min:lon_max,lat_min:lat_max,k))) / (sumA*n) );

        rms_v_a(k) = sqrt(  sum(sum(diff2_v_a(lon_min:lon_max,lat_min:lat_max,k))) / (sumA*n) );
        rms_v_b(k) = sqrt(  sum(sum(diff2_v_b(lon_min:lon_max,lat_min:lat_max,k))) / (sumA*n) );
        rms_v_c(k) = sqrt(  sum(sum(diff2_v_c(lon_min:lon_max,lat_min:lat_max,k))) / (sumA*n) );
        
        rms_t_a(k) = sqrt(  sum(sum(diff2_t_a(lon_min:lon_max,lat_min:lat_max,k))) / (sumA*n) );
        rms_t_b(k) = sqrt(  sum(sum(diff2_t_b(lon_min:lon_max,lat_min:lat_max,k))) / (sumA*n) );
        rms_t_c(k) = sqrt(  sum(sum(diff2_t_c(lon_min:lon_max,lat_min:lat_max,k))) / (sumA*n) );
        
        rms_q_a(k) = sqrt(  sum(sum(diff2_q_a(lon_min:lon_max,lat_min:lat_max,k))) / (sumA*n) );
        rms_q_b(k) = sqrt(  sum(sum(diff2_q_b(lon_min:lon_max,lat_min:lat_max,k))) / (sumA*n) );
        rms_q_c(k) = sqrt(  sum(sum(diff2_q_c(lon_min:lon_max,lat_min:lat_max,k))) / (sumA*n) );
        
        rms_p_a(k) = sqrt(  sum(sum(diff2_p_a(lon_min:lon_max,lat_min:lat_max,k))) / (sumA*n) );
        rms_p_b(k) = sqrt(  sum(sum(diff2_p_b(lon_min:lon_max,lat_min:lat_max,k))) / (sumA*n) );
        rms_p_c(k) = sqrt(  sum(sum(diff2_p_c(lon_min:lon_max,lat_min:lat_max,k))) / (sumA*n) );

        rms_qi_a(k) = sqrt(  sum(sum(diff2_qi_a(lon_min:lon_max,lat_min:lat_max,k))) / (sumA*n) );
        rms_qi_b(k) = sqrt(  sum(sum(diff2_qi_b(lon_min:lon_max,lat_min:lat_max,k))) / (sumA*n) );
        rms_qi_c(k) = sqrt(  sum(sum(diff2_qi_c(lon_min:lon_max,lat_min:lat_max,k))) / (sumA*n) );
        
        rms_ql_a(k) = sqrt(  sum(sum(diff2_ql_a(lon_min:lon_max,lat_min:lat_max,k))) / (sumA*n) );
        rms_ql_b(k) = sqrt(  sum(sum(diff2_ql_b(lon_min:lon_max,lat_min:lat_max,k))) / (sumA*n) );
        rms_ql_c(k) = sqrt(  sum(sum(diff2_ql_c(lon_min:lon_max,lat_min:lat_max,k))) / (sumA*n) );
        
        rms_o3_a(k) = sqrt(  sum(sum(diff2_o3_a(lon_min:lon_max,lat_min:lat_max,k))) / (sumA*n) );
        rms_o3_b(k) = sqrt(  sum(sum(diff2_o3_b(lon_min:lon_max,lat_min:lat_max,k))) / (sumA*n) );
        rms_o3_c(k) = sqrt(  sum(sum(diff2_o3_c(lon_min:lon_max,lat_min:lat_max,k))) / (sumA*n) );
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
    
    grey = 0.75;
    fontsize = 11;
    lin_wid = 1.2;
    line_wid_det = 0.6;

    
    
    figure
    set(gcf,'position',[788   486   491   433])
 
    subplot(1,2,1)
    hold on
    hl11 = line(rms_ql_b(lev_min:lev_max),p_ref(lev_min+1:lev_max+1),'Color','k','LineWidth',lin_wid,'LineStyle','--');
    hl12 = line(rms_ql_c(lev_min:lev_max),p_ref(lev_min+1:lev_max+1),'Color','k','LineWidth',lin_wid);
           
    set(gca,'Ydir','reverse')
    ylim([50 1000])
    legend('Dry Physics','Moist Physics')

    set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
    xlabel('(a) q_l\prime (kgkg^{-1})','FontSize',fontsize,'FontName','TimesNewRoman')
    ylabel('Pressure (hPa)','FontSize',fontsize,'FontName','TimesNewRoman')
    box on

    set(gca,'XTick',[0 0.5e-7 1e-7])
    set(gca,'XTickLabel',{'0' '0.5' '1e-7'})
    
    pos = get(gca,'position');
    set(gca,'position',[pos(1) pos(2) pos(3) pos(4)])
    
    ax1 = gca;
    ax3 = axes('Position',get(ax1,'Position'),...
           'XAxisLocation','top',...
           'YAxisLocation','right',...
           'YTick',[],...
           'XTick',[-8e-9 -6e-9 -4e-9 -2e-9 0],...
           'XTickLabel',{'-8e-9' '-6' '-4' '-2' '0'},...
           'YDir','reverse',...
           'ylim',[50 1000],...
           'Color','none',...
           'XColor','k','YColor','k');
       
    subplot(1,2,2)
    hold on
    hl21 = line(rms_qi_b(lev_min:lev_max),p_ref(lev_min+1:lev_max+1),'Color','k','LineWidth',lin_wid,'LineStyle','--');
    hl22 = line(rms_qi_c(lev_min:lev_max),p_ref(lev_min+1:lev_max+1),'Color','k','LineWidth',lin_wid);
    
    set(gca,'Ydir','reverse')
    set(gca,'YAxisLocation','right')
    ylim([50 1000])

    set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
    xlabel('(b) q_i\prime (kgkg^{-1})','FontSize',fontsize,'FontName','TimesNewRoman')
    ylabel('Pressure (hPa)','FontSize',fontsize,'FontName','TimesNewRoman')
    box on

    pos = get(gca,'position');
    set(gca,'position',[0.95*pos(1) pos(2) pos(3) pos(4)])
    
    set(gca,'XTick',[0 4e-8 8e-8])
    set(gca,'XTickLabel',{'0' '4' '8e-8'})
    
    ax2 = gca;
    ax4 = axes('Position',get(ax2,'Position'),...
           'XAxisLocation','top',...
           'YAxisLocation','right',...
           'YTick',[],...
           'XTick',[-2e-8 -1.5e-8 -1e-8 -0.5e-8 0],...
           'XTickLabel',{'-2e-8' '-1.5' '-1' '-0.5' '0'},...
           'YDir','reverse',...
           'ylim',[50 1000],...
           'Color','none',...
           'XColor','k','YColor','k');

    
       
    
    hl2 = line(rms_ql_c(lev_min:lev_max)-rms_ql_b(lev_min:lev_max),p_ref(lev_min+1:lev_max+1),'Color',[grey grey grey],'LineWidth',lin_wid,'Parent',ax3);
    hl3 = line(rms_qi_c(lev_min:lev_max)-rms_qi_b(lev_min:lev_max),p_ref(lev_min+1:lev_max+1),'Color',[grey grey grey],'LineWidth',lin_wid,'Parent',ax4);
    h = legend('Moist-Dry');
    set(h,'color','w')

    
end



