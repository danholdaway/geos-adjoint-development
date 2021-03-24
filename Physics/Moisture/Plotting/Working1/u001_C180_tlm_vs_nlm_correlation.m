close all
clear
clc

%Choose day to examine
date_start = datenum(2013, 01, 16, 00, 00, 00);

%12th

%Choose number of days to examine, if more than 1 then average plotted.
num_days = 1;

%Vertical resolution
lm = 72;


for i = 0:num_days-1
    
    %Choose lead time, e.g. 24, 36 or 48 hours
    leadtime = 24;

    %Choose Region,
    % 1 = Global
    % 2 = Tropics 23S to 23N, below 100hPa
    % 3 = Northern hemisphere
    % 4 = Southern hemisphere
    region = 1;

    %Get required dates
    daten = date_start + i;
        
    fprintf('Date is %s \n', datestr(daten, 'yyyymmddHHMM'))
    
    dater = datestr(daten - 0.1250, 'yyyymmddHHMM');
    date00 = datestr(daten + 0, 'yyyymmddHHMM');
    dateplot = datestr(daten + leadtime/24, 'yyyymmddHHMM');

    cd /home/drholdaw/Lin_Moist_Physics/Inputs/
    pref

    %Load Free (background) State.
    dir = ['/discover/nobackup/drholdaw/u001_C180/prog/Y',dater(1:4),'/M',dater(5:6),'/D',dater(7:8),'/H21-NLM-free/']
    cd(dir)
    file = ['u001_C180.prog.eta.',dater(1:8),'_21z+',dateplot(1:8),'_',dateplot(9:10),'z.nc4']

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

    %Load Perturbed (analysis) state
    dir = ['/discover/nobackup/drholdaw/u001_C180/prog/Y',dater(1:4),'/M',dater(5:6),'/D',dater(7:8),'/H21-NLM-replay/']
    cd(dir)
    file = ['u001_C180.prog.eta.',dater(1:8),'_21z+',dateplot(1:8),'_',dateplot(9:10),'z.nc4']

    u_replay = ncread(file,'u');
    v_replay = ncread(file,'v');
    t_replay = ncread(file,'tv');
    q_replay = ncread(file,'sphu');
    p_replay = ncread(file,'delp');
    qi_replay = ncread(file,'qitot');
    ql_replay = ncread(file,'qltot');
    o3_replay = ncread(file,'ozone');

    %Load TLM state to compare.
    dir = ['/discover/nobackup/drholdaw/u001_C180/prog/Y',dater(1:4),'/M',dater(5:6),'/D',dater(7:8),'/H21-TLM-GQ1-PH1/']
    cd(dir)
    file = ['u001_C180.fvpert.eta.',dater(1:8),'_21z+',dateplot(1:8),'_',dateplot(9:10),'z.nc4']

    u_tl1 = ncread(file,'u');
    v_tl1 = ncread(file,'v');
    t_tl1 = ncread(file,'tv');
    q_tl1 = ncread(file,'sphu');
    p_tl1 = ncread(file,'delp');
    qi_tl1 = ncread(file,'qitot');
    ql_tl1 = ncread(file,'qltot');
    o3_tl1 = ncread(file,'ozone');

    %Load TLM state to compare.
    dir = ['/discover/nobackup/drholdaw/tmp.22365/sens.20130118.000000'];
    cd(dir)
    file = ['u001_C180.fvpert.eta.20130117_0000z.nc4'];

    u_tl2 = ncread(file,'u');
    v_tl2 = ncread(file,'v');
    t_tl2 = ncread(file,'tv');
    q_tl2 = ncread(file,'sphu');
    p_tl2 = ncread(file,'delp');
    qi_tl2 = ncread(file,'qitot');
    ql_tl2 = ncread(file,'qltot');
    o3_tl2 = ncread(file,'ozone');

%     %Load TLM state to compare.
%     dir = ['/discover/nobackup/drholdaw/u001_C180/prog/Y',dater(1:4),'/M',dater(5:6),'/D',dater(7:8),'/H21-TLM-GQ1-PH2/'];
%     cd(dir)
%     file = ['u001_C180.fvpert.eta.',dater(1:8),'_21z+',dateplot(1:8),'_',dateplot(9:10),'z.nc4'];
% 
%     u_tl3 = ncread(file,'u');
%     v_tl3 = ncread(file,'v');
%     t_tl3 = ncread(file,'tv');
%     q_tl3 = ncread(file,'sphu');
%     p_tl3 = ncread(file,'delp');
%     qi_tl3 = ncread(file,'qitot');
%     ql_tl3 = ncread(file,'qltot');
%     o3_tl3 = ncread(file,'ozone');
% 
%     %Load TLM state to compare.
%     dir = ['/discover/nobackup/drholdaw/u001_C180/prog/Y',dater(1:4),'/M',dater(5:6),'/D',dater(7:8),'/H21-TLM-GQ2-PH0/'];
%     cd(dir)
%     file = ['u001_C180.fvpert.eta.',date00(1:8),'_00z+',dateplot(1:8),'_',dateplot(9:10),'z.nc4'];
% 
%     u_tl4 = ncread(file,'u');
%     v_tl4 = ncread(file,'v');
%     t_tl4 = ncread(file,'tv');
%     q_tl4 = ncread(file,'sphu');
%     p_tl4 = ncread(file,'delp');
%     qi_tl4 = ncread(file,'qitot');
%     ql_tl4 = ncread(file,'qltot');
%     o3_tl4 = ncread(file,'ozone');
% 
%     %Load TLM state to compare.
%     dir = ['/discover/nobackup/drholdaw/u001_C180/prog/Y',dater(1:4),'/M',dater(5:6),'/D',dater(7:8),'/H21-TLM-GQ2-PH1/'];
%     cd(dir)
%     file = ['u001_C180.fvpert.eta.',date00(1:8),'_00z+',dateplot(1:8),'_',dateplot(9:10),'z.nc4'];
% 
%     u_tl5 = ncread(file,'u');
%     v_tl5 = ncread(file,'v');
%     t_tl5 = ncread(file,'tv');
%     q_tl5 = ncread(file,'sphu');
%     p_tl5 = ncread(file,'delp');
%     qi_tl5 = ncread(file,'qitot');
%     ql_tl5 = ncread(file,'qltot');
%     o3_tl5 = ncread(file,'ozone');
% 
%     %Load TLM state to compare.
%     dir = ['/discover/nobackup/drholdaw/u001_C180/prog/Y',dater(1:4),'/M',dater(5:6),'/D',dater(7:8),'/H21-TLM-GQ2-PH2/'];
%     cd(dir)
%     file = ['u001_C180.fvpert.eta.',date00(1:8),'_00z+',dateplot(1:8),'_',dateplot(9:10),'z.nc4'];
% 
%     u_tl6 = ncread(file,'u');
%     v_tl6 = ncread(file,'v');
%     t_tl6 = ncread(file,'tv');
%     q_tl6 = ncread(file,'sphu');
%     p_tl6 = ncread(file,'delp');
%     qi_tl6 = ncread(file,'qitot');
%     ql_tl6 = ncread(file,'qltot');
%     o3_tl6 = ncread(file,'ozone');
% 
%     cd /home/drholdaw/Lin_Moist_Physics/DAS_plotting/Cloud' Fraction'/
% 
    %Compute NL perturbation trajectory.
    u_nlm = u_replay - u_free;
    v_nlm = v_replay - v_free;
    t_nlm = t_replay - t_free;
    q_nlm = q_replay - q_free;
    p_nlm = p_replay - p_free;
    qi_nlm = qi_replay - qi_free;
    ql_nlm = ql_replay - ql_free;
    o3_nlm = o3_replay - o3_free;


    im = length(lon);
    jm = length(lat);
    lm = length(lev);

    fprintf('>Beginning correlation calculation \n')

    %Weighting 
    r_earth = 6378.1;
    lat_pole = lat+lat(1)*-1;
    A = zeros(im,jm,lm);
    for k = 1:lm
        A(1,:,k) = (2*pi^2*r_earth^2/(length(lon)*length(lat))) * sind(lat_pole);
        for j = 2:length(lon)
            A(j,:,k) = A(1,:,k); %Make into rank-2 array for simplicity
        end
    end

    dA = sum(sum(A(:,:,1)));


    if region == 1
        LonMin = 1;
        LonMax = im;
        LatMin = 1;
        LatMax = jm;
        LevMin = 30;
        LevMax = lm;
    elseif region == 2
        LonMin = 1;
        LonMax = im;
        LatMin = 135;
        LatMax = 227;
        LevMin = 30;
        LevMax = lm;
    elseif region == 3
        LonMin = 1;
        LonMax = im;
        LatMin = ceil(jm/2);
        LatMax = jm;
        LevMin = 30;
        LevMax = lm;
    elseif region == 4
        LonMin = 1;
        LonMax = im;
        LatMin = 1;
        LatMax = floor(jm/2);
        LevMin = 30;
        LevMax = lm;    
    end

    %Variance, covariance and correlation
    %Second index is number of variables.
    var_nlm = zeros(lm,8);
    var_tl1 = zeros(lm,8);
    var_tl2 = zeros(lm,8);
    var_tl3 = zeros(lm,8);
    %var_tl4 = zeros(lm,8);
    %var_tl5 = zeros(lm,8);
    %var_tl6 = zeros(lm,8);

    cov_tl1 = zeros(lm,8);
    cov_tl2 = zeros(lm,8);
    %cov_tl3 = zeros(lm,8);
    %cov_tl4 = zeros(lm,8);
    %cov_tl5 = zeros(lm,8);
    %cov_tl6 = zeros(lm,8);
        
    %Initialize correlations
    corr_tl1 = zeros(lm,8);
    corr_tl2 = zeros(lm,8);
    %corr_tl3 = zeros(lm,8);
    %corr_tl4 = zeros(lm,8);
    %corr_tl5 = zeros(lm,8);
    %corr_tl6 = zeros(lm,8);
    
    % U winds
    var_nlm(LevMin:LevMax,1) = sum(sum(u_nlm(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).^2.*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;
    var_tl1(LevMin:LevMax,1) = sum(sum(u_tl1(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).^2.*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;
    var_tl2(LevMin:LevMax,1) = sum(sum(u_tl2(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).^2.*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;
    %var_tl3(LevMin:LevMax,1) = sum(sum(u_tl3(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).^2.*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;
    %var_tl4(LevMin:LevMax,1) = sum(sum(u_tl4(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).^2.*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;
    %var_tl5(LevMin:LevMax,1) = sum(sum(u_tl5(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).^2.*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;
    %var_tl6(LevMin:LevMax,1) = sum(sum(u_tl6(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).^2.*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;

    cov_tl1(LevMin:LevMax,1) = sum(sum(u_nlm(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).*u_tl1(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;
    cov_tl2(LevMin:LevMax,1) = sum(sum(u_nlm(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).*u_tl2(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;
    %cov_tl3(LevMin:LevMax,1) = sum(sum(u_nlm(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).*u_tl3(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;
    %cov_tl4(LevMin:LevMax,1) = sum(sum(u_nlm(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).*u_tl4(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;
    %cov_tl5(LevMin:LevMax,1) = sum(sum(u_nlm(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).*u_tl5(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;
    %cov_tl6(LevMin:LevMax,1) = sum(sum(u_nlm(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).*u_tl6(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;

    corr_tl1(LevMin:LevMax,1) = cov_tl1(LevMin:LevMax,1)./sqrt(var_nlm(LevMin:LevMax,1).*var_tl1(LevMin:LevMax,1));
    corr_tl2(LevMin:LevMax,1) = cov_tl2(LevMin:LevMax,1)./sqrt(var_nlm(LevMin:LevMax,1).*var_tl2(LevMin:LevMax,1));
    %corr_tl3(LevMin:LevMax,1) = cov_tl3(LevMin:LevMax,1)./sqrt(var_nlm(LevMin:LevMax,1).*%var_tl3(LevMin:LevMax,1));
    %corr_tl4(LevMin:LevMax,1) = cov_tl4(LevMin:LevMax,1)./sqrt(var_nlm(LevMin:LevMax,1).*%var_tl4(LevMin:LevMax,1));
    %corr_tl5(LevMin:LevMax,1) = cov_tl5(LevMin:LevMax,1)./sqrt(var_nlm(LevMin:LevMax,1).*%var_tl5(LevMin:LevMax,1));
    %corr_tl6(LevMin:LevMax,1) = cov_tl6(LevMin:LevMax,1)./sqrt(var_nlm(LevMin:LevMax,1).*%var_tl6(LevMin:LevMax,1));

    % V winds
    var_nlm(LevMin:LevMax,2) = sum(sum(v_nlm(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).^2.*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;
    var_tl1(LevMin:LevMax,2) = sum(sum(v_tl1(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).^2.*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;
    var_tl2(LevMin:LevMax,2) = sum(sum(v_tl2(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).^2.*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;
    %var_tl3(LevMin:LevMax,2) = sum(sum(v_tl3(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).^2.*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;
    %var_tl4(LevMin:LevMax,2) = sum(sum(v_tl4(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).^2.*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;
    %var_tl5(LevMin:LevMax,2) = sum(sum(v_tl5(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).^2.*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;
    %var_tl6(LevMin:LevMax,2) = sum(sum(v_tl6(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).^2.*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;

    cov_tl1(LevMin:LevMax,2) = sum(sum(v_nlm(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).*v_tl1(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;
    cov_tl2(LevMin:LevMax,2) = sum(sum(v_nlm(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).*v_tl2(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;
    %cov_tl3(LevMin:LevMax,2) = sum(sum(v_nlm(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).*v_tl3(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;
    %cov_tl4(LevMin:LevMax,2) = sum(sum(v_nlm(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).*v_tl4(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;
    %cov_tl5(LevMin:LevMax,2) = sum(sum(v_nlm(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).*v_tl5(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;
    %cov_tl6(LevMin:LevMax,2) = sum(sum(v_nlm(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).*v_tl6(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;

    corr_tl1(LevMin:LevMax,2) = cov_tl1(LevMin:LevMax,2)./sqrt(var_nlm(LevMin:LevMax,2).*var_tl1(LevMin:LevMax,2));
    corr_tl2(LevMin:LevMax,2) = cov_tl2(LevMin:LevMax,2)./sqrt(var_nlm(LevMin:LevMax,2).*var_tl2(LevMin:LevMax,2));
    %corr_tl3(LevMin:LevMax,2) = cov_tl3(LevMin:LevMax,2)./sqrt(var_nlm(LevMin:LevMax,2).*%var_tl3(LevMin:LevMax,2));
    %corr_tl4(LevMin:LevMax,2) = cov_tl4(LevMin:LevMax,2)./sqrt(var_nlm(LevMin:LevMax,2).*%var_tl4(LevMin:LevMax,2));
    %corr_tl5(LevMin:LevMax,2) = cov_tl5(LevMin:LevMax,2)./sqrt(var_nlm(LevMin:LevMax,2).*%var_tl5(LevMin:LevMax,2));
    %corr_tl6(LevMin:LevMax,2) = cov_tl6(LevMin:LevMax,2)./sqrt(var_nlm(LevMin:LevMax,2).*%var_tl6(LevMin:LevMax,2));

    % Temperature
    var_nlm(LevMin:LevMax,3) = sum(sum(t_nlm(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).^2.*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;
    var_tl1(LevMin:LevMax,3) = sum(sum(t_tl1(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).^2.*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;
    var_tl2(LevMin:LevMax,3) = sum(sum(t_tl2(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).^2.*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;
    %var_tl3(LevMin:LevMax,3) = sum(sum(t_tl3(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).^2.*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;
    %var_tl4(LevMin:LevMax,3) = sum(sum(t_tl4(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).^2.*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;
    %var_tl5(LevMin:LevMax,3) = sum(sum(t_tl5(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).^2.*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;
    %var_tl6(LevMin:LevMax,3) = sum(sum(t_tl6(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).^2.*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;

    cov_tl1(LevMin:LevMax,3) = sum(sum(t_nlm(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).*t_tl1(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;
    cov_tl2(LevMin:LevMax,3) = sum(sum(t_nlm(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).*t_tl2(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;
    %cov_tl3(LevMin:LevMax,3) = sum(sum(t_nlm(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).*t_tl3(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;
    %cov_tl4(LevMin:LevMax,3) = sum(sum(t_nlm(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).*t_tl4(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;
    %cov_tl5(LevMin:LevMax,3) = sum(sum(t_nlm(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).*t_tl5(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;
    %cov_tl6(LevMin:LevMax,3) = sum(sum(t_nlm(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).*t_tl6(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;

    corr_tl1(LevMin:LevMax,3) = cov_tl1(LevMin:LevMax,3)./sqrt(var_nlm(LevMin:LevMax,3).*var_tl1(LevMin:LevMax,3));
    corr_tl2(LevMin:LevMax,3) = cov_tl2(LevMin:LevMax,3)./sqrt(var_nlm(LevMin:LevMax,3).*var_tl2(LevMin:LevMax,3));
    %corr_tl3(LevMin:LevMax,3) = cov_tl3(LevMin:LevMax,3)./sqrt(var_nlm(LevMin:LevMax,3).*%var_tl3(LevMin:LevMax,3));
    %corr_tl4(LevMin:LevMax,3) = cov_tl4(LevMin:LevMax,3)./sqrt(var_nlm(LevMin:LevMax,3).*%var_tl4(LevMin:LevMax,3));
    %corr_tl5(LevMin:LevMax,3) = cov_tl5(LevMin:LevMax,3)./sqrt(var_nlm(LevMin:LevMax,3).*%var_tl5(LevMin:LevMax,3));
    %corr_tl6(LevMin:LevMax,3) = cov_tl6(LevMin:LevMax,3)./sqrt(var_nlm(LevMin:LevMax,3).*%var_tl6(LevMin:LevMax,3));

    % Specific Humidity
    var_nlm(LevMin:LevMax,4) = sum(sum(q_nlm(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).^2.*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;
    var_tl1(LevMin:LevMax,4) = sum(sum(q_tl1(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).^2.*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;
    var_tl2(LevMin:LevMax,4) = sum(sum(q_tl2(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).^2.*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;
    %var_tl3(LevMin:LevMax,4) = sum(sum(q_tl3(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).^2.*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;
    %var_tl4(LevMin:LevMax,4) = sum(sum(q_tl4(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).^2.*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;
    %var_tl5(LevMin:LevMax,4) = sum(sum(q_tl5(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).^2.*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;
    %var_tl6(LevMin:LevMax,4) = sum(sum(q_tl6(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).^2.*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;

    cov_tl1(LevMin:LevMax,4) = sum(sum(q_nlm(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).*q_tl1(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;
    cov_tl2(LevMin:LevMax,4) = sum(sum(q_nlm(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).*q_tl2(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;
    %cov_tl3(LevMin:LevMax,4) = sum(sum(q_nlm(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).*q_tl3(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;
    %cov_tl4(LevMin:LevMax,4) = sum(sum(q_nlm(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).*q_tl4(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;
    %cov_tl5(LevMin:LevMax,4) = sum(sum(q_nlm(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).*q_tl5(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;
    %cov_tl6(LevMin:LevMax,4) = sum(sum(q_nlm(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).*q_tl6(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;

    corr_tl1(LevMin:LevMax,4) = cov_tl1(LevMin:LevMax,4)./sqrt(var_nlm(LevMin:LevMax,4).*var_tl1(LevMin:LevMax,4));
    corr_tl2(LevMin:LevMax,4) = cov_tl2(LevMin:LevMax,4)./sqrt(var_nlm(LevMin:LevMax,4).*var_tl2(LevMin:LevMax,4));
    %corr_tl3(LevMin:LevMax,4) = cov_tl3(LevMin:LevMax,4)./sqrt(var_nlm(LevMin:LevMax,4).*%var_tl3(LevMin:LevMax,4));
    %corr_tl4(LevMin:LevMax,4) = cov_tl4(LevMin:LevMax,4)./sqrt(var_nlm(LevMin:LevMax,4).*%var_tl4(LevMin:LevMax,4));
    %corr_tl5(LevMin:LevMax,4) = cov_tl5(LevMin:LevMax,4)./sqrt(var_nlm(LevMin:LevMax,4).*%var_tl5(LevMin:LevMax,4));
    %corr_tl6(LevMin:LevMax,4) = cov_tl6(LevMin:LevMax,4)./sqrt(var_nlm(LevMin:LevMax,4).*%var_tl6(LevMin:LevMax,4));

    % Pressure
    var_nlm(LevMin:LevMax,5) = sum(sum(p_nlm(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).^2.*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;
    var_tl1(LevMin:LevMax,5) = sum(sum(p_tl1(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).^2.*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;
    var_tl2(LevMin:LevMax,5) = sum(sum(p_tl2(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).^2.*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;
    %var_tl3(LevMin:LevMax,5) = sum(sum(p_tl3(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).^2.*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;
    %var_tl4(LevMin:LevMax,5) = sum(sum(p_tl4(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).^2.*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;
    %var_tl5(LevMin:LevMax,5) = sum(sum(p_tl5(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).^2.*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;
    %var_tl6(LevMin:LevMax,5) = sum(sum(p_tl6(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).^2.*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;

    cov_tl1(LevMin:LevMax,5) = sum(sum(p_nlm(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).*p_tl1(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;
    cov_tl2(LevMin:LevMax,5) = sum(sum(p_nlm(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).*p_tl2(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;
    %cov_tl3(LevMin:LevMax,5) = sum(sum(p_nlm(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).*p_tl3(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;
    %cov_tl4(LevMin:LevMax,5) = sum(sum(p_nlm(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).*p_tl4(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;
    %cov_tl5(LevMin:LevMax,5) = sum(sum(p_nlm(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).*p_tl5(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;
    %cov_tl6(LevMin:LevMax,5) = sum(sum(p_nlm(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).*p_tl6(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;

    corr_tl1(LevMin:LevMax,5) = cov_tl1(LevMin:LevMax,5)./sqrt(var_nlm(LevMin:LevMax,5).*var_tl1(LevMin:LevMax,5));
    corr_tl2(LevMin:LevMax,5) = cov_tl2(LevMin:LevMax,5)./sqrt(var_nlm(LevMin:LevMax,5).*var_tl2(LevMin:LevMax,5));
    %corr_tl3(LevMin:LevMax,5) = cov_tl3(LevMin:LevMax,5)./sqrt(var_nlm(LevMin:LevMax,5).*%var_tl3(LevMin:LevMax,5));
    %corr_tl4(LevMin:LevMax,5) = cov_tl4(LevMin:LevMax,5)./sqrt(var_nlm(LevMin:LevMax,5).*%var_tl4(LevMin:LevMax,5));
    %corr_tl5(LevMin:LevMax,5) = cov_tl5(LevMin:LevMax,5)./sqrt(var_nlm(LevMin:LevMax,5).*%var_tl5(LevMin:LevMax,5));
    %corr_tl6(LevMin:LevMax,5) = cov_tl6(LevMin:LevMax,5)./sqrt(var_nlm(LevMin:LevMax,5).*%var_tl6(LevMin:LevMax,5));

    % Cloud Liquid Water
    var_nlm(LevMin:LevMax,6) = sum(sum(ql_nlm(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).^2.*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;
    var_tl1(LevMin:LevMax,6) = sum(sum(ql_tl1(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).^2.*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;
    var_tl2(LevMin:LevMax,6) = sum(sum(ql_tl2(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).^2.*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;
    %var_tl3(LevMin:LevMax,6) = sum(sum(ql_tl3(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).^2.*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;
    %var_tl4(LevMin:LevMax,6) = sum(sum(ql_tl4(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).^2.*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;
    %var_tl5(LevMin:LevMax,6) = sum(sum(ql_tl5(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).^2.*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;
    %var_tl6(LevMin:LevMax,6) = sum(sum(ql_tl6(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).^2.*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;

    cov_tl1(LevMin:LevMax,6) = sum(sum(ql_nlm(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).*ql_tl1(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;
    cov_tl2(LevMin:LevMax,6) = sum(sum(ql_nlm(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).*ql_tl2(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;
    %cov_tl3(LevMin:LevMax,6) = sum(sum(ql_nlm(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).*ql_tl3(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;
    %cov_tl4(LevMin:LevMax,6) = sum(sum(ql_nlm(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).*ql_tl4(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;
    %cov_tl5(LevMin:LevMax,6) = sum(sum(ql_nlm(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).*ql_tl5(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;
    %cov_tl6(LevMin:LevMax,6) = sum(sum(ql_nlm(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).*ql_tl6(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;

    corr_tl1(LevMin:LevMax,6) = cov_tl1(LevMin:LevMax,6)./sqrt(var_nlm(LevMin:LevMax,6).*var_tl1(LevMin:LevMax,6));
    corr_tl2(LevMin:LevMax,6) = cov_tl2(LevMin:LevMax,6)./sqrt(var_nlm(LevMin:LevMax,6).*var_tl2(LevMin:LevMax,6));
    %corr_tl3(LevMin:LevMax,6) = cov_tl3(LevMin:LevMax,6)./sqrt(var_nlm(LevMin:LevMax,6).*%var_tl3(LevMin:LevMax,6));
    %corr_tl4(LevMin:LevMax,6) = cov_tl4(LevMin:LevMax,6)./sqrt(var_nlm(LevMin:LevMax,6).*%var_tl4(LevMin:LevMax,6));
    %corr_tl5(LevMin:LevMax,6) = cov_tl5(LevMin:LevMax,6)./sqrt(var_nlm(LevMin:LevMax,6).*%var_tl5(LevMin:LevMax,6));
    %corr_tl6(LevMin:LevMax,6) = cov_tl6(LevMin:LevMax,6)./sqrt(var_nlm(LevMin:LevMax,6).*%var_tl6(LevMin:LevMax,6));

    % Cloud Liquid Ice
    var_nlm(LevMin:LevMax,7) = sum(sum(qi_nlm(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).^2.*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;
    var_tl1(LevMin:LevMax,7) = sum(sum(qi_tl1(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).^2.*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;
    var_tl2(LevMin:LevMax,7) = sum(sum(qi_tl2(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).^2.*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;
    %var_tl3(LevMin:LevMax,7) = sum(sum(qi_tl3(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).^2.*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;
    %var_tl4(LevMin:LevMax,7) = sum(sum(qi_tl4(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).^2.*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;
    %var_tl5(LevMin:LevMax,7) = sum(sum(qi_tl5(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).^2.*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;
    %var_tl6(LevMin:LevMax,7) = sum(sum(qi_tl6(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).^2.*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;

    cov_tl1(LevMin:LevMax,7) = sum(sum(qi_nlm(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).*qi_tl1(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;
    cov_tl2(LevMin:LevMax,7) = sum(sum(qi_nlm(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).*qi_tl2(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;
    %cov_tl3(LevMin:LevMax,7) = sum(sum(qi_nlm(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).*qi_tl3(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;
    %cov_tl4(LevMin:LevMax,7) = sum(sum(qi_nlm(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).*qi_tl4(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;
    %cov_tl5(LevMin:LevMax,7) = sum(sum(qi_nlm(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).*qi_tl5(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;
    %cov_tl6(LevMin:LevMax,7) = sum(sum(qi_nlm(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).*qi_tl6(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;

    corr_tl1(LevMin:LevMax,7) = cov_tl1(LevMin:LevMax,7)./sqrt(var_nlm(LevMin:LevMax,7).*var_tl1(LevMin:LevMax,7));
    corr_tl2(LevMin:LevMax,7) = cov_tl2(LevMin:LevMax,7)./sqrt(var_nlm(LevMin:LevMax,7).*var_tl2(LevMin:LevMax,7));
    %corr_tl3(LevMin:LevMax,7) = cov_tl3(LevMin:LevMax,7)./sqrt(var_nlm(LevMin:LevMax,7).*%var_tl3(LevMin:LevMax,7));
    %corr_tl4(LevMin:LevMax,7) = cov_tl4(LevMin:LevMax,7)./sqrt(var_nlm(LevMin:LevMax,7).*%var_tl4(LevMin:LevMax,7));
    %corr_tl5(LevMin:LevMax,7) = cov_tl5(LevMin:LevMax,7)./sqrt(var_nlm(LevMin:LevMax,7).*%var_tl5(LevMin:LevMax,7));
    %corr_tl6(LevMin:LevMax,7) = cov_tl6(LevMin:LevMax,7)./sqrt(var_nlm(LevMin:LevMax,7).*%var_tl6(LevMin:LevMax,7));

    % Ozone
    var_nlm(LevMin:LevMax,8) = sum(sum(o3_nlm(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).^2.*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;
    var_tl1(LevMin:LevMax,8) = sum(sum(o3_tl1(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).^2.*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;
    var_tl2(LevMin:LevMax,8) = sum(sum(o3_tl2(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).^2.*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;
    %var_tl3(LevMin:LevMax,8) = sum(sum(o3_tl3(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).^2.*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;
    %var_tl4(LevMin:LevMax,8) = sum(sum(o3_tl4(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).^2.*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;
    %var_tl5(LevMin:LevMax,8) = sum(sum(o3_tl5(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).^2.*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;
    %var_tl6(LevMin:LevMax,8) = sum(sum(o3_tl6(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).^2.*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;

    cov_tl1(LevMin:LevMax,8) = sum(sum(o3_nlm(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).*o3_tl1(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;
    cov_tl2(LevMin:LevMax,8) = sum(sum(o3_nlm(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).*o3_tl2(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;
    %cov_tl3(LevMin:LevMax,8) = sum(sum(o3_nlm(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).*o3_tl3(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;
    %cov_tl4(LevMin:LevMax,8) = sum(sum(o3_nlm(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).*o3_tl4(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;
    %cov_tl5(LevMin:LevMax,8) = sum(sum(o3_nlm(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).*o3_tl5(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;
    %cov_tl6(LevMin:LevMax,8) = sum(sum(o3_nlm(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).*o3_tl6(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax).*A(LonMin:LonMax,LatMin:LatMax,LevMin:LevMax))) / dA;

    corr_tl1(LevMin:LevMax,8) = cov_tl1(LevMin:LevMax,8)./sqrt(var_nlm(LevMin:LevMax,8).*var_tl1(LevMin:LevMax,8));
    corr_tl2(LevMin:LevMax,8) = cov_tl2(LevMin:LevMax,8)./sqrt(var_nlm(LevMin:LevMax,8).*var_tl2(LevMin:LevMax,8));
    %corr_tl3(LevMin:LevMax,8) = cov_tl3(LevMin:LevMax,8)./sqrt(var_nlm(LevMin:LevMax,8).*%var_tl3(LevMin:LevMax,8));
    %corr_tl4(LevMin:LevMax,8) = cov_tl4(LevMin:LevMax,8)./sqrt(var_nlm(LevMin:LevMax,8).*%var_tl4(LevMin:LevMax,8));
    %corr_tl5(LevMin:LevMax,8) = cov_tl5(LevMin:LevMax,8)./sqrt(var_nlm(LevMin:LevMax,8).*%%var_tl5(LevMin:LevMax,8));
    %corr_tl6(LevMin:LevMax,8) = cov_tl6(LevMin:LevMax,8)./sqrt(var_nlm(LevMin:LevMax,8).*%var_tl6(LevMin:LevMax,8));

    fprintf('  >Done correlation calculation \n\n\n')
    
end
    

lin_wid = 1.75;
fontsize = 11;
fontsize1 = 12;


figure
set(gcf,'position',[3 343 1276 576])

subplot(1,8,1)
plot(corr_tl1(LevMin:LevMax,1),p_ref(LevMin:LevMax),'k','LineWidth',lin_wid)
hold on
plot(corr_tl2(LevMin:LevMax,1),p_ref(LevMin:LevMax),'r','LineWidth',lin_wid)
% plot(corr_tl3(LevMin:LevMax,1),p_ref(LevMin:LevMax),'b--','LineWidth',lin_wid)
set(gca,'YDir','reverse')
xlim([0 1])
ylim([p_ref(LevMin+1) p_ref(LevMax)])
title('u (ms^{-1})','FontSize',fontsize1,'FontName','TimesNewRoman')
ylabel('Height (hPa)','FontSize',fontsize1,'FontName','TimesNewRoman')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

subplot(1,8,2)
plot(corr_tl1(LevMin:LevMax,2),p_ref(LevMin:LevMax),'k','LineWidth',lin_wid)
hold on
plot(corr_tl2(LevMin:LevMax,2),p_ref(LevMin:LevMax),'r','LineWidth',lin_wid)
% plot(corr_tl3(LevMin:LevMax,2),p_ref(LevMin:LevMax),'b--','LineWidth',lin_wid)
set(gca,'YDir','reverse')
xlim([0 1])
ylim([p_ref(LevMin+1) p_ref(LevMax)])
title('v (ms^{-1})','FontSize',fontsize1,'FontName','TimesNewRoman')
set(gca,'YTickLabel',[])
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

subplot(1,8,3)
plot(corr_tl1(LevMin:LevMax,3),p_ref(LevMin:LevMax),'k','LineWidth',lin_wid)
hold on
plot(corr_tl2(LevMin:LevMax,3),p_ref(LevMin:LevMax),'r','LineWidth',lin_wid)
% plot(corr_tl3(LevMin:LevMax,3),p_ref(LevMin:LevMax),'b--','LineWidth',lin_wid)
set(gca,'YDir','reverse')
xlim([0 1])
ylim([p_ref(LevMin+1) p_ref(LevMax)])
title('T_v (K)','FontSize',fontsize1,'FontName','TimesNewRoman')
set(gca,'YTickLabel',[])
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

subplot(1,8,4)
plot(corr_tl1(LevMin:LevMax,4),p_ref(LevMin:LevMax),'k','LineWidth',lin_wid)
hold on
plot(corr_tl2(LevMin:LevMax,4),p_ref(LevMin:LevMax),'r','LineWidth',lin_wid)
% plot(corr_tl3(LevMin:LevMax,4),p_ref(LevMin:LevMax),'b--','LineWidth',lin_wid)
set(gca,'YDir','reverse')
xlim([0 0.5])
ylim([p_ref(LevMin+1) p_ref(LevMax)])
title('q (kgkg^{-1})','FontSize',fontsize1,'FontName','TimesNewRoman')
set(gca,'YTickLabel',[])
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

subplot(1,8,5)
plot(corr_tl1(LevMin:LevMax,5),p_ref(LevMin:LevMax),'k','LineWidth',lin_wid)
hold on
plot(corr_tl2(LevMin:LevMax,5),p_ref(LevMin:LevMax),'r','LineWidth',lin_wid)
% plot(corr_tl3(LevMin:LevMax,5),p_ref(LevMin:LevMax),'b--','LineWidth',lin_wid)
set(gca,'YDir','reverse')
xlim([0 1])
ylim([p_ref(LevMin+1) p_ref(LevMax)])
title('p (Pa)','FontSize',fontsize1,'FontName','TimesNewRoman')
set(gca,'YTickLabel',[])
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

subplot(1,8,6)
plot(corr_tl1(LevMin:LevMax,6),p_ref(LevMin:LevMax),'k','LineWidth',lin_wid)
hold on
plot(corr_tl2(LevMin:LevMax,6),p_ref(LevMin:LevMax),'r','LineWidth',lin_wid)
% plot(corr_tl3(LevMin:LevMax,6),p_ref(LevMin:LevMax),'b--','LineWidth',lin_wid)
set(gca,'YDir','reverse')
xlim([0 0.5])
ylim([p_ref(LevMin+1) p_ref(LevMax)])
title('q_l (kgkg^{-1})','FontSize',fontsize1,'FontName','TimesNewRoman')
set(gca,'YTickLabel',[])
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

subplot(1,8,7)
plot(corr_tl1(LevMin:LevMax,7),p_ref(LevMin:LevMax),'k','LineWidth',lin_wid)
hold on
plot(corr_tl2(LevMin:LevMax,7),p_ref(LevMin:LevMax),'r','LineWidth',lin_wid)
% plot(corr_tl3(LevMin:LevMax,7),p_ref(LevMin:LevMax),'b--','LineWidth',lin_wid)
set(gca,'YDir','reverse')
xlim([0 0.5])
ylim([p_ref(LevMin+1) p_ref(LevMax)])
title('q_i (kgkg^{-1})','FontSize',fontsize1,'FontName','TimesNewRoman')
set(gca,'YTickLabel',[])
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

subplot(1,8,8)
plot(corr_tl1(LevMin:LevMax,8),p_ref(LevMin:LevMax),'k','LineWidth',lin_wid)
hold on
plot(corr_tl2(LevMin:LevMax,8),p_ref(LevMin:LevMax),'r','LineWidth',lin_wid)
% plot(corr_tl3(LevMin:LevMax,8),p_ref(LevMin:LevMax),'b--','LineWidth',lin_wid)
set(gca,'YDir','reverse')
xlim([0 0.5])
ylim([p_ref(LevMin+1) p_ref(LevMax)])
title('o3 (ppm)','FontSize',fontsize1,'FontName','TimesNewRoman')
ylabel('Height (hPa)','FontSize',fontsize1,'FontName','TimesNewRoman')
legend('Dry Old','Dry New','Moist')

set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'YAxisLocation','right')




% figure
% set(gcf,'position',[3 343 1276 576])
% 
% subplot(1,8,1)
% plot(corr_tl4(LevMin:LevMax,1),p_ref(LevMin:LevMax),'k','LineWidth',lin_wid)
% hold on
% plot(corr_tl5(LevMin:LevMax,1),p_ref(LevMin:LevMax),'r','LineWidth',lin_wid)
% plot(corr_tl6(LevMin:LevMax,1),p_ref(LevMin:LevMax),'b--','LineWidth',lin_wid)
% set(gca,'YDir','reverse')
% xlim([0 1])
% ylim([p_ref(LevMin+1) p_ref(LevMax)])
% title('u (ms^{-1})','FontSize',fontsize1,'FontName','TimesNewRoman')
% ylabel('Height (hPa)','FontSize',fontsize1,'FontName','TimesNewRoman')
% set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
% 
% subplot(1,8,2)
% plot(corr_tl4(LevMin:LevMax,2),p_ref(LevMin:LevMax),'k','LineWidth',lin_wid)
% hold on
% plot(corr_tl5(LevMin:LevMax,2),p_ref(LevMin:LevMax),'r','LineWidth',lin_wid)
% plot(corr_tl6(LevMin:LevMax,2),p_ref(LevMin:LevMax),'b--','LineWidth',lin_wid)
% set(gca,'YDir','reverse')
% xlim([0 1])
% ylim([p_ref(LevMin+1) p_ref(LevMax)])
% title('v (ms^{-1})','FontSize',fontsize1,'FontName','TimesNewRoman')
% set(gca,'YTickLabel',[])
% set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
% 
% subplot(1,8,3)
% plot(corr_tl4(LevMin:LevMax,3),p_ref(LevMin:LevMax),'k','LineWidth',lin_wid)
% hold on
% plot(corr_tl5(LevMin:LevMax,3),p_ref(LevMin:LevMax),'r','LineWidth',lin_wid)
% plot(corr_tl6(LevMin:LevMax,3),p_ref(LevMin:LevMax),'b--','LineWidth',lin_wid)
% set(gca,'YDir','reverse')
% xlim([0 1])
% ylim([p_ref(LevMin+1) p_ref(LevMax)])
% title('T_v (K)','FontSize',fontsize1,'FontName','TimesNewRoman')
% set(gca,'YTickLabel',[])
% set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
% 
% subplot(1,8,4)
% plot(corr_tl4(LevMin:LevMax,4),p_ref(LevMin:LevMax),'k','LineWidth',lin_wid)
% hold on
% plot(corr_tl5(LevMin:LevMax,4),p_ref(LevMin:LevMax),'r','LineWidth',lin_wid)
% plot(corr_tl6(LevMin:LevMax,4),p_ref(LevMin:LevMax),'b--','LineWidth',lin_wid)
% set(gca,'YDir','reverse')
% xlim([0 0.5])
% ylim([p_ref(LevMin+1) p_ref(LevMax)])
% title('q (kgkg^{-1})','FontSize',fontsize1,'FontName','TimesNewRoman')
% set(gca,'YTickLabel',[])
% set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
% 
% subplot(1,8,5)
% plot(corr_tl4(LevMin:LevMax,5),p_ref(LevMin:LevMax),'k','LineWidth',lin_wid)
% hold on
% plot(corr_tl5(LevMin:LevMax,5),p_ref(LevMin:LevMax),'r','LineWidth',lin_wid)
% plot(corr_tl6(LevMin:LevMax,5),p_ref(LevMin:LevMax),'b--','LineWidth',lin_wid)
% set(gca,'YDir','reverse')
% xlim([0 1])
% ylim([p_ref(LevMin+1) p_ref(LevMax)])
% title('p (Pa)','FontSize',fontsize1,'FontName','TimesNewRoman')
% set(gca,'YTickLabel',[])
% set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
% 
% subplot(1,8,6)
% plot(corr_tl4(LevMin:LevMax,6),p_ref(LevMin:LevMax),'k','LineWidth',lin_wid)
% hold on
% plot(corr_tl5(LevMin:LevMax,6),p_ref(LevMin:LevMax),'r','LineWidth',lin_wid)
% plot(corr_tl6(LevMin:LevMax,6),p_ref(LevMin:LevMax),'b--','LineWidth',lin_wid)
% set(gca,'YDir','reverse')
% xlim([0 0.5])
% ylim([p_ref(LevMin+1) p_ref(LevMax)])
% title('q_l (kgkg^{-1})','FontSize',fontsize1,'FontName','TimesNewRoman')
% set(gca,'YTickLabel',[])
% set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
% 
% subplot(1,8,7)
% plot(corr_tl4(LevMin:LevMax,7),p_ref(LevMin:LevMax),'k','LineWidth',lin_wid)
% hold on
% plot(corr_tl5(LevMin:LevMax,7),p_ref(LevMin:LevMax),'r','LineWidth',lin_wid)
% plot(corr_tl6(LevMin:LevMax,7),p_ref(LevMin:LevMax),'b--','LineWidth',lin_wid)
% set(gca,'YDir','reverse')
% xlim([0 0.5])
% ylim([p_ref(LevMin+1) p_ref(LevMax)])
% title('q_i (kgkg^{-1})','FontSize',fontsize1,'FontName','TimesNewRoman')
% set(gca,'YTickLabel',[])
% set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
% 
% subplot(1,8,8)
% plot(corr_tl4(LevMin:LevMax,8),p_ref(LevMin:LevMax),'k','LineWidth',lin_wid)
% hold on
% plot(corr_tl5(LevMin:LevMax,8),p_ref(LevMin:LevMax),'r','LineWidth',lin_wid)
% plot(corr_tl6(LevMin:LevMax,8),p_ref(LevMin:LevMax),'b--','LineWidth',lin_wid)
% set(gca,'YDir','reverse')
% xlim([0 0.5])
% ylim([p_ref(LevMin+1) p_ref(LevMax)])
% title('o3 (ppm)','FontSize',fontsize1,'FontName','TimesNewRoman')
% ylabel('Height (hPa)','FontSize',fontsize1,'FontName','TimesNewRoman')
% legend('Dry Old','Dry New','Moist')
% 
% set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
% set(gca,'YAxisLocation','right')
