close all
clear
clc

%Choose day to examine
date_start = datenum(2014, 02, 02, 00, 00, 00);

%Choose lead time, e.g. 24, 36 or 48 hours
leadtime = 24;
    
%Choose number of days to examine, if more than 1 then average plotted.
num_days = 1;

%Resolution
im = 576;
jm = 361;
lm = 72;

LevMin = 30;
LevMax = lm;

%Choose Region,
% 1 = Global
% 2 = Tropics 23S to 23N, below 100hPa
% 3 = Northern hemisphere
% 4 = Southern hemisphere
region = 1;

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

    
%Get required dates
daten = date_start;

fprintf('Date is %s \n', datestr(daten, 'yyyymmddHHMM'))

dater = datestr(daten - 0.1250, 'yyyymmddHHMM');
date00 = datestr(daten + 0, 'yyyymmddHHMM');
dateplot = datestr(daten + leadtime/24, 'yyyymmddHHMM');

cd /home/drholdaw/Lin_Moist_Physics/Inputs/
pref

%Load Free (background) State.
dir = ['/discover/nobackup/drholdaw/v001_C180/prog/Y',dater(1:4),'/M',dater(5:6),'/D',dater(7:8),'/H21-NLM-free/'];
cd(dir)
file = ['v001_C180.prog.eta.',dater(1:8),'_21z+',dateplot(1:8),'_',dateplot(9:10),'z.nc4'];


lon = ncread(file,'lon');
lat = ncread(file,'lat');
lev = ncread(file,'lev');


u_free = ncread(file,'u',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
v_free = ncread(file,'v',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
t_free = ncread(file,'tv',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
q_free = ncread(file,'sphu',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);

%Load Perturbed (analysis) state
dir = ['/discover/nobackup/drholdaw/v001_C180/prog/Y',dater(1:4),'/M',dater(5:6),'/D',dater(7:8),'/H21-NLM-replay/'];
cd(dir)
file = ['v001_C180.prog.eta.',dater(1:8),'_21z+',dateplot(1:8),'_',dateplot(9:10),'z.nc4'];

u_replay = ncread(file,'u',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
v_replay = ncread(file,'v',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
t_replay = ncread(file,'tv',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
q_replay = ncread(file,'sphu',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);

%Load TLM state to compare.
dir = ['/discover/nobackup/drholdaw/v001_C180/prog/Y',dater(1:4),'/M',dater(5:6),'/D',dater(7:8),'/H21-TLM-GQ1-SHAP0-MOIST0-NEWDYN0/'];
cd(dir)
file = ['v001_C180.fvpert.eta.',dater(1:8),'_21z+',dateplot(1:8),'_',dateplot(9:10),'z.nc4'];

u_tl1 = ncread(file,'U',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
v_tl1 = ncread(file,'V',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
t_tl1 = ncread(file,'TV',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
q_tl1 = ncread(file,'QV',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);

%Load TLM state to compare.
dir = ['/discover/nobackup/drholdaw/v001_C180/prog/Y',dater(1:4),'/M',dater(5:6),'/D',dater(7:8),'/H21-TLM-GQ1-SHAP0-MOIST0-NEWDYN1/'];
cd(dir)
file = ['v001_C180.fvpert.eta.',dater(1:8),'_21z+',dateplot(1:8),'_',dateplot(9:10),'z.nc4'];

u_tl2 = ncread(file,'U',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
v_tl2 = ncread(file,'V',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
t_tl2 = ncread(file,'TV',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
q_tl2 = ncread(file,'QV',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);

%Load TLM state to compare.
dir = ['/discover/nobackup/drholdaw/v001_C180/prog/Y',dater(1:4),'/M',dater(5:6),'/D',dater(7:8),'/H21-TLM-GQ1-SHAP1-MOIST0-NEWDYN1/'];
cd(dir)
file = ['v001_C180.fvpert.eta.',dater(1:8),'_21z+',dateplot(1:8),'_',dateplot(9:10),'z.nc4'];

u_tl3 = ncread(file,'U',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
v_tl3 = ncread(file,'V',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
t_tl3 = ncread(file,'TV',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
q_tl3 = ncread(file,'QV',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);

cd /home/drholdaw/Lin_Moist_Physics/DAS_plotting/Cloud' Fraction'/

%Compute NL perturbation trajectory.
u_nlm = u_replay - u_free;
v_nlm = v_replay - v_free;
t_nlm = t_replay - t_free;
q_nlm = q_replay - q_free;



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

%Variance, covariance and correlation
%Second index is number of variables.
var_nlm = zeros(lm,8);
var_tl1 = zeros(lm,8);
var_tl2 = zeros(lm,8);
var_tl3 = zeros(lm,8);

cov_tl1 = zeros(lm,8);
cov_tl2 = zeros(lm,8);
cov_tl3 = zeros(lm,8);

%Initialize correlations
corr_tl1 = zeros(lm,8);
corr_tl2 = zeros(lm,8);
corr_tl3 = zeros(lm,8);

% U winds
var_nlm(1:LevMax-LevMin+1,1) = sum(sum(u_nlm(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;
var_tl1(1:LevMax-LevMin+1,1) = sum(sum(u_tl1(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;
var_tl2(1:LevMax-LevMin+1,1) = sum(sum(u_tl2(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;
var_tl3(1:LevMax-LevMin+1,1) = sum(sum(u_tl3(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;

cov_tl1(1:LevMax-LevMin+1,1) = sum(sum(u_nlm(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*u_tl1(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;
cov_tl2(1:LevMax-LevMin+1,1) = sum(sum(u_nlm(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*u_tl2(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;
cov_tl3(1:LevMax-LevMin+1,1) = sum(sum(u_nlm(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*u_tl3(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;

corr_tl1(1:LevMax-LevMin+1,1) = cov_tl1(1:LevMax-LevMin+1,1)./sqrt(var_nlm(1:LevMax-LevMin+1,1).*var_tl1(1:LevMax-LevMin+1,1));
corr_tl2(1:LevMax-LevMin+1,1) = cov_tl2(1:LevMax-LevMin+1,1)./sqrt(var_nlm(1:LevMax-LevMin+1,1).*var_tl2(1:LevMax-LevMin+1,1));
corr_tl3(1:LevMax-LevMin+1,1) = cov_tl3(1:LevMax-LevMin+1,1)./sqrt(var_nlm(1:LevMax-LevMin+1,1).*var_tl3(1:LevMax-LevMin+1,1));

% V winds
var_nlm(1:LevMax-LevMin+1,2) = sum(sum(v_nlm(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;
var_tl1(1:LevMax-LevMin+1,2) = sum(sum(v_tl1(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;
var_tl2(1:LevMax-LevMin+1,2) = sum(sum(v_tl2(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;
var_tl3(1:LevMax-LevMin+1,2) = sum(sum(v_tl3(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;

cov_tl1(1:LevMax-LevMin+1,2) = sum(sum(v_nlm(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*v_tl1(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;
cov_tl2(1:LevMax-LevMin+1,2) = sum(sum(v_nlm(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*v_tl2(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;
cov_tl3(1:LevMax-LevMin+1,2) = sum(sum(v_nlm(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*v_tl3(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;

corr_tl1(1:LevMax-LevMin+1,2) = cov_tl1(1:LevMax-LevMin+1,2)./sqrt(var_nlm(1:LevMax-LevMin+1,2).*var_tl1(1:LevMax-LevMin+1,2));
corr_tl2(1:LevMax-LevMin+1,2) = cov_tl2(1:LevMax-LevMin+1,2)./sqrt(var_nlm(1:LevMax-LevMin+1,2).*var_tl2(1:LevMax-LevMin+1,2));
corr_tl3(1:LevMax-LevMin+1,2) = cov_tl3(1:LevMax-LevMin+1,2)./sqrt(var_nlm(1:LevMax-LevMin+1,2).*var_tl3(1:LevMax-LevMin+1,2));

% Temperature
var_nlm(1:LevMax-LevMin+1,3) = sum(sum(t_nlm(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;
var_tl1(1:LevMax-LevMin+1,3) = sum(sum(t_tl1(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;
var_tl2(1:LevMax-LevMin+1,3) = sum(sum(t_tl2(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;
var_tl3(1:LevMax-LevMin+1,3) = sum(sum(t_tl3(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;

cov_tl1(1:LevMax-LevMin+1,3) = sum(sum(t_nlm(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*t_tl1(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;
cov_tl2(1:LevMax-LevMin+1,3) = sum(sum(t_nlm(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*t_tl2(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;
cov_tl3(1:LevMax-LevMin+1,3) = sum(sum(t_nlm(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*t_tl3(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;

corr_tl1(1:LevMax-LevMin+1,3) = cov_tl1(1:LevMax-LevMin+1,3)./sqrt(var_nlm(1:LevMax-LevMin+1,3).*var_tl1(1:LevMax-LevMin+1,3));
corr_tl2(1:LevMax-LevMin+1,3) = cov_tl2(1:LevMax-LevMin+1,3)./sqrt(var_nlm(1:LevMax-LevMin+1,3).*var_tl2(1:LevMax-LevMin+1,3));
corr_tl3(1:LevMax-LevMin+1,3) = cov_tl3(1:LevMax-LevMin+1,3)./sqrt(var_nlm(1:LevMax-LevMin+1,3).*var_tl3(1:LevMax-LevMin+1,3));

% Specific Humidity
var_nlm(1:LevMax-LevMin+1,4) = sum(sum(q_nlm(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;
var_tl1(1:LevMax-LevMin+1,4) = sum(sum(q_tl1(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;
var_tl2(1:LevMax-LevMin+1,4) = sum(sum(q_tl2(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;
var_tl3(1:LevMax-LevMin+1,4) = sum(sum(q_tl3(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;

cov_tl1(1:LevMax-LevMin+1,4) = sum(sum(q_nlm(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*q_tl1(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;
cov_tl2(1:LevMax-LevMin+1,4) = sum(sum(q_nlm(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*q_tl2(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;
cov_tl3(1:LevMax-LevMin+1,4) = sum(sum(q_nlm(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*q_tl3(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;

corr_tl1(1:LevMax-LevMin+1,4) = cov_tl1(1:LevMax-LevMin+1,4)./sqrt(var_nlm(1:LevMax-LevMin+1,4).*var_tl1(1:LevMax-LevMin+1,4));
corr_tl2(1:LevMax-LevMin+1,4) = cov_tl2(1:LevMax-LevMin+1,4)./sqrt(var_nlm(1:LevMax-LevMin+1,4).*var_tl2(1:LevMax-LevMin+1,4));
corr_tl3(1:LevMax-LevMin+1,4) = cov_tl3(1:LevMax-LevMin+1,4)./sqrt(var_nlm(1:LevMax-LevMin+1,4).*var_tl3(1:LevMax-LevMin+1,4));

fprintf('  >Done correlation calculation \n\n\n')
    
    

lin_wid = 1.75;
fontsize = 11;
fontsize1 = 12;

p_ref = 0.5*(p_ref(1:end-1) + p_ref(2:end));

figure
set(gcf,'position',[3 343 1276 576])

subplot(1,4,1)
plot(corr_tl1(1:LevMax-LevMin+1,1),p_ref(LevMin:LevMax),'b','LineWidth',lin_wid)
hold on
plot(corr_tl2(1:LevMax-LevMin+1,1),p_ref(LevMin:LevMax),'r','LineWidth',lin_wid)
plot(corr_tl3(1:LevMax-LevMin+1,1),p_ref(LevMin:LevMax),'g','LineWidth',lin_wid)
set(gca,'YDir','reverse')
xlim([0 1])
ylim([p_ref(LevMin+1) p_ref(LevMax)])
title('u (ms^{-1})','FontSize',fontsize1,'FontName','TimesNewRoman')
ylabel('Height (hPa)','FontSize',fontsize1,'FontName','TimesNewRoman')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

subplot(1,4,2)
plot(corr_tl1(1:LevMax-LevMin+1,2),p_ref(LevMin:LevMax),'b','LineWidth',lin_wid)
hold on
plot(corr_tl2(1:LevMax-LevMin+1,2),p_ref(LevMin:LevMax),'r','LineWidth',lin_wid)
plot(corr_tl3(1:LevMax-LevMin+1,2),p_ref(LevMin:LevMax),'g','LineWidth',lin_wid)
set(gca,'YDir','reverse')
xlim([0 1])
ylim([p_ref(LevMin+1) p_ref(LevMax)])
title('v (ms^{-1})','FontSize',fontsize1,'FontName','TimesNewRoman')
set(gca,'YTickLabel',[])
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

subplot(1,4,3)
plot(corr_tl1(1:LevMax-LevMin+1,3),p_ref(LevMin:LevMax),'b','LineWidth',lin_wid)
hold on
plot(corr_tl2(1:LevMax-LevMin+1,3),p_ref(LevMin:LevMax),'r','LineWidth',lin_wid)
plot(corr_tl3(1:LevMax-LevMin+1,3),p_ref(LevMin:LevMax),'g','LineWidth',lin_wid)
set(gca,'YDir','reverse')
xlim([0 1])
ylim([p_ref(LevMin+1) p_ref(LevMax)])
title('T_v (K)','FontSize',fontsize1,'FontName','TimesNewRoman')
set(gca,'YTickLabel',[])
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

subplot(1,4,4)
plot(corr_tl1(1:LevMax-LevMin+1,4),p_ref(LevMin:LevMax),'b','LineWidth',lin_wid)
hold on
plot(corr_tl2(1:LevMax-LevMin+1,4),p_ref(LevMin:LevMax),'r','LineWidth',lin_wid)
plot(corr_tl3(1:LevMax-LevMin+1,4),p_ref(LevMin:LevMax),'g','LineWidth',lin_wid)
set(gca,'YDir','reverse')
xlim([0 0.5])
ylim([p_ref(LevMin+1) p_ref(LevMax)])
title('q (kgkg^{-1})','FontSize',fontsize1,'FontName','TimesNewRoman')
set(gca,'YTickLabel',[])
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

legend('Control','New Dyn','Shapiro')
