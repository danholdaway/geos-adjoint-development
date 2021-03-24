close all
clear
clc

%Choose day to examine
date_start = datenum(2013, 01, 03, 00, 00, 00);

%Choose number of days to examine, if more than 1 then average plotted.
num_days = 1;

%Level to plot
plotlev = 68;


load coast
coast_lat = lat; clear lat
coast_lon = long; clear long

i = 0;
    
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
cd /discover/nobackup/jgkim/4dan/Feb01_14
file = 'v000_C180.prog.eta.20140203_00z.nc4';

lon = ncread(file,'lon');
lat = ncread(file,'lat');
lev = ncread(file,'lev');

u_free = ncread(file,'u');
t_free = ncread(file,'tv');

%Load Perturbed (analysis) state
cd /discover/nobackup/jgkim/4dan/Feb01_14
file = 'v000_C180.prog.eta.20140203_00z_pert.nc4';

u_replay = ncread(file,'u');
t_replay = ncread(file,'tv');

%Load TLM state to compare.
cd /discover/nobackup/jgkim/4dan/Feb01_14
file = 'fvpert.eta.20140203_00z.nc4';

u_tl1 = ncread(file,'u');
t_tl1 = ncread(file,'tv');

%Load TLM state to compare.
cd /discover/nobackup/jgkim/4dan/Feb01_14
file = 'fvpert.eta.20140203_00z_2001.nc4';

u_tl2 = ncread(file,'u');
t_tl2 = ncread(file,'tv');

cd /home/drholdaw/Lin_Moist_Physics/DAS_plotting/Cloud' Fraction'/

        
u_nlm = u_replay - u_free;
t_nlm = t_replay - t_free;

%Vertical resolution
im = length(lon);
jm = length(lat);
lm = length(lev);

% Uncomment to switch to RMS of the error.
% u_tl1 = u_tl1 - u_nlm;
% u_tl2 = u_tl2 - u_nlm;

u_nlm_rms  = zeros(1,lm);
u_tl1_rms = zeros(1,lm);
u_tl2_rms = zeros(1,lm);

t_nlm_rms  = zeros(1,lm);
t_tl1_rms = zeros(1,lm);
t_tl2_rms = zeros(1,lm);

for k = 1:lm

    u_nlm_rms(k)  = sqrt(sum(sum(u_nlm(:,:,k).^2))/(im*jm));
    u_tl1_rms(k) = sqrt(sum(sum(u_tl1(:,:,k).^2))/(im*jm));
    u_tl2_rms(k) = sqrt(sum(sum(u_tl2(:,:,k).^2))/(im*jm));

    t_nlm_rms(k)  = sqrt(sum(sum(t_nlm(:,:,k).^2))/(im*jm));
    t_tl1_rms(k) = sqrt(sum(sum(t_tl1(:,:,k).^2))/(im*jm));
    t_tl2_rms(k) = sqrt(sum(sum(t_tl2(:,:,k).^2))/(im*jm));

end

figure
set(gcf,'position',[1315 26 1200 913])
subplot(1,2,1)
plot(u_nlm_rms,1:lm,'k','LineWidth',1.8)
hold on
plot(u_tl1_rms,1:lm,'b','LineWidth',1.8)
plot(u_tl2_rms,1:lm,'r','LineWidth',1.8)
plot(0:9,40*ones(1,10),'k--')
set(gca,'Ydir','reverse')
ylim([1 72])

pref = 0.5*(p_ref(1:end-1) + p_ref(2:end));

subplot(1,2,2)
plot(u_nlm_rms,pref,'k','LineWidth',1.8)
hold on
plot(u_tl1_rms,pref,'b','LineWidth',1.8)
plot(u_tl2_rms,pref,'r','LineWidth',1.8)
plot(0:9,pref(40)*ones(1,10),'k--')
set(gca,'Ydir','reverse')
ylim([pref(1) pref(72)])


figure
set(gcf,'position',[1315 26 1200 913])
subplot(1,2,1)
plot(t_nlm_rms,1:lm,'k','LineWidth',1.8)
hold on
plot(t_tl1_rms,1:lm,'b','LineWidth',1.8)
plot(t_tl2_rms,1:lm,'r','LineWidth',1.8)
plot(0:9,40*ones(1,10),'k--')
set(gca,'Ydir','reverse')
ylim([1 72])

pref = 0.5*(p_ref(1:end-1) + p_ref(2:end));

subplot(1,2,2)
plot(t_nlm_rms,pref,'k','LineWidth',1.8)
hold on
plot(t_tl1_rms,pref,'b','LineWidth',1.8)
plot(t_tl2_rms,pref,'r','LineWidth',1.8)
plot(0:9,pref(40)*ones(1,10),'k--')
set(gca,'Ydir','reverse')
ylim([pref(1) pref(72)])



