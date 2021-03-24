close all
clear
clc

%Choose day to examine
date_start = datenum(2014, 02, 01, 00, 00, 00);

%Choose number of days to examine, if more than 1 then average plotted.
num_days = 1;

%Level to plot
plotlev = 63;

%Vertical resolution
lm = 72;

load coast
coast_lat = lat; clear lat
coast_lon = long; clear long

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
    dir = ['/discover/nobackup/drholdaw/v001_C180/prog/Y',dater(1:4),'/M',dater(5:6),'/D',dater(7:8),'/H21-NLM-free/']
    cd(dir)
    file = ['v001_C180.prog.eta.',dater(1:8),'_21z+',dateplot(1:8),'_',dateplot(9:10),'z.nc4']

    lon = ncread(file,'lon');
    lat = ncread(file,'lat');
    lev = ncread(file,'lev');

    u_free = ncread(file,'u');

    %Load Perturbed (analysis) state
    dir = ['/discover/nobackup/drholdaw/v001_C180/prog/Y',dater(1:4),'/M',dater(5:6),'/D',dater(7:8),'/H21-NLM-replay/']
    cd(dir)
    file = ['v001_C180.prog.eta.',dater(1:8),'_21z+',dateplot(1:8),'_',dateplot(9:10),'z.nc4'];

    u_replay = ncread(file,'u');

    %Load TLM state to compare.
    dir = ['/discover/nobackup/drholdaw/v001_C180/prog/Y',dater(1:4),'/M',dater(5:6),'/D',dater(7:8),'/H21-TLM-GQ1-SHAP0-MOIST0-NEWDYN0/'];
    cd(dir)
    file = ['v001_C180.fvpert.eta.',dater(1:8),'_21z+',dateplot(1:8),'_',dateplot(9:10),'z.nc4'];

    u_tl1 = ncread(file,'U');

    %Load TLM state to compare.
    dir = ['/discover/nobackup/drholdaw/v001_C180/prog/Y',dater(1:4),'/M',dater(5:6),'/D',dater(7:8),'/H21-TLM-GQ1-SHAP0-MOIST0-NEWDYN1/'];
    cd(dir)
    file = ['v001_C180.fvpert.eta.',dater(1:8),'_21z+',dateplot(1:8),'_',dateplot(9:10),'z.nc4'];

    u_tl2 = ncread(file,'U');

%     %Load TLM state to compare.
%     dir = ['/discover/nobackup/drholdaw/v001_C180/prog/Y',dater(1:4),'/M',dater(5:6),'/D',dater(7:8),'/H21-TLM-GQ1-PH2/'];
%     cd(dir)
%     file = ['v001_C180.fvpert.eta.',dater(1:8),'_21z+',dateplot(1:8),'_',dateplot(9:10),'z.nc4'];
% 
%     u_tl3 = ncread(file,'u');
% 
%     %Load TLM state to compare.
%     dir = ['/discover/nobackup/drholdaw/v001_C180/prog/Y',dater(1:4),'/M',dater(5:6),'/D',dater(7:8),'/H21-TLM-GQ2-PH0/'];
%     cd(dir)
%     file = ['v001_C180.fvpert.eta.',date00(1:8),'_00z+',dateplot(1:8),'_',dateplot(9:10),'z.nc4'];
% 
%     u_tl4 = ncread(file,'u');
% 
%     %Load TLM state to compare.
%     dir = ['/discover/nobackup/drholdaw/v001_C180/prog/Y',dater(1:4),'/M',dater(5:6),'/D',dater(7:8),'/H21-TLM-GQ2-PH1/'];
%     cd(dir)
%     file = ['v001_C180.fvpert.eta.',date00(1:8),'_00z+',dateplot(1:8),'_',dateplot(9:10),'z.nc4'];
% 
%     u_tl5 = ncread(file,'u');
% 
%     %Load TLM state to compare.
%     dir = ['/discover/nobackup/drholdaw/v001_C180/prog/Y',dater(1:4),'/M',dater(5:6),'/D',dater(7:8),'/H21-TLM-GQ2-PH2/'];
%     cd(dir)
%     file = ['v001_C180.fvpert.eta.',date00(1:8),'_00z+',dateplot(1:8),'_',dateplot(9:10),'z.nc4'];
% 
%     u_tl6 = ncread(file,'u');

    cd /home/drholdaw/Lin_Moist_Physics/DAS_plotting/Cloud' Fraction'/

    %Compute NL perturbation trajectory.
    u_nlm = u_replay - u_free;


    im = length(lon);
    jm = length(lat);
    lm = length(lev);

        
end
    
lin_wid = 1.75;
fontsize = 11;
fontsize1 = 12;

u_free_plot = u_free(:,:,plotlev);
u_repl_plot = u_replay(:,:,plotlev);
u_nlm_plot = u_nlm(:,:,plotlev);
u_tlm1_plot = u_tl1(:,:,plotlev);
u_tlm2_plot = u_tl2(:,:,plotlev);
% u_tlm3_plot = u_tl3(:,:,plotlev);
% u_tlm4_plot = u_tl4(:,:,plotlev);
% u_tlm5_plot = u_tl5(:,:,plotlev);
% u_tlm6_plot = u_tl6(:,:,plotlev);

% figure
% contourf(lon,lat,u_free_plot')
% 
% figure
% contourf(lon,lat,u_repl_plot')

grey = [0.75 0.75 0.75];
cmap = colormap;
cmap(32,:) = [1 1 1];
cmap(33,:) = [1 1 1];
close


figure
set(gcf,'position',[719    30   560   889])
subplot(3,1,1)
[~,h] = contourf(lon,lat,u_nlm_plot','LineStyle','none');
hold on
plot(coast_lon,coast_lat,'Color',grey)
ls = get(h,'LevelStep');
set(h,'LevelStep',ls/10)
colormap(cmap)
caxis([-max(abs(u_nlm_plot(:))) max(abs(u_nlm_plot(:)))])
colorbar
subplot(3,1,2)
[~,h] = contourf(lon,lat,u_tlm1_plot','LineStyle','none');
hold on
plot(coast_lon,coast_lat,'Color',grey)
ls = get(h,'LevelStep');
set(h,'LevelStep',ls/10)
colormap(cmap)
caxis([-max(abs(u_tlm2_plot(:))) max(abs(u_tlm2_plot(:)))])
colorbar
subplot(3,1,3)
[~,h] = contourf(lon,lat,u_tlm2_plot','LineStyle','none');
hold on
plot(coast_lon,coast_lat,'Color',grey)
ls = get(h,'LevelStep');
set(h,'LevelStep',ls/10)
colormap(cmap)
caxis([-max(abs(u_tlm2_plot(:))) max(abs(u_tlm2_plot(:)))])
colorbar
