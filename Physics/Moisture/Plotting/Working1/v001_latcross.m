close all
clear
clc
mydir = pwd;

%Choose day to examine
date_start = datenum(2014, 02, 01, 00, 00, 00);
date_end   = datenum(2014, 02, 07, 00, 00, 00);

num_days = date_end - date_start + 1;

%Choose lead time, e.g. 24, 36 or 48 hours
leadtime = 24;

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

rmsd_tl1(1:LevMax-LevMin+1,1:4) = 0.0;
rmsd_tl2(1:LevMax-LevMin+1,1:4) = 0.0;
rmsd_tl3(1:LevMax-LevMin+1,1:4) = 0.0;

for i = 0:num_days-1
    
    %Get required dates
    daten = date_start + i;
        
    fprintf('Reading data for %s \n', datestr(daten, 'yyyymmddHHMM'))
    
    dater = datestr(daten - 0.1250, 'yyyymmddHHMM');
    date00 = datestr(daten + 0, 'yyyymmddHHMM');
    dateplot = datestr(daten + leadtime/24, 'yyyymmddHHMM');

    cd /home/drholdaw/LinearisedPhysics/Inputs/
    pref

    %Load Free (background) State.
    dir = ['/archive/u/drholdaw/v001_C180/prog/Y',dater(1:4),'/M',dater(5:6),'/D',dater(7:8),'/H21-NLM-free/'];
    cd(dir)
    file = ['v001_C180.prog.eta.',dater(1:8),'_21z+',dateplot(1:8),'_',dateplot(9:10),'z.nc4'];

    if i == 0
        lon = ncread(file,'lon');
        lat = ncread(file,'lat');
        lev = ncread(file,'lev');
    end
    
    u_free = ncread(file,'u',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
    v_free = ncread(file,'v',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
    t_free = ncread(file,'tv',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
    q_free = ncread(file,'sphu',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);

    %Load Perturbed (analysis) state
    dir = ['/archive/u/drholdaw/v001_C180/prog/Y',dater(1:4),'/M',dater(5:6),'/D',dater(7:8),'/H21-NLM-replay/'];
    cd(dir)
    file = ['v001_C180.prog.eta.',dater(1:8),'_21z+',dateplot(1:8),'_',dateplot(9:10),'z.nc4'];

    u_replay = ncread(file,'u',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
    v_replay = ncread(file,'v',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
    t_replay = ncread(file,'tv',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
    q_replay = ncread(file,'sphu',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
    
    %Load TLM state to compare.
    dir = ['/archive/u/drholdaw/v001_C180/prog/Y',dater(1:4),'/M',dater(5:6),'/D',dater(7:8),'/H21-TLM-GQ1-SHAP0-MOIST0-NEWDYN0/'];
    cd(dir)
    file = ['v001_C180.fvpert.eta.',dater(1:8),'_21z+',dateplot(1:8),'_',dateplot(9:10),'z.nc4'];

    u_tl1 = ncread(file,'U',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
    v_tl1 = ncread(file,'V',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
    t_tl1 = ncread(file,'TV',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
    q_tl1 = ncread(file,'QV',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);

    %Load TLM state to compare.
    dir = ['/archive/u/drholdaw/v001_C180/prog/Y',dater(1:4),'/M',dater(5:6),'/D',dater(7:8),'/H21-TLM-GQ1-SHAP0-MOIST0-NEWDYN1/'];
    cd(dir)
    file = ['v001_C180.fvpert.eta.',dater(1:8),'_21z+',dateplot(1:8),'_',dateplot(9:10),'z.nc4'];

    u_tl2 = ncread(file,'U',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
    v_tl2 = ncread(file,'V',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
    t_tl2 = ncread(file,'TV',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
    q_tl2 = ncread(file,'QV',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);

    %Load TLM state to compare.
    dir = ['/archive/u/drholdaw/v001_C180/prog/Y',dater(1:4),'/M',dater(5:6),'/D',dater(7:8),'/H21-TLM-GQ1-SHAP1-MOIST0-NEWDYN1/'];
    cd(dir)
    file = ['v001_C180.fvpert.eta.',dater(1:8),'_21z+',dateplot(1:8),'_',dateplot(9:10),'z.nc4'];

    u_tl3 = ncread(file,'U',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
    v_tl3 = ncread(file,'V',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
    t_tl3 = ncread(file,'TV',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
    q_tl3 = ncread(file,'QV',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);

    cd(mydir)

    %Compute NL perturbation trajectory.
    u_nlm = u_replay - u_free;
    v_nlm = v_replay - v_free;
    t_nlm = t_replay - t_free;
    q_nlm = q_replay - q_free;
   

    
end
 
cint = 20;


u_nlm_latcross = zeros(jm, LevMax-LevMin+1);
u_nlm_latcross(:,:) = mean(u_nlm,1);

% Uncomment below to plot the actual wind field
% u_nlm_latcross(:,:) = mean(u_free,1);

u_tl1_latcross = zeros(jm, LevMax-LevMin+1);
u_tl1_latcross(:,:) = mean(u_tl1,1);

u_tl2_latcross = zeros(jm, LevMax-LevMin+1);
u_tl2_latcross(:,:) = mean(u_tl2,1);

u_tl3_latcross = zeros(jm, LevMax-LevMin+1);
u_tl3_latcross(:,:) = mean(u_tl3,1);

u_tl1_latcross_diff = zeros(jm, LevMax-LevMin+1);
u_tl1_latcross_diff(:,:) = mean(u_tl1 - u_nlm,1);

u_tl2_latcross_diff = zeros(jm, LevMax-LevMin+1);
u_tl2_latcross_diff(:,:) = mean(u_tl2 - u_nlm,1);

u_tl3_latcross_diff = zeros(jm, LevMax-LevMin+1);
u_tl3_latcross_diff(:,:) = mean(u_tl3 - u_nlm,1);

maxu = max([abs(u_nlm_latcross(:)) ; abs(u_tl1_latcross(:)) ; abs(u_tl2_latcross(:)) ; abs(u_tl3_latcross(:)) ] );
maxu_diff = max([abs(u_nlm_latcross(:)) ; abs(u_tl1_latcross_diff(:)) ; abs(u_tl2_latcross_diff(:)) ; abs(u_tl3_latcross_diff(:)) ] );
cint_u = 2*maxu/cint;
cint_u_diff = 2*maxu_diff/cint;

t_nlm_latcross = zeros(jm, LevMax-LevMin+1);
t_nlm_latcross(:,:) = mean(t_nlm,1);

t_tl1_latcross = zeros(jm, LevMax-LevMin+1);
t_tl1_latcross(:,:) = mean(t_tl1,1);

t_tl2_latcross = zeros(jm, LevMax-LevMin+1);
t_tl2_latcross(:,:) = mean(t_tl2,1);

t_tl3_latcross = zeros(jm, LevMax-LevMin+1);
t_tl3_latcross(:,:) = mean(t_tl3,1);

t_tl1_latcross_diff = zeros(jm, LevMax-LevMin+1);
t_tl1_latcross_diff(:,:) = mean(t_tl1 - t_nlm,1);

t_tl2_latcross_diff = zeros(jm, LevMax-LevMin+1);
t_tl2_latcross_diff(:,:) = mean(t_tl2 - t_nlm,1);

t_tl3_latcross_diff = zeros(jm, LevMax-LevMin+1);
t_tl3_latcross_diff(:,:) = mean(t_tl3 - t_nlm,1);

maxt = max([abs(t_nlm_latcross(:)) ; abs(t_tl1_latcross(:)) ; abs(t_tl2_latcross(:)) ; abs(t_tl3_latcross(:)) ] );
maxt_diff = max([abs(t_nlm_latcross(:)) ; abs(t_tl1_latcross_diff(:)) ; abs(t_tl2_latcross_diff(:)) ; abs(t_tl3_latcross_diff(:)) ] );
cint_t = 2*maxt/cint;
cint_t_diff = 2*maxt_diff/cint;



%Make center of colormap white
cmap = colormap;
cmap(31,:) = [1 1 1];
cmap(32,:) = [1 1 1];
cmap(33,:) = [1 1 1];
cmap(34,:) = [1 1 1];
close


figure
set(gcf,'position',[3 343 1276 576])

subplot(2,1,1)
[~,h] = contourf(lat,LevMin:LevMax,t_nlm_latcross','LineStyle','none');
caxis([-maxt maxt])
colorbar
colormap(cmap)
set(gca,'YDir','Reverse')
set(h,'LevelStep',cint_t);
title('T_v - nonlinear difference')

subplot(2,1,2)
[~,h] = contourf(lat,LevMin:LevMax,t_tl1_latcross','LineStyle','none');
caxis([-maxt maxt])
colorbar
colormap(cmap)
set(gca,'YDir','Reverse')
set(h,'LevelStep',cint_t);
title('T_v - tlm (default)')



% 
% figure
% subplot(4,1,1)
% [~,h] = contourf(lat,LevMin:LevMax,u_nlm_latcross','LineStyle','none');
% caxis([-maxu maxu])
% colorbar
% colormap(cmap)
% set(gca,'YDir','Reverse')
% set(h,'LevelStep',cint_u);
% title('u wind - nonlinear difference')
% 
% subplot(4,1,2)
% [~,h] = contourf(lat,LevMin:LevMax,u_tl1_latcross','LineStyle','none');
% caxis([-maxu maxu])
% colorbar
% colormap(cmap)
% set(gca,'YDir','Reverse')
% set(h,'LevelStep',cint_u);
% title('u wind - tlm (default)')
% 
% subplot(4,1,3)
% [~,h] = contourf(lat,LevMin:LevMax,u_tl2_latcross','LineStyle','none');
% caxis([-maxu maxu])
% colorbar
% colormap(cmap)
% set(gca,'YDir','Reverse')
% set(h,'LevelStep',cint_u);
% title('u wind - tlm (new dynamics)')
% 
% subplot(4,1,4)
% [~,h] = contourf(lat,LevMin:LevMax,u_tl3_latcross','LineStyle','none');
% caxis([-maxu maxu])
% colorbar
% colormap(cmap)
% set(gca,'YDir','Reverse')
% set(h,'LevelStep',cint_u);
% title('u wind - tlm (new dynamics + shapiro filter)')
% 
% 
% 
% figure
% subplot(4,1,1)
% [~,h] = contourf(lat,LevMin:LevMax,u_nlm_latcross','LineStyle','none');
% caxis([-maxu_diff maxu_diff])
% colorbar
% colormap(cmap)
% set(gca,'YDir','Reverse')
% set(h,'LevelStep',cint_u_diff);
% title('u wind - nonlinear difference')
% 
% subplot(4,1,2)
% [~,h] = contourf(lat,LevMin:LevMax,u_tl1_latcross_diff','LineStyle','none');
% caxis([-maxu_diff maxu_diff])
% colorbar
% colormap(cmap)
% set(gca,'YDir','Reverse')
% set(h,'LevelStep',cint_u_diff);
% title('u wind - mean(nlm - tlm default)')
% 
% subplot(4,1,3)
% [~,h] = contourf(lat,LevMin:LevMax,u_tl2_latcross_diff','LineStyle','none');
% caxis([-maxu_diff maxu_diff])
% colorbar
% colormap(cmap)
% set(gca,'YDir','Reverse')
% set(h,'LevelStep',cint_u_diff);
% title('u wind - mean(nlm - tlm new dynamics)')
% 
% subplot(4,1,4)
% [~,h] = contourf(lat,LevMin:LevMax,u_tl3_latcross_diff','LineStyle','none');
% caxis([-maxu_diff maxu_diff])
% colorbar
% colormap(cmap)
% set(gca,'YDir','Reverse')
% set(h,'LevelStep',cint_u_diff);
% title('u wind - mean(nlm - tlm new dynamics + shapiro)')
% 
% 
% 
% 
% figure
% subplot(4,1,1)
% [~,h] = contourf(lat,LevMin:LevMax,t_nlm_latcross','LineStyle','none');
% caxis([-maxt maxt])
% colorbar
% colormap(cmap)
% set(gca,'YDir','Reverse')
% set(h,'LevelStep',cint_t);
% title('T_v - nonlinear difference')
% 
% subplot(4,1,2)
% [~,h] = contourf(lat,LevMin:LevMax,t_tl1_latcross','LineStyle','none');
% caxis([-maxt maxt])
% colorbar
% colormap(cmap)
% set(gca,'YDir','Reverse')
% set(h,'LevelStep',cint_t);
% title('T_v - tlm (default)')
% 
% subplot(4,1,3)
% [~,h] = contourf(lat,LevMin:LevMax,t_tl2_latcross','LineStyle','none');
% caxis([-maxt maxt])
% colorbar
% colormap(cmap)
% set(gca,'YDir','Reverse')
% set(h,'LevelStep',cint_t);
% title('T_v - tlm (new dynamics)')
% 
% subplot(4,1,4)
% [~,h] = contourf(lat,LevMin:LevMax,t_tl3_latcross','LineStyle','none');
% caxis([-maxt maxt])
% colorbar
% colormap(cmap)
% set(gca,'YDir','Reverse')
% set(h,'LevelStep',cint_t);
% title('T_v - tlm (new dynamics + shapiro filter)')
% 
% 
% 
% figure
% subplot(4,1,1)
% [~,h] = contourf(lat,LevMin:LevMax,t_nlm_latcross','LineStyle','none');
% caxis([-maxt_diff maxt_diff])
% colorbar
% colormap(cmap)
% set(gca,'YDir','Reverse')
% set(h,'LevelStep',cint_t_diff);
% title('T_v - nonlinear difference')
% 
% subplot(4,1,2)
% [~,h] = contourf(lat,LevMin:LevMax,t_tl1_latcross_diff','LineStyle','none');
% caxis([-maxt_diff maxt_diff])
% colorbar
% colormap(cmap)
% set(gca,'YDir','Reverse')
% set(h,'LevelStep',cint_t_diff);
% title('T_v - mean(nlm - tlm default)')
% 
% subplot(4,1,3)
% [~,h] = contourf(lat,LevMin:LevMax,t_tl2_latcross_diff','LineStyle','none');
% caxis([-maxt_diff maxt_diff])
% colorbar
% colormap(cmap)
% set(gca,'YDir','Reverse')
% set(h,'LevelStep',cint_t_diff);
% title('T_v - mean(nlm - tlm new dynamics)')
% 
% subplot(4,1,4)
% [~,h] = contourf(lat,LevMin:LevMax,t_tl3_latcross_diff','LineStyle','none');
% caxis([-maxt_diff maxt_diff])
% colorbar
% colormap(cmap)
% set(gca,'YDir','Reverse')
% set(h,'LevelStep',cint_t_diff);
% title('T_v - mean(nlm - tlm new dynamics + shapiro)')
