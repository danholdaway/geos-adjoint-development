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


cint = 20;

t_nlm_latcross_sum = zeros(jm, LevMax-LevMin+1);
t_tl1_latcross_sum = zeros(jm, LevMax-LevMin+1);

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

    t_free = ncread(file,'tv',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);

    %Load Perturbed (analysis) state
%     dir = ['/archive/u/drholdaw/v001_C180/prog/Y',dater(1:4),'/M',dater(5:6),'/D',dater(7:8),'/H21-NLM-replay/'];
%     cd(dir)
%     file = ['v001_C180.prog.eta.',dater(1:8),'_21z+',dateplot(1:8),'_',dateplot(9:10),'z.nc4'];

        dir = ['/archive/u/drholdaw/v001_C180/prog/Y',dater(1:4),'/M',dater(5:6),'/D',dater(7:8),'/H21-TLM-GQ1-SHAP0-MOIST0-NEWDYN0/'];
    cd(dir)
    file = ['v001_C180.fvpert.eta.',dater(1:8),'_21z+',dateplot(1:8),'_',dateplot(9:10),'z.nc4'];
    
    t_replay = ncread(file,'TV',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
    
    %Load TLM state to compare.
    dir = ['/archive/u/drholdaw/v001_C180/prog/Y',dater(1:4),'/M',dater(5:6),'/D',dater(7:8),'/H21-TLM-GQ1-SHAP1-MOIST0-NEWDYN1/'];
    cd(dir)
    file = ['v001_C180.fvpert.eta.',dater(1:8),'_21z+',dateplot(1:8),'_',dateplot(9:10),'z.nc4'];

    t_tl1 = ncread(file,'TV',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
    
    
    cd(mydir)

    %Compute NL perturbation trajectory.

    t_nlm = t_replay;
   

    t_nlm_latcross = zeros(jm, LevMax-LevMin+1);
    t_nlm_latcross(:,:) = mean(t_nlm,1);

    t_nlm_latcross_sum(:,:) = t_nlm_latcross_sum(:,:) + t_nlm_latcross(:,:);
    
    t_tl1_latcross = zeros(jm, LevMax-LevMin+1);
    t_tl1_latcross(:,:) = mean(t_tl1,1);
    
    t_tl1_latcross_sum(:,:) = t_tl1_latcross_sum(:,:) + t_tl1_latcross(:,:);
    
end
 
t_nlm_latcross_sum(:,:) = t_nlm_latcross_sum(:,:)/num_days;
t_tl1_latcross_sum(:,:) = t_tl1_latcross_sum(:,:)/num_days;


t_nlm_latcross_sum_plot = t_nlm_latcross_sum(:,:);
t_tl1_latcross_sum_plot = t_tl1_latcross_sum(:,:);

maxt = max([abs(t_nlm_latcross_sum_plot(:)) ; abs(t_tl1_latcross_sum_plot(:)) ] );
cint_t = 2*0.9031/30;



%Make center of colormap white
cmap = colormap;
cmap(31,:) = [1 1 1];
cmap(32,:) = [1 1 1];
cmap(33,:) = [1 1 1];
cmap(34,:) = [1 1 1];
close


fontsize = 14;

figure
set(gcf,'position',[3 343 1276 576])

subplot(2,1,1)
[~,h] = contourf(lat,LevMin:LevMax,t_nlm_latcross_sum_plot','LineStyle','none');
caxis([-maxt maxt])
colorbar
colormap(cmap)
set(gca,'YDir','Reverse')
set(h,'LevelStep',cint_t);
title('T_v (K) tangent linear model perturbation old dynamics','FontSize',fontsize,'FontName','TimesNewRoman')
xlabel('Latitude (\circ)','FontSize',fontsize,'FontName','TimesNewRoman')
ylabel('Model Level','FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
    
subplot(2,1,2)
[~,h] = contourf(lat,LevMin:LevMax,t_tl1_latcross_sum_plot','LineStyle','none');
caxis([-maxt maxt])
colorbar
colormap(cmap)
set(gca,'YDir','Reverse')
set(h,'LevelStep',cint_t);
title('T_v (K) tangent linear model perturbation new dynamics','FontSize',fontsize,'FontName','TimesNewRoman')
xlabel('Latitude (\circ)','FontSize',fontsize,'FontName','TimesNewRoman')
ylabel('Model Level','FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')



