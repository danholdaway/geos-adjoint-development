close all
clear
clc
mydir = pwd;

%Choose day to examine
date_start = datenum(2014, 02, 01, 00, 00, 00);
date_end   = datenum(2014, 02, 28, 00, 00, 00);

num_days = date_end - date_start + 1;

%Choose lead time, e.g. 24, 36 or 48 hours
leadtime = 24;

%Resolution
im = 576;
jm = 361;
lm = 72;

LevMin = 1;
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

    fprintf(' >Beginning RMS difference calculation \n')
    
    n = (LonMax-LonMin+1)*(LatMax-LatMin+1);

    % U winds
    d(1:LevMax-LevMin+1) = sum(sum((u_nlm(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1) - u_tl1(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1)).^2)) / n ;
    rmsd_tl1(1:LevMax-LevMin+1,1) = rmsd_tl1(1:LevMax-LevMin+1,1) + d(1:LevMax-LevMin+1)';
    d(1:LevMax-LevMin+1) = sum(sum((u_nlm(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1) - u_tl2(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1)).^2)) / n ;
    rmsd_tl2(1:LevMax-LevMin+1,1) = rmsd_tl2(1:LevMax-LevMin+1,1) + d(1:LevMax-LevMin+1)';
    d(1:LevMax-LevMin+1) = sum(sum((u_nlm(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1) - u_tl3(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1)).^2)) / n ;
    rmsd_tl3(1:LevMax-LevMin+1,1) = rmsd_tl3(1:LevMax-LevMin+1,1) + d(1:LevMax-LevMin+1)';

    % V winds
    d(1:LevMax-LevMin+1) = sum(sum((v_nlm(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1) - v_tl1(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1)).^2)) / n ;
    rmsd_tl1(1:LevMax-LevMin+1,2) = rmsd_tl1(1:LevMax-LevMin+1,2) + d(1:LevMax-LevMin+1)';
    d(1:LevMax-LevMin+1) = sum(sum((v_nlm(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1) - v_tl2(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1)).^2)) / n ;
    rmsd_tl2(1:LevMax-LevMin+1,2) = rmsd_tl2(1:LevMax-LevMin+1,2) + d(1:LevMax-LevMin+1)';
    d(1:LevMax-LevMin+1) = sum(sum((v_nlm(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1) - v_tl3(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1)).^2)) / n ;
    rmsd_tl3(1:LevMax-LevMin+1,2) = rmsd_tl3(1:LevMax-LevMin+1,2) + d(1:LevMax-LevMin+1)';

    % Temperature
    d(1:LevMax-LevMin+1) = sum(sum((t_nlm(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1) - t_tl1(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1)).^2)) / n ;
    rmsd_tl1(1:LevMax-LevMin+1,3) = rmsd_tl1(1:LevMax-LevMin+1,3) + d(1:LevMax-LevMin+1)';
    d(1:LevMax-LevMin+1) = sum(sum((t_nlm(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1) - t_tl2(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1)).^2)) / n ;
    rmsd_tl2(1:LevMax-LevMin+1,3) = rmsd_tl2(1:LevMax-LevMin+1,3) + d(1:LevMax-LevMin+1)';
    d(1:LevMax-LevMin+1) = sum(sum((t_nlm(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1) - t_tl3(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1)).^2)) / n ;
    rmsd_tl3(1:LevMax-LevMin+1,3) = rmsd_tl3(1:LevMax-LevMin+1,3) + d(1:LevMax-LevMin+1)';

    % Specific Humidity
    d(1:LevMax-LevMin+1) = sum(sum((q_nlm(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1) - q_tl1(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1)).^2)) / n ;
    rmsd_tl1(1:LevMax-LevMin+1,4) = rmsd_tl1(1:LevMax-LevMin+1,4) + d(1:LevMax-LevMin+1)';
    d(1:LevMax-LevMin+1) = sum(sum((q_nlm(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1) - q_tl2(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1)).^2)) / n ;
    rmsd_tl2(1:LevMax-LevMin+1,4) = rmsd_tl2(1:LevMax-LevMin+1,4) + d(1:LevMax-LevMin+1)';
    d(1:LevMax-LevMin+1) = sum(sum((q_nlm(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1) - q_tl3(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1)).^2)) / n ;
    rmsd_tl3(1:LevMax-LevMin+1,4) = rmsd_tl3(1:LevMax-LevMin+1,4) + d(1:LevMax-LevMin+1)';
 
    fprintf(' >Done RMS difference calculation \n\n')
    
end
    
rmsd_tl1(1:LevMax-LevMin+1,1) = rmsd_tl1(1:LevMax-LevMin+1,1)/num_days;
rmsd_tl2(1:LevMax-LevMin+1,1) = rmsd_tl2(1:LevMax-LevMin+1,1)/num_days;
rmsd_tl3(1:LevMax-LevMin+1,1) = rmsd_tl3(1:LevMax-LevMin+1,1)/num_days;

rmsd_tl1(1:LevMax-LevMin+1,2) = rmsd_tl1(1:LevMax-LevMin+1,2)/num_days;
rmsd_tl2(1:LevMax-LevMin+1,2) = rmsd_tl2(1:LevMax-LevMin+1,2)/num_days;
rmsd_tl3(1:LevMax-LevMin+1,2) = rmsd_tl3(1:LevMax-LevMin+1,2)/num_days;

rmsd_tl1(1:LevMax-LevMin+1,3) = rmsd_tl1(1:LevMax-LevMin+1,3)/num_days;
rmsd_tl2(1:LevMax-LevMin+1,3) = rmsd_tl2(1:LevMax-LevMin+1,3)/num_days;
rmsd_tl3(1:LevMax-LevMin+1,3) = rmsd_tl3(1:LevMax-LevMin+1,3)/num_days;

rmsd_tl1(1:LevMax-LevMin+1,4) = rmsd_tl1(1:LevMax-LevMin+1,4)/num_days;
rmsd_tl2(1:LevMax-LevMin+1,4) = rmsd_tl2(1:LevMax-LevMin+1,4)/num_days;
rmsd_tl3(1:LevMax-LevMin+1,4) = rmsd_tl3(1:LevMax-LevMin+1,4)/num_days;

lin_wid = 1.75;
fontsize = 11;
fontsize1 = 12;

p_ref = 0.5*(p_ref(1:end-1) + p_ref(2:end));

figure
set(gcf,'position',[3 343 1276 576])

subplot(1,4,1)
plot(rmsd_tl1(1:LevMax-LevMin+1,1),p_ref(LevMin:LevMax),'b','LineWidth',lin_wid)
hold on
plot(rmsd_tl2(1:LevMax-LevMin+1,1),p_ref(LevMin:LevMax),'r','LineWidth',lin_wid)
plot(rmsd_tl3(1:LevMax-LevMin+1,1),p_ref(LevMin:LevMax),'g','LineWidth',lin_wid)
set(gca,'YDir','reverse')
ylim([p_ref(LevMin+1) p_ref(LevMax)])
title('u','FontSize',fontsize1,'FontName','TimesNewRoman')
ylabel('Height (hPa)','FontSize',fontsize1,'FontName','TimesNewRoman')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

subplot(1,4,2)
plot(rmsd_tl1(1:LevMax-LevMin+1,2),p_ref(LevMin:LevMax),'b','LineWidth',lin_wid)
hold on
plot(rmsd_tl2(1:LevMax-LevMin+1,2),p_ref(LevMin:LevMax),'r','LineWidth',lin_wid)
plot(rmsd_tl3(1:LevMax-LevMin+1,2),p_ref(LevMin:LevMax),'g','LineWidth',lin_wid)
set(gca,'YDir','reverse')
ylim([p_ref(LevMin+1) p_ref(LevMax)])
title('v','FontSize',fontsize1,'FontName','TimesNewRoman')
set(gca,'YTickLabel',[])
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

subplot(1,4,3)
plot(rmsd_tl1(1:LevMax-LevMin+1,3),p_ref(LevMin:LevMax),'b','LineWidth',lin_wid)
hold on
plot(rmsd_tl2(1:LevMax-LevMin+1,3),p_ref(LevMin:LevMax),'r','LineWidth',lin_wid)
plot(rmsd_tl3(1:LevMax-LevMin+1,3),p_ref(LevMin:LevMax),'g','LineWidth',lin_wid)
set(gca,'YDir','reverse')
ylim([p_ref(LevMin+1) p_ref(LevMax)])
title('T_v','FontSize',fontsize1,'FontName','TimesNewRoman')
set(gca,'YTickLabel',[])
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

subplot(1,4,4)
plot(rmsd_tl1(1:LevMax-LevMin+1,4),p_ref(LevMin:LevMax),'b','LineWidth',lin_wid)
hold on
plot(rmsd_tl2(1:LevMax-LevMin+1,4),p_ref(LevMin:LevMax),'r','LineWidth',lin_wid)
plot(rmsd_tl3(1:LevMax-LevMin+1,4),p_ref(LevMin:LevMax),'g','LineWidth',lin_wid)
set(gca,'YDir','reverse')
ylim([p_ref(LevMin+1) p_ref(LevMax)])
title('q','FontSize',fontsize1,'FontName','TimesNewRoman')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
ylabel('Height (hPa)','FontSize',fontsize1,'FontName','TimesNewRoman')
set(gca,'YAxisLocation','right')

legend('Control','New Dyn','New Dyn + Shapiro')



% figure
% set(gcf,'position',[3 343 1276 576])
% 
% subplot(1,4,1)
% plot(rmsd_tl1(1:LevMax-LevMin+1,1),LevMin:LevMax,'b','LineWidth',lin_wid)
% hold on
% plot(rmsd_tl2(1:LevMax-LevMin+1,1),LevMin:LevMax,'r','LineWidth',lin_wid)
% plot(rmsd_tl3(1:LevMax-LevMin+1,1),LevMin:LevMax,'g','LineWidth',lin_wid)
% set(gca,'YDir','reverse')
% ylim([LevMin LevMax])
% title('u (ms^{-1})','FontSize',fontsize1,'FontName','TimesNewRoman')
% ylabel('Model Level','FontSize',fontsize1,'FontName','TimesNewRoman')
% set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
% 
% subplot(1,4,2)
% plot(rmsd_tl1(1:LevMax-LevMin+1,2),LevMin:LevMax,'b','LineWidth',lin_wid)
% hold on
% plot(rmsd_tl2(1:LevMax-LevMin+1,2),LevMin:LevMax,'r','LineWidth',lin_wid)
% plot(rmsd_tl3(1:LevMax-LevMin+1,2),LevMin:LevMax,'g','LineWidth',lin_wid)
% set(gca,'YDir','reverse')
% ylim([LevMin LevMax])
% title('v (ms^{-1})','FontSize',fontsize1,'FontName','TimesNewRoman')
% set(gca,'YTickLabel',[])
% set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
% 
% subplot(1,4,3)
% plot(rmsd_tl1(1:LevMax-LevMin+1,3),LevMin:LevMax,'b','LineWidth',lin_wid)
% hold on
% plot(rmsd_tl2(1:LevMax-LevMin+1,3),LevMin:LevMax,'r','LineWidth',lin_wid)
% plot(rmsd_tl3(1:LevMax-LevMin+1,3),LevMin:LevMax,'g','LineWidth',lin_wid)
% set(gca,'YDir','reverse')
% ylim([LevMin LevMax])
% title('T_v (K)','FontSize',fontsize1,'FontName','TimesNewRoman')
% set(gca,'YTickLabel',[])
% set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
% 
% subplot(1,4,4)
% plot(rmsd_tl1(1:LevMax-LevMin+1,4),LevMin:LevMax,'b','LineWidth',lin_wid)
% hold on
% plot(rmsd_tl2(1:LevMax-LevMin+1,4),LevMin:LevMax,'r','LineWidth',lin_wid)
% plot(rmsd_tl3(1:LevMax-LevMin+1,4),LevMin:LevMax,'g','LineWidth',lin_wid)
% set(gca,'YDir','reverse')
% ylim([LevMin LevMax])
% title('q (kgkg^{-1})','FontSize',fontsize1,'FontName','TimesNewRoman')
% set(gca,'YTickLabel',[])
% set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
% 
% legend('Control','New Dyn','New Dyn + Shapiro')
