close all
clear
clc
mydir = pwd;

% dir = 'wrk.sandy';
% dir = 'wrk.dust';
dir = 'btmp.25956';

% exp = 'sens_sandy';
% exp = 'C360_513_snamma';
exp = 'v000_C180';


%Set directories of experiment
dir_free = dir;
dir_repl = dir;
dir_tlm1 = dir;
dir_tlm2 = dir;
dir_tlm3 = dir;

% Set experiment names
exp_free = exp;
exp_repl = exp;
exp_tlm1 = exp;
exp_tlm2 = exp;
exp_tlm3 = exp;

%Choose physics options for TLM
phy_tlm1 = '22000';
phy_tlm2 = '2Tk000';
phy_tlm3 = '2Tl000';

%Choose shapiro filter options for TLM
shp_tlm1 = '1';
shp_tlm2 = '1';
shp_tlm3 = '1';

%Choose 'other' options for TLM/count2

opt_tlm1 = '';
opt_tlm2 = '';
opt_tlm3 = '';

%Choose highest model level
LevMin = 1;

%Choose start date
% datea = '20060827';
% datea = '20060827';
datea = '20140201';

%Choose start time
timea = '030000';

%Choose length of forecasts
lead = 3;

%Plotting options
lin_wid = 1.75;
fontsize = 11;
fontsize1 = 12;

%Choose Region,
% 1 = Global
% 2 = Tropics 23S to 23N, below 100hPa
% 3 = Northern hemisphere
% 4 = Southern hemisphere
region = 1;

if region == 1 
    fprintf(' \nWill compute RMSD for the whole globe \n\n')
elseif region == 2
    fprintf('Will compute RMSD for the tropics \n\n')  
elseif region == 3
    fprintf('Will compute RMSD for the northern hemisphere \n\n')
elseif region == 4
    fprintf('Will compute RMSD for the southern hemisphere \n\n')
else
    fprintf('Not a valid region selection \n\n')
end

cd /home/drholdaw/LinearisedPhysics/Inputs/
pref
p = 0.5*(p_ref(1:end-1) + p_ref(2:end));

date_start = datenum(str2double(datea(1:4)), str2double(datea(5:6)), str2double(datea(7:8)), ...
                     str2double(timea(1:2)), str2double(timea(3:4)), str2double(timea(5:6)));

datetimef = datestr(date_start + lead/24     , 'yyyymmddHHMMSS');
datetimev = datestr(date_start + 1           , 'yyyymmddHHMMSS'); 

fprintf(' Reading in the states \n')

%Load Free (background) State.
dir = ['/discover/nobackup/drholdaw/',dir_free,'/prog/prog_free_nosinks/']; cd(dir)
file = [exp_free,'.prog.eta.',datetimef(1:8),'_',datetimef(9:12),'z.nc4'];

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
dir = ['/discover/nobackup/drholdaw/',dir_repl,'/prog/prog_replay_nosinks/']; cd(dir)
file = [exp_repl,'.prog.eta.',datetimef(1:8),'_',datetimef(9:12),'z.nc4'];

u_replay = ncread(file,'u');
v_replay = ncread(file,'v');
t_replay = ncread(file,'tv');
q_replay = ncread(file,'sphu');
p_replay = ncread(file,'delp');
qi_replay = ncread(file,'qitot');
ql_replay = ncread(file,'qltot');
o3_replay = ncread(file,'ozone');

%Load first Tlm state to compare.
dir = ['/discover/nobackup/drholdaw/',dir_tlm1,'/sens.',datetimev(1:8),'.000000']; cd(dir)
file = [exp_tlm1,'.fvpert.eta.',datetimef(1:8),'_',datetimef(9:12),'z_ST03z_P',phy_tlm1,'_S',shp_tlm1,opt_tlm1,'.nc4'];

u_tl1 = ncread(file,'U');
v_tl1 = ncread(file,'V');
t_tl1 = ncread(file,'TV');
q_tl1 = ncread(file,'QV');
p_tl1 = ncread(file,'DP');
qi_tl1 = ncread(file,'QI');
ql_tl1 = ncread(file,'QL');
o3_tl1 = ncread(file,'O3');

%Load second Tlm state to compare.
dir = ['/discover/nobackup/drholdaw/',dir_tlm2,'/sens.',datetimev(1:8),'.000000']; cd(dir)
file = [exp_tlm2,'.fvpert.eta.',datetimef(1:8),'_',datetimef(9:12),'z_ST03z_P',phy_tlm2,'_S',shp_tlm2,opt_tlm2,'.nc4'];

u_tl2 = ncread(file,'U');
v_tl2 = ncread(file,'V');
t_tl2 = ncread(file,'TV');
q_tl2 = ncread(file,'QV');
p_tl2 = ncread(file,'DP');
qi_tl2 = ncread(file,'QI');
ql_tl2 = ncread(file,'QL');
o3_tl2 = ncread(file,'O3');

%Load third Tlm state to compare.
dir = ['/discover/nobackup/drholdaw/',dir_tlm3,'/sens.',datetimev(1:8),'.000000']; cd(dir)
file = [exp_tlm3,'.fvpert.eta.',datetimef(1:8),'_',datetimef(9:12),'z_ST03z_P',phy_tlm3,'_S',shp_tlm3,opt_tlm3,'.nc4'];

u_tl3 = ncread(file,'U');
v_tl3 = ncread(file,'V');
t_tl3 = ncread(file,'TV');
q_tl3 = ncread(file,'QV');
p_tl3 = ncread(file,'DP');
qi_tl3 = ncread(file,'QI');
ql_tl3 = ncread(file,'QL');
o3_tl3 = ncread(file,'O3');

cd(mydir)

fprintf(' Done reading in the states \n\n')

fprintf('  RMSD computation \n')

im = size(u_free,1);
jm = size(u_free,2);
lm = size(u_free,3);

if region == 1
    LonMin = 1;
    LonMax = im;
    LatMin = 1;
    LatMax = jm;
elseif region == 2
    LonMin = 1;
    LonMax = im;
    LatMin = find(lat == -24);
    LatMax = find(lat ==  24);
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

%Compute NL perturbation trajectory.
u_nlm = u_replay - u_free;
v_nlm = v_replay - v_free;
t_nlm = t_replay - t_free;
q_nlm = q_replay - q_free;
p_nlm = p_replay - p_free;
qi_nlm = qi_replay - qi_free;
ql_nlm = ql_replay - ql_free;
o3_nlm = o3_replay - o3_free;

rmsd_tl1 = zeros(lm,8);
rmsd_tl2 = zeros(lm,8);
rmsd_tl3 = zeros(lm,8);

n = (LonMax-LonMin+1)*(LatMax-LatMin+1);

%Loop over model all levels.
for k = 1:lm

    rmsd_tl1(k,1) = sqrt(  sum(sum( (u_nlm(LonMin:LonMax,LatMin:LatMax,k) - u_tl1(LonMin:LonMax,LatMin:LatMax,k)).^2 )) / n );
    rmsd_tl2(k,1) = sqrt(  sum(sum( (u_nlm(LonMin:LonMax,LatMin:LatMax,k) - u_tl2(LonMin:LonMax,LatMin:LatMax,k)).^2 )) / n );
    rmsd_tl3(k,1) = sqrt(  sum(sum( (u_nlm(LonMin:LonMax,LatMin:LatMax,k) - u_tl3(LonMin:LonMax,LatMin:LatMax,k)).^2 )) / n );

    rmsd_tl1(k,2) = sqrt(  sum(sum( (v_nlm(LonMin:LonMax,LatMin:LatMax,k) - v_tl1(LonMin:LonMax,LatMin:LatMax,k)).^2 )) / n );
    rmsd_tl2(k,2) = sqrt(  sum(sum( (v_nlm(LonMin:LonMax,LatMin:LatMax,k) - v_tl2(LonMin:LonMax,LatMin:LatMax,k)).^2 )) / n );
    rmsd_tl3(k,2) = sqrt(  sum(sum( (v_nlm(LonMin:LonMax,LatMin:LatMax,k) - v_tl3(LonMin:LonMax,LatMin:LatMax,k)).^2 )) / n );

    rmsd_tl1(k,3) = sqrt(  sum(sum( (t_nlm(LonMin:LonMax,LatMin:LatMax,k) - t_tl1(LonMin:LonMax,LatMin:LatMax,k)).^2 )) / n );
    rmsd_tl2(k,3) = sqrt(  sum(sum( (t_nlm(LonMin:LonMax,LatMin:LatMax,k) - t_tl2(LonMin:LonMax,LatMin:LatMax,k)).^2 )) / n );
    rmsd_tl3(k,3) = sqrt(  sum(sum( (t_nlm(LonMin:LonMax,LatMin:LatMax,k) - t_tl3(LonMin:LonMax,LatMin:LatMax,k)).^2 )) / n );

    rmsd_tl1(k,4) = sqrt(  sum(sum( (q_nlm(LonMin:LonMax,LatMin:LatMax,k) - q_tl1(LonMin:LonMax,LatMin:LatMax,k)).^2 )) / n );
    rmsd_tl2(k,4) = sqrt(  sum(sum( (q_nlm(LonMin:LonMax,LatMin:LatMax,k) - q_tl2(LonMin:LonMax,LatMin:LatMax,k)).^2 )) / n );
    rmsd_tl3(k,4) = sqrt(  sum(sum( (q_nlm(LonMin:LonMax,LatMin:LatMax,k) - q_tl3(LonMin:LonMax,LatMin:LatMax,k)).^2 )) / n );

    rmsd_tl1(k,5) = sqrt(  sum(sum( (p_nlm(LonMin:LonMax,LatMin:LatMax,k) - p_tl1(LonMin:LonMax,LatMin:LatMax,k)).^2 )) / n );
    rmsd_tl2(k,5) = sqrt(  sum(sum( (p_nlm(LonMin:LonMax,LatMin:LatMax,k) - p_tl2(LonMin:LonMax,LatMin:LatMax,k)).^2 )) / n );
    rmsd_tl3(k,5) = sqrt(  sum(sum( (p_nlm(LonMin:LonMax,LatMin:LatMax,k) - p_tl3(LonMin:LonMax,LatMin:LatMax,k)).^2 )) / n );

    rmsd_tl1(k,6) = sqrt(  sum(sum( (qi_nlm(LonMin:LonMax,LatMin:LatMax,k) - qi_tl1(LonMin:LonMax,LatMin:LatMax,k)).^2 )) / n );
    rmsd_tl2(k,6) = sqrt(  sum(sum( (qi_nlm(LonMin:LonMax,LatMin:LatMax,k) - qi_tl2(LonMin:LonMax,LatMin:LatMax,k)).^2 )) / n );
    rmsd_tl3(k,6) = sqrt(  sum(sum( (qi_nlm(LonMin:LonMax,LatMin:LatMax,k) - qi_tl3(LonMin:LonMax,LatMin:LatMax,k)).^2 )) / n );

    rmsd_tl1(k,7) = sqrt(  sum(sum( (ql_nlm(LonMin:LonMax,LatMin:LatMax,k) - ql_tl1(LonMin:LonMax,LatMin:LatMax,k)).^2 )) / n );
    rmsd_tl2(k,7) = sqrt(  sum(sum( (ql_nlm(LonMin:LonMax,LatMin:LatMax,k) - ql_tl2(LonMin:LonMax,LatMin:LatMax,k)).^2 )) / n );
    rmsd_tl3(k,7) = sqrt(  sum(sum( (ql_nlm(LonMin:LonMax,LatMin:LatMax,k) - ql_tl3(LonMin:LonMax,LatMin:LatMax,k)).^2 )) / n );

    rmsd_tl1(k,8) = sqrt(  sum(sum( (o3_nlm(LonMin:LonMax,LatMin:LatMax,k) - o3_tl1(LonMin:LonMax,LatMin:LatMax,k)).^2 )) / n );
    rmsd_tl2(k,8) = sqrt(  sum(sum( (o3_nlm(LonMin:LonMax,LatMin:LatMax,k) - o3_tl2(LonMin:LonMax,LatMin:LatMax,k)).^2 )) / n );
    rmsd_tl3(k,8) = sqrt(  sum(sum( (o3_nlm(LonMin:LonMax,LatMin:LatMax,k) - o3_tl3(LonMin:LonMax,LatMin:LatMax,k)).^2 )) / n );
end

fprintf('  Done RMSD computation \n\n')



figure
set(gcf,'position',[3 343 1276 576])

subplot(1,8,1)
plot(rmsd_tl1(1:lm-LevMin+1,1),p(LevMin:lm),'b','LineWidth',lin_wid)
hold on
plot(rmsd_tl2(1:lm-LevMin+1,1),p(LevMin:lm),'r','LineWidth',lin_wid)
plot(rmsd_tl3(1:lm-LevMin+1,1),p(LevMin:lm),'g--','LineWidth',lin_wid)
set(gca,'YDir','reverse')
ylim([p(LevMin) p(lm)])
title('u (ms^{-1})','FontSize',fontsize1,'FontName','TimesNewRoman')
ylabel('Height (hPa)','FontSize',fontsize1,'FontName','TimesNewRoman')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

subplot(1,8,2)
plot(rmsd_tl1(1:lm-LevMin+1,2),p(LevMin:lm),'b','LineWidth',lin_wid)
hold on
plot(rmsd_tl2(1:lm-LevMin+1,2),p(LevMin:lm),'r','LineWidth',lin_wid)
plot(rmsd_tl3(1:lm-LevMin+1,2),p(LevMin:lm),'g--','LineWidth',lin_wid)
set(gca,'YDir','reverse')
ylim([p(LevMin) p(lm)])
title('v (ms^{-1})','FontSize',fontsize1,'FontName','TimesNewRoman')
set(gca,'YTickLabel',[])
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

subplot(1,8,3)
plot(rmsd_tl1(1:lm-LevMin+1,3),p(LevMin:lm),'b','LineWidth',lin_wid)
hold on
plot(rmsd_tl2(1:lm-LevMin+1,3),p(LevMin:lm),'r','LineWidth',lin_wid)
plot(rmsd_tl3(1:lm-LevMin+1,3),p(LevMin:lm),'g--','LineWidth',lin_wid)
set(gca,'YDir','reverse')
ylim([p(LevMin) p(lm)])
title('T_v (K)','FontSize',fontsize1,'FontName','TimesNewRoman')
set(gca,'YTickLabel',[])
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

subplot(1,8,4)
plot(rmsd_tl1(1:lm-LevMin+1,4),p(LevMin:lm),'b','LineWidth',lin_wid)
hold on
plot(rmsd_tl2(1:lm-LevMin+1,4),p(LevMin:lm),'r','LineWidth',lin_wid)
plot(rmsd_tl3(1:lm-LevMin+1,4),p(LevMin:lm),'g--','LineWidth',lin_wid)
set(gca,'YDir','reverse')
ylim([p(LevMin) p(lm)])
title('q (kgkg^{-1})','FontSize',fontsize1,'FontName','TimesNewRoman')
set(gca,'YTickLabel',[])
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

subplot(1,8,5)
plot(rmsd_tl1(1:lm-LevMin+1,5),p(LevMin:lm),'b','LineWidth',lin_wid)
hold on
plot(rmsd_tl2(1:lm-LevMin+1,5),p(LevMin:lm),'r','LineWidth',lin_wid)
plot(rmsd_tl3(1:lm-LevMin+1,5),p(LevMin:lm),'g--','LineWidth',lin_wid)
set(gca,'YDir','reverse')
ylim([p(LevMin) p(lm)])
title('p (Pa)','FontSize',fontsize1,'FontName','TimesNewRoman')
set(gca,'YTickLabel',[])
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

subplot(1,8,6)
plot(rmsd_tl1(1:lm-LevMin+1,6),p(LevMin:lm),'b','LineWidth',lin_wid)
hold on
plot(rmsd_tl2(1:lm-LevMin+1,6),p(LevMin:lm),'r','LineWidth',lin_wid)
plot(rmsd_tl3(1:lm-LevMin+1,6),p(LevMin:lm),'g--','LineWidth',lin_wid)
set(gca,'YDir','reverse')
ylim([p(LevMin) p(lm)])
title('q_i (kgkg^{-1})','FontSize',fontsize1,'FontName','TimesNewRoman')
set(gca,'YTickLabel',[])
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

subplot(1,8,7)
plot(rmsd_tl1(1:lm-LevMin+1,7),p(LevMin:lm),'b','LineWidth',lin_wid)
hold on
plot(rmsd_tl2(1:lm-LevMin+1,7),p(LevMin:lm),'r','LineWidth',lin_wid)
plot(rmsd_tl3(1:lm-LevMin+1,7),p(LevMin:lm),'g--','LineWidth',lin_wid)
set(gca,'YDir','reverse')
ylim([p(LevMin) p(lm)])
title('q_l (kgkg^{-1})','FontSize',fontsize1,'FontName','TimesNewRoman')
set(gca,'YTickLabel',[])
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

subplot(1,8,8)
plot(rmsd_tl1(1:lm-LevMin+1,8),p(LevMin:lm),'b','LineWidth',lin_wid)
hold on
plot(rmsd_tl2(1:lm-LevMin+1,8),p(LevMin:lm),'r','LineWidth',lin_wid)
plot(rmsd_tl3(1:lm-LevMin+1,8),p(LevMin:lm),'g--','LineWidth',lin_wid)
set(gca,'YDir','reverse')
ylim([p(0+1) p(lm)])
title('o3 (ppm)','FontSize',fontsize1,'FontName','TimesNewRoman')
ylabel('Height (hPa)','FontSize',fontsize1,'FontName','TimesNewRoman')

set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'YAxisLocation','right')





figure
set(gcf,'position',[3 343 1276 576])

subplot(1,8,1)
plot(rmsd_tl1(1:lm-LevMin+1,1),LevMin:lm,'b','LineWidth',lin_wid)
hold on
plot(rmsd_tl2(1:lm-LevMin+1,1),LevMin:lm,'r','LineWidth',lin_wid)
plot(rmsd_tl3(1:lm-LevMin+1,1),LevMin:lm,'g--','LineWidth',lin_wid)
set(gca,'YDir','reverse')
ylim([LevMin lm])
title('u (ms^{-1})','FontSize',fontsize1,'FontName','TimesNewRoman')
ylabel('Model level','FontSize',fontsize1,'FontName','TimesNewRoman')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

subplot(1,8,2)
plot(rmsd_tl1(1:lm-LevMin+1,2),LevMin:lm,'b','LineWidth',lin_wid)
hold on
plot(rmsd_tl2(1:lm-LevMin+1,2),LevMin:lm,'r','LineWidth',lin_wid)
plot(rmsd_tl3(1:lm-LevMin+1,2),LevMin:lm,'g--','LineWidth',lin_wid)
set(gca,'YDir','reverse')
ylim([LevMin lm])
title('v (ms^{-1})','FontSize',fontsize1,'FontName','TimesNewRoman')
set(gca,'YTickLabel',[])
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

subplot(1,8,3)
plot(rmsd_tl1(1:lm-LevMin+1,3),LevMin:lm,'b','LineWidth',lin_wid)
hold on
plot(rmsd_tl2(1:lm-LevMin+1,3),LevMin:lm,'r','LineWidth',lin_wid)
plot(rmsd_tl3(1:lm-LevMin+1,3),LevMin:lm,'g--','LineWidth',lin_wid)
set(gca,'YDir','reverse')
ylim([LevMin lm])
title('T_v (K)','FontSize',fontsize1,'FontName','TimesNewRoman')
set(gca,'YTickLabel',[])
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

subplot(1,8,4)
plot(rmsd_tl1(1:lm-LevMin+1,4),LevMin:lm,'b','LineWidth',lin_wid)
hold on
plot(rmsd_tl2(1:lm-LevMin+1,4),LevMin:lm,'r','LineWidth',lin_wid)
plot(rmsd_tl3(1:lm-LevMin+1,4),LevMin:lm,'g--','LineWidth',lin_wid)
set(gca,'YDir','reverse')
ylim([LevMin lm])
title('q (kgkg^{-1})','FontSize',fontsize1,'FontName','TimesNewRoman')
set(gca,'YTickLabel',[])
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

subplot(1,8,5)
plot(rmsd_tl1(1:lm-LevMin+1,5),LevMin:lm,'b','LineWidth',lin_wid)
hold on
plot(rmsd_tl2(1:lm-LevMin+1,5),LevMin:lm,'r','LineWidth',lin_wid)
plot(rmsd_tl3(1:lm-LevMin+1,5),LevMin:lm,'g--','LineWidth',lin_wid)
set(gca,'YDir','reverse')
ylim([LevMin lm])
title('p (Pa)','FontSize',fontsize1,'FontName','TimesNewRoman')
set(gca,'YTickLabel',[])
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

subplot(1,8,6)
plot(rmsd_tl1(1:lm-LevMin+1,6),LevMin:lm,'b','LineWidth',lin_wid)
hold on
plot(rmsd_tl2(1:lm-LevMin+1,6),LevMin:lm,'r','LineWidth',lin_wid)
plot(rmsd_tl3(1:lm-LevMin+1,6),LevMin:lm,'g--','LineWidth',lin_wid)
set(gca,'YDir','reverse')
ylim([LevMin lm])
title('q_i (kgkg^{-1})','FontSize',fontsize1,'FontName','TimesNewRoman')
set(gca,'YTickLabel',[])
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

subplot(1,8,7)
plot(rmsd_tl1(1:lm-LevMin+1,7),LevMin:lm,'b','LineWidth',lin_wid)
hold on
plot(rmsd_tl2(1:lm-LevMin+1,7),LevMin:lm,'r','LineWidth',lin_wid)
plot(rmsd_tl3(1:lm-LevMin+1,7),LevMin:lm,'g--','LineWidth',lin_wid)
set(gca,'YDir','reverse')
ylim([LevMin lm])
title('q_l (kgkg^{-1})','FontSize',fontsize1,'FontName','TimesNewRoman')
set(gca,'YTickLabel',[])
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

subplot(1,8,8)
plot(rmsd_tl1(1:lm-LevMin+1,8),LevMin:lm,'b','LineWidth',lin_wid)
hold on
plot(rmsd_tl2(1:lm-LevMin+1,8),LevMin:lm,'r','LineWidth',lin_wid)
plot(rmsd_tl3(1:lm-LevMin+1,8),LevMin:lm,'g--','LineWidth',lin_wid)
set(gca,'YDir','reverse')
ylim([LevMin lm])
title('o3 (ppm)','FontSize',fontsize1,'FontName','TimesNewRoman')
ylabel('Model level','FontSize',fontsize1,'FontName','TimesNewRoman')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'YAxisLocation','right')

figure(1)
