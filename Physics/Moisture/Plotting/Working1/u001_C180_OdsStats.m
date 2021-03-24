close all
clear
clc

% Pick starting and ending date.
startdate = '20130101';
enddate   = '20130128';

% Names of forecast experiment.
expidf = 'u001_C180_DRIVER' ;

% Pick forecast lead time
leadtime = 72;

% Names of adjoint experiments.
expid1 = ['u001_C180_PH0_',num2str(leadtime)];
expid2 = ['u001_C180_PH1_',num2str(leadtime)];
expid3 = ['u001_C180_PH2_',num2str(leadtime)];
expid4 = ['u001_C180_PH2w_',num2str(leadtime)];
expid5 = ['u001_C180_PH0_GQ_',num2str(leadtime)];
expid6 = ['u001_C180_PH1_GQ_',num2str(leadtime)];
expid7 = ['u001_C180_PH2_GQ_',num2str(leadtime)];
expid8 = ['u001_C180_PH2w_GQ_',num2str(leadtime)];

% Number of days included.
numdays = datenum(enddate, 'yyyymmdd') - datenum(startdate, 'yyyymmdd') + 1;

% Initiialize the norms
dry_norm_15z = zeros(numdays,1);
dry_norm_21z = zeros(numdays,1);
wet_norm_15z = zeros(numdays,1);
wet_norm_21z = zeros(numdays,1);

for i = 1:numdays
    
    % Set analysis, verification and restart times.
    datea = datestr(datenum(startdate, 'yyyymmdd') + (i-1), 'yyyymmdd');
    datev = datestr(datenum(datea, 'yyyymmdd') + (leadtime/24), 'yyyymmdd');
    dater15 = datestr(datenum(datea, 'yyyymmdd') - 9/24, 'yyyymmdd');
    dater21 = datestr(datenum(datea, 'yyyymmdd') - 3/25, 'yyyymmdd');

    
    %Read in 15z forecast error norm
    path = ['/archive/u/drholdaw/',expidf,'/prog/Y',dater15(1:4),'/M',dater15(5:6),'/D',dater15(7:8),'/H15'];
    cd (path)
    filename = [expidf,'.Jnormf_twe.',dater15(1:4),dater15(5:6),dater15(7:8),'_15z+',datev(1:4),datev(5:6),datev(7:8),'_00z.txt'];
    fid = fopen(filename,'r');
    D_15 = textscan(fid,'%s %f %s %f %s %s %s \n %s %s %s %s %s %s %s %s \n %f %f %f %f %f %f \n %f %f ');
    fclose(fid);

    dry_norm_15z(i) = D_15{1,16} - D_15{1,21};
    wet_norm_15z(i) = D_15{1,16};
        
    %Read in 21z forecast error norm
    path = ['/archive/u/drholdaw/',expidf,'/prog/Y',dater21(1:4),'/M',dater21(5:6),'/D',dater21(7:8),'/H21'];
    cd (path)
    filename = [expidf,'.Jnormf_twe.',dater21(1:4),dater21(5:6),dater21(7:8),'_21z+',datev(1:4),datev(5:6),datev(7:8),'_00z.txt'];
    fid = fopen(filename,'r');
    D_21 = textscan(fid,'%s %f %s %f %s %s %s \n %s %s %s %s %s %s %s %s \n %f %f %f %f %f %f \n %f %f ');
    fclose(fid);

    dry_norm_21z(i) = D_21{1,16} - D_21{1,21};
    wet_norm_21z(i) = D_21{1,16};

    clear D_21 D_15
                    
    cd /home/drholdaw/Lin_Moist_Physics/DAS_plotting/Cloud' Fraction'/
    
end

dry_norm_nl_error = dry_norm_21z - dry_norm_15z;
wet_norm_nl_error = wet_norm_21z - wet_norm_15z;


% Read in linear model errors from odsstats
dir = ['/discover/nobackup/drholdaw/u001_C180_DRIVER/DRIVER/OdsStats/',expid1,'/'];
cd(dir)
fid = fopen('accum.sum_odsstats.txt');
C = textscan(fid,'%f %f %s %f %f');
fclose(fid);
u001_C180_PH0_error = C{5}; clear C
u001_C180_PH0_error = u001_C180_PH0_error(1:numdays);

dir = ['/discover/nobackup/drholdaw/u001_C180_DRIVER/DRIVER/OdsStats/',expid2,'/'];
cd(dir)
fid = fopen('accum.sum_odsstats.txt');
C = textscan(fid,'%f %f %s %f %f');
fclose(fid);
u001_C180_PH1_error = C{5}; clear C
u001_C180_PH1_error = u001_C180_PH1_error(1:numdays);

dir = ['/discover/nobackup/drholdaw/u001_C180_DRIVER/DRIVER/OdsStats/',expid3,'/'];
cd(dir)
fid = fopen('accum.sum_odsstats.txt');
C = textscan(fid,'%f %f %s %f %f');
fclose(fid);
u001_C180_PH2_error = C{5}; clear C
u001_C180_PH2_error = u001_C180_PH2_error(1:numdays);

dir = ['/discover/nobackup/drholdaw/u001_C180_DRIVER/DRIVER/OdsStats/',expid4,'/'];
cd(dir)
fid = fopen('accum.sum_odsstats.txt');
C = textscan(fid,'%f %f %s %f %f');
fclose(fid);
u001_C180_PH2w_error = C{5}; clear C
u001_C180_PH2w_error = u001_C180_PH2w_error(1:numdays);

dir = ['/discover/nobackup/drholdaw/u001_C180_DRIVER/DRIVER/OdsStats/',expid5,'/'];
cd(dir)
fid = fopen('accum.sum_odsstats.txt');
C = textscan(fid,'%f %f %s %f %f');
fclose(fid);
u001_C180_PH0_GQ_error = C{5}; clear C
u001_C180_PH0_GQ_error = u001_C180_PH0_GQ_error(1:numdays);

dir = ['/discover/nobackup/drholdaw/u001_C180_DRIVER/DRIVER/OdsStats/',expid6,'/'];
cd(dir)
fid = fopen('accum.sum_odsstats.txt');
C = textscan(fid,'%f %f %s %f %f');
fclose(fid);
u001_C180_PH1_GQ_error = C{5}; clear C
u001_C180_PH1_GQ_error = u001_C180_PH1_GQ_error(1:numdays);

dir = ['/discover/nobackup/drholdaw/u001_C180_DRIVER/DRIVER/OdsStats/',expid7,'/'];
cd(dir)
fid = fopen('accum.sum_odsstats.txt');
C = textscan(fid,'%f %f %s %f %f');
fclose(fid);
u001_C180_PH2_GQ_error = C{5}; clear C
u001_C180_PH2_GQ_error = u001_C180_PH2_GQ_error(1:numdays);

dir = ['/discover/nobackup/drholdaw/u001_C180_DRIVER/DRIVER/OdsStats/',expid8,'/'];
cd(dir)
fid = fopen('accum.sum_odsstats.txt');
C = textscan(fid,'%f %f %s %f %f');
fclose(fid);
u001_C180_PH2w_GQ_error = C{5}; clear C
u001_C180_PH2w_GQ_error = u001_C180_PH2w_GQ_error(1:numdays);


cd /home/drholdaw/Lin_Moist_Physics/DAS_plotting/Cloud' Fraction'/

fprintf('Percentage PH0    = %g \n',    100*mean(u001_C180_PH0_error./dry_norm_nl_error))
fprintf('Percentage PH0_GQ = %g \n\n',    100*mean(u001_C180_PH0_GQ_error./dry_norm_nl_error))

fprintf('Percentage PH1    = %g \n',    100*mean(u001_C180_PH1_error./dry_norm_nl_error))
fprintf('Percentage PH1_GQ = %g \n\n',    100*mean(u001_C180_PH1_GQ_error./dry_norm_nl_error))

fprintf('Percentage PH2    = %g \n',    100*mean(u001_C180_PH2_error./dry_norm_nl_error))
fprintf('Percentage PH2_GQ = %g \n\n',    100*mean(u001_C180_PH2_GQ_error./dry_norm_nl_error))

fprintf('Percentage PH2w    = %g \n',    100*mean(u001_C180_PH2w_error./wet_norm_nl_error))
fprintf('Percentage PH2w_GQ = %g \n\n',    100*mean(u001_C180_PH2w_GQ_error./wet_norm_nl_error))

s = std(u001_C180_PH0_error)
s = std(u001_C180_PH0_GQ_error)

line_wid = 1.0;
line_wid1 = 1.5;
fontsize = 16;


figure
set(gcf,'position',[1281 1 1280 943])

subplot(2,2,1)
plot(dry_norm_15z,'k--','LineWidth', line_wid1)
hold on
plot(dry_norm_21z,'k-', 'LineWidth', line_wid1)
plot(dry_norm_nl_error,'k', 'LineWidth', line_wid)
plot(u001_C180_PH0_error,'b', 'LineWidth', line_wid)
plot(u001_C180_PH0_GQ_error,'r', 'LineWidth', line_wid)

hline = refline([0 0]);
set(hline,'Color','k', 'LineWidth', line_wid)
title('Dry Physics - OLD BL','FontSize',fontsize,'FontName','TimesNewRoman')
ylabel('Energy (Jkg^{-1})','FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
xlim([1 numdays])

subplot(2,2,2)
plot(dry_norm_15z,'k--','LineWidth', line_wid1)
hold on
plot(dry_norm_21z,'k-', 'LineWidth', line_wid1)
plot(dry_norm_nl_error,'k', 'LineWidth', line_wid)
plot(u001_C180_PH1_error,'b', 'LineWidth', line_wid)
plot(u001_C180_PH1_GQ_error,'r', 'LineWidth', line_wid)

hline = refline([0 0]);
set(hline,'Color','k', 'LineWidth', line_wid)
title('Dry Physics - New BL','FontSize',fontsize,'FontName','TimesNewRoman')
ylabel('Energy (Jkg^{-1})','FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
xlim([1 numdays])

subplot(2,2,3)
plot(dry_norm_15z,'k--','LineWidth', line_wid1)
hold on
plot(dry_norm_21z,'k-', 'LineWidth', line_wid1)
plot(dry_norm_nl_error,'k', 'LineWidth', line_wid)
plot(u001_C180_PH2_error,'b', 'LineWidth', line_wid)
plot(u001_C180_PH2_GQ_error,'r', 'LineWidth', line_wid)

hline = refline([0 0]);
set(hline,'Color','k', 'LineWidth', line_wid)
title('Moist Physics - Dry Norm','FontSize',fontsize,'FontName','TimesNewRoman')
ylabel('Energy (Jkg^{-1})','FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
xlim([1 numdays])

subplot(2,2,4)
plot(wet_norm_15z,'k--','LineWidth', line_wid1)
hold on
plot(wet_norm_21z,'k-', 'LineWidth', line_wid1)
plot(wet_norm_nl_error,'k', 'LineWidth', line_wid)
plot(u001_C180_PH2w_error,'b', 'LineWidth', line_wid)
plot(u001_C180_PH2w_GQ_error,'r', 'LineWidth', line_wid)

hline = refline([0 0]);
set(hline,'Color','k', 'LineWidth', line_wid)
title('Moist Physics - Moist Norm','FontSize',fontsize,'FontName','TimesNewRoman')
ylabel('Energy (Jkg^{-1})','FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
xlim([1 numdays])


% ttest_PH0 = ttest2(u001_C180_PH0_error, u001_C180_PH0_GQ_error)
% ttest_PH1 = ttest2(u001_C180_PH0_error, u001_C180_PH0_GQ_error)
% ttest_PH2 = ttest2(u001_C180_PH0_error, u001_C180_PH0_GQ_error)
% ttest_PH2w = ttest2(u001_C180_PH0_error, u001_C180_PH0_GQ_error)

