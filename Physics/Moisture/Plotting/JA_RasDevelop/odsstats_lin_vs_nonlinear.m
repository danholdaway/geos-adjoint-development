close all
clear
clc

startdate = '20120316';
enddate   = '20120416';

expid1 = 'iau590a';
expid2 = 'iau590a_mp';

daterange = (datenum(startdate, 'yyyymmdd') : 1 : datenum(enddate, 'yyyymmdd') )';

numdays = length(daterange);

%Add two days to get end point of 24h forecast.
dates = [datestr(daterange, 'yyyymmdd') ; datestr(daterange(end)+1, 'yyyymmdd') ; datestr(daterange(end)+2, 'yyyymmdd')];

norm_dp_dn_15z = zeros(numdays,1);
norm_dp_dn_21z = zeros(numdays,1);

norm_dp_mn_15z = zeros(numdays,1);
norm_dp_mn_21z = zeros(numdays,1);

norm_mp_dn_15z = zeros(numdays,1);
norm_mp_dn_21z = zeros(numdays,1);

norm_mp_mn_15z = zeros(numdays,1);
norm_mp_mn_21z = zeros(numdays,1);

for i = 1:numdays
       
    path = ['/archive/u/drholdaw/p13/',expid1,'/patch13/prog_txe/Y',dates(i,1:4),'/M',dates(i,5:6),'/D',dates(i,7:8),'/H15'];
    cd (path)
    filename = [expid1,'.Jnormf.',dates(i,1:4),dates(i,5:6),dates(i,7:8),'_15z+',dates(i+2,1:4),dates(i+2,5:6),dates(i+2,7:8),'_00z.txt'];
    fid = fopen(filename,'r');
    D_15 = textscan(fid,'%s %f %s %f %s %s %s \n %s %s %s %s %s %s \n %f %f %f %f %f %f \n %s %f %s %f %s %s %s \n %s %s %s %s %s %s \n %f %f %f %f %f %f');
    fclose(fid);

    norm_dp_dn_15z(i) = D_15{1,14};
%     norm_dp_mn_15z(i) = D_15{1,33};
    
    path = ['/archive/u/drholdaw/p13/',expid1,'/patch13/prog_txe/Y',dates(i,1:4),'/M',dates(i,5:6),'/D',dates(i,7:8),'/H21'];
    cd (path)
    filename = [expid1,'.Jnormf.',dates(i,1:4),dates(i,5:6),dates(i,7:8),'_21z+',dates(i+2,1:4),dates(i+2,5:6),dates(i+2,7:8),'_00z.txt'];
    fid = fopen(filename,'r');
    D_21 = textscan(fid,'%s %f %s %f %s %s %s \n %s %s %s %s %s %s \n %f %f %f %f %f %f \n %s %f %s %f %s %s %s \n %s %s %s %s %s %s \n %f %f %f %f %f %f');
    fclose(fid);

    norm_dp_dn_21z(i) = D_21{1,14};
%     norm_dp_mn_21z(i) = D_21{1,33};
% 
%     clear D_21 D_15
           
    cd /discover/nobackup/drholdaw/odsstats/
   
    path = ['/archive/u/drholdaw/p13/',expid2,'/patch13/prog_txe/Y',dates(i,1:4),'/M',dates(i,5:6),'/D',dates(i,7:8),'/H15'];
    cd (path)
    filename = [expid2,'.Jnormf.',dates(i,1:4),dates(i,5:6),dates(i,7:8),'_15z+',dates(i+2,1:4),dates(i+2,5:6),dates(i+2,7:8),'_00z.txt'];
    fid = fopen(filename,'r');
    D_15 = textscan(fid,'%s %f %s %f %s %s %s \n %s %s %s %s %s %s \n %f %f %f %f %f %f \n %s %f %s %f %s %s %s \n %s %s %s %s %s %s \n %f %f %f %f %f %f');
    fclose(fid);

    norm_mp_dn_15z(i) = D_15{1,14};
    norm_mp_mn_15z(i) = D_15{1,33};
    
    path = ['/archive/u/drholdaw/p13/',expid2,'/patch13/prog_txe/Y',dates(i,1:4),'/M',dates(i,5:6),'/D',dates(i,7:8),'/H21'];
    cd (path)
    filename = [expid2,'.Jnormf.',dates(i,1:4),dates(i,5:6),dates(i,7:8),'_21z+',dates(i+2,1:4),dates(i+2,5:6),dates(i+2,7:8),'_00z.txt'];
    fid = fopen(filename,'r');
    D_21 = textscan(fid,'%s %f %s %f %s %s %s \n %s %s %s %s %s %s \n %f %f %f %f %f %f \n %s %f %s %f %s %s %s \n %s %s %s %s %s %s \n %f %f %f %f %f %f');
    fclose(fid);

    norm_mp_dn_21z(i) = D_21{1,14};
    norm_mp_mn_21z(i) = D_21{1,33};

    clear D_21 D_15
           
    cd /home/drholdaw/Lin_Moist_Physics/Moist_Dev_MWR/
    
end

nonlin_error_dp_dn = norm_dp_dn_21z - norm_dp_dn_15z;
nonlin_error_dp_mn = norm_dp_mn_21z - norm_dp_mn_15z;
nonlin_error_mp_dn = norm_mp_dn_21z - norm_mp_dn_15z;
nonlin_error_mp_mn = norm_mp_mn_21z - norm_mp_mn_15z;

cd /discover/nobackup/drholdaw/ExperimentData/Journal_Articles/moist_dev_mwr/20130205

fid = fopen(['odsstats_iau590_dp_dn.txt']);
C = textscan(fid,'%f %f %s %f %f');
fclose(fid);
linear_error_dp_dn = C{5}; clear C

fid = fopen(['odsstats_iau590_dp_mn.txt']);
C = textscan(fid,'%f %f %s %f %f');
fclose(fid);
linear_error_dp_mn = C{5}; clear C

fid = fopen(['odsstats_iau590_mp_dn.txt']);
C = textscan(fid,'%f %f %s %f %f');
fclose(fid);
linear_error_mp_dn = C{5}; clear C

fid = fopen(['odsstats_iau590_mp_mn.txt']);
C = textscan(fid,'%f %f %s %f %f');
fclose(fid);
linear_error_mp_mn = C{5}; clear C

cd /home/drholdaw/Lin_Moist_Physics/Moist_Dev_MWR/

fprintf('Percentage of error Dry Pysics Dry Norm = %g \n', 100*mean(linear_error_dp_dn./nonlin_error_dp_dn))
fprintf('Percentage of error Dry Pysics Wet Norm = %g \n', 100*mean(linear_error_dp_mn./nonlin_error_mp_mn))
fprintf('Percentage of error Wet Pysics Dry Norm = %g \n', 100*mean(linear_error_mp_dn./nonlin_error_mp_dn))
fprintf('Percentage of error Wet Pysics Wet Norm = %g \n', 100*mean(linear_error_mp_mn./nonlin_error_mp_mn))

dates_plot = str2num(dates(2:end-1,end-1:end));



line_wid = 1.0;
line_wid1 = 1.5;
fontsize = 10;

months_extra = ['-Mar-12'; '-Mar-12'; '-Mar-12'; '-Apr-12'; '-Apr-12'; '-Apr-12'];

figure('Position',[125   398   439   290])

plot(norm_mp_dn_15z,'k--','LineWidth', line_wid1)
hold on
plot(norm_mp_dn_21z,'k-', 'LineWidth', line_wid1)
plot(nonlin_error_mp_dn,'k--', 'LineWidth', line_wid)
plot(linear_error_dp_dn,'LineWidth', line_wid,'Color',[0.5 0.5 0.5])
plot(linear_error_mp_dn,'k', 'LineWidth', line_wid)
% plot(linear_error_mp_dn,'k', 'LineWidth', line_wid)
hline = refline([0 0]);
set(hline,'Color','k', 'LineWidth', line_wid)
title('24 Hour Forecast Error','FontSize',fontsize,'FontName','TimesNewRoman')
xlim([1 length(dates_plot)])
set(gca,'XTick',[1 7 13 19 25 32])
xticks = str2num(get(gca,'XTickLabel'));
set(gca,'XTickLabel',{[num2str(dates_plot(xticks)) months_extra]})
ylabel('Energy (Jkg^{-1})','FontSize',fontsize,'FontName','TimesNewRoman')
xlabel('Date','FontSize',fontsize,'FontName','TimesNewRoman')
ylim([-4 12])
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
legend('e(x_b)','e(x_a)','Nonlinear Error','Dry Lin Error','Moist Lin Error','Location','East')

saveas(gcf,'odsstats_lin_vs_nonlin_dn.eps', 'psc2')

figure('Position',[725   398   439   290])

plot(norm_mp_mn_15z,'k--','LineWidth', line_wid1)
hold on
plot(norm_mp_mn_21z,'k-', 'LineWidth', line_wid1)
plot(nonlin_error_mp_mn,'k--', 'LineWidth', line_wid)
plot(linear_error_dp_mn,'LineWidth', line_wid,'Color',[0.5 0.5 0.5])
plot(linear_error_mp_mn,'k', 'LineWidth', line_wid)
% plot(linear_error_mp_dn,'k', 'LineWidth', line_wid)
hline = refline([0 0]);
set(hline,'Color','k', 'LineWidth', line_wid)
title('24 Hour Forecast Error','FontSize',fontsize,'FontName','TimesNewRoman')
xlim([1 length(dates_plot)])
set(gca,'XTick',[1 7 13 19 25 32])
xticks = str2num(get(gca,'XTickLabel'));
set(gca,'XTickLabel',{[num2str(dates_plot(xticks)) months_extra]})
ylabel('Energy (Jkg^{-1})','FontSize',fontsize,'FontName','TimesNewRoman')
xlabel('Date','FontSize',fontsize,'FontName','TimesNewRoman')
ylim([-4 12])
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
legend('e(x_b)','e(x_a)','Nonlinear Error','Dry Lin Error','Moist Lin Error','Location','East')

saveas(gcf,'odsstats_lin_vs_nonlin_mn.eps', 'psc2')


ttest_dry = ttest2(linear_error_dp_dn, linear_error_mp_dn)
ttest_wet = ttest2(linear_error_dp_mn, linear_error_mp_mn)
