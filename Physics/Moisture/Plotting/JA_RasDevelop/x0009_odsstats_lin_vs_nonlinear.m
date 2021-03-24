close all
clear
clc

%Dates are for the starting point (take one away)
startdate = '20120630';
enddate   = '20120730';

expid1 = 'x0009a';
expid2 = 'x0009b';
expid3 = 'x0009c';
expid4 = 'x0009d';

daterange = (datenum(startdate, 'yyyymmdd') : 1 : datenum(enddate, 'yyyymmdd') )';

numdays = length(daterange);

%Add two days to get end point of 24h forecast.
dates = [datestr(daterange, 'yyyymmdd') ; datestr(daterange(end)+1, 'yyyymmdd') ; datestr(daterange(end)+2, 'yyyymmdd')];

%Dry Norm
norm_dn00_15z = zeros(numdays,1);
norm_dn00_21z = zeros(numdays,1);

%Moist Norm (epsilon = 0.3)
norm_mn03_15z = zeros(numdays,1);
norm_mn03_21z = zeros(numdays,1);

%Moist Norm (epsilon = 1.0)
norm_mn10_15z = zeros(numdays,1);
norm_mn10_21z = zeros(numdays,1);


for i = 1:numdays
       
    %READ IN DRY AND MOIST (epsilon = 0.3) NORMS - x0009a
    path = ['/archive/u/drholdaw/',expid1,'/prog/Y',dates(i,1:4),'/M',dates(i,5:6),'/D',dates(i,7:8),'/H15'];
    cd (path)
    filename = [expid1,'.Jnormf.',dates(i,1:4),dates(i,5:6),dates(i,7:8),'_15z+',dates(i+2,1:4),dates(i+2,5:6),dates(i+2,7:8),'_00z.txt'];
    fid = fopen(filename,'r');
    D_15 = textscan(fid,'%s %f %s %f %s %s %s \n %s %s %s %s %s %s \n %f %f %f %f %f %f \n %s %f %s %f %s %s %s \n %s %s %s %s %s %s \n %f %f %f %f %f %f');
    fclose(fid);

    norm_dn00_15z(i) = D_15{1,14};
    norm_mn03_15z(i) = D_15{1,33};
    
    path = ['/archive/u/drholdaw/',expid1,'/prog/Y',dates(i,1:4),'/M',dates(i,5:6),'/D',dates(i,7:8),'/H21'];
    cd (path)
    filename = [expid1,'.Jnormf.',dates(i,1:4),dates(i,5:6),dates(i,7:8),'_21z+',dates(i+2,1:4),dates(i+2,5:6),dates(i+2,7:8),'_00z.txt'];
    fid = fopen(filename,'r');
    D_21 = textscan(fid,'%s %f %s %f %s %s %s \n %s %s %s %s %s %s \n %f %f %f %f %f %f \n %s %f %s %f %s %s %s \n %s %s %s %s %s %s \n %f %f %f %f %f %f');
    fclose(fid);

    norm_dn00_21z(i) = D_21{1,14};
    norm_mn03_21z(i) = D_21{1,33};

    clear D_21 D_15
                          
    %READ IN MOIST (epsilon = 1.0) NORM - x0009c
    path = ['/archive/u/drholdaw/',expid3,'/prog/Y',dates(i,1:4),'/M',dates(i,5:6),'/D',dates(i,7:8),'/H15'];
    cd (path)
    filename = [expid3,'.Jnormf.',dates(i,1:4),dates(i,5:6),dates(i,7:8),'_15z+',dates(i+2,1:4),dates(i+2,5:6),dates(i+2,7:8),'_00z.txt'];
    fid = fopen(filename,'r');
    D_15 = textscan(fid,'%s %f %s %f %s %s %s \n %s %s %s %s %s %s \n %f %f %f %f %f %f \n %s %f %s %f %s %s %s \n %s %s %s %s %s %s \n %f %f %f %f %f %f');
    fclose(fid);

    norm_mn10_15z(i) = D_15{1,14};
    
    path = ['/archive/u/drholdaw/',expid3,'/prog/Y',dates(i,1:4),'/M',dates(i,5:6),'/D',dates(i,7:8),'/H21'];
    cd (path)
    filename = [expid3,'.Jnormf.',dates(i,1:4),dates(i,5:6),dates(i,7:8),'_21z+',dates(i+2,1:4),dates(i+2,5:6),dates(i+2,7:8),'_00z.txt'];
    fid = fopen(filename,'r');
    D_21 = textscan(fid,'%s %f %s %f %s %s %s \n %s %s %s %s %s %s \n %f %f %f %f %f %f \n %s %f %s %f %s %s %s \n %s %s %s %s %s %s \n %f %f %f %f %f %f');
    fclose(fid);

    norm_mn10_21z(i) = D_21{1,14};

    clear D_21 D_15
                  
    cd /home/drholdaw/Lin_Moist_Physics/Moist_Dev_MWR/
    
end

nonlin_error_dn00 = norm_dn00_21z - norm_dn00_15z;
nonlin_error_mn03 = norm_mn03_21z - norm_mn03_15z;
nonlin_error_mn10 = norm_mn10_21z - norm_mn10_15z;


%NOW READ IN LINEAR MODEL ERORRS
cd /discover/nobackup/drholdaw/odsstats/x0009/

% Dry physics with dry norm
fid = fopen('x0009_nd00_pd.txt');
C = textscan(fid,'%f %f %s %f %f');
fclose(fid);
linear_error_dp_dn00 = C{5}; clear C

% Dry physics with moist (ep = 0.3) norm
fid = fopen('x0009_nm03_pd.txt');
C = textscan(fid,'%f %f %s %f %f');
fclose(fid);
linear_error_dp_mn03 = C{5}; clear C

% Dry physics with moist (ep = 1.0) norm
fid = fopen('x0009_nm10_pd.txt');
C = textscan(fid,'%f %f %s %f %f');
fclose(fid);
linear_error_dp_mn10 = C{5}; clear C

% Moist physics with dry norm
fid = fopen('x0009_nd00_pm.txt');
C = textscan(fid,'%f %f %s %f %f');
fclose(fid);
linear_error_mp_dn00 = C{5}; clear C

% % Moist physics with moist (ep = 0.3) norm
fid = fopen('x0009_nm03_pm.txt');
C = textscan(fid,'%f %f %s %f %f');
fclose(fid);
linear_error_mp_mn03 = C{5}; clear C

% Moist physics with moist (ep = 1.0) norm
fid = fopen('x0009_nm10_pm.txt');
C = textscan(fid,'%f %f %s %f %f');
fclose(fid);
linear_error_mp_mn10 = C{5}; clear C

cd /home/drholdaw/Lin_Moist_Physics/Moist_Dev_MWR/

fprintf('Percentage Dry   Pysics Dry (0.0) Norm = %g \n', 100*mean(linear_error_dp_dn00./nonlin_error_dn00))
fprintf('Percentage Moist Pysics Dry (0.0) Norm = %g \n \n', 100*mean(linear_error_mp_dn00./nonlin_error_dn00))
% 
fprintf('Percentage Dry   Pysics Moist (0.3) Norm = %g \n', 100*mean(linear_error_dp_mn03./nonlin_error_mn03))
fprintf('Percentage Moist Pysics Moist (0.3) Norm = %g \n \n', 100*mean(linear_error_mp_mn03./nonlin_error_mn03))
% 
fprintf('Percentage Dry   Pysics Moist (1.0) Norm = %g \n', 100*mean(linear_error_dp_mn10./nonlin_error_mn10))
fprintf('Percentage Moist Pysics Moist (1.0) Norm = %g \n \n', 100*mean(linear_error_mp_mn10./nonlin_error_mn10))


dates_plot = str2num(dates(2:end-1,end-1:end));


line_wid = 1.0;
line_wid1 = 1.5;
fontsize = 10;

months_extra = ['01-Jul-12'; '07-Jul-12'; '13-Jul-12'; '19-Jul-12'; '25-Jul-12'; '31-Jul-12'];

figure('Position',[125   398   439   290])

plot(norm_dn00_15z,'k--','LineWidth', line_wid1)
hold on
plot(norm_dn00_21z,'k-', 'LineWidth', line_wid1)
plot(nonlin_error_dn00,'k--', 'LineWidth', line_wid)

plot(linear_error_dp_dn00,'LineWidth', line_wid,'Color',[0.5 0.5 0.5])
plot(linear_error_mp_dn00,'k', 'LineWidth', line_wid)

hline = refline([0 0]);
set(hline,'Color','k', 'LineWidth', line_wid)

title('24 Hour Forecast Error','FontSize',fontsize,'FontName','TimesNewRoman')
ylabel('Energy (Jkg^{-1})','FontSize',fontsize,'FontName','TimesNewRoman')
xlabel('Date','FontSize',fontsize,'FontName','TimesNewRoman')
legend('15z','21z','Nonlinear Error','Dry Lin Error','Moist Lin Error','Location','East')

xlim([1 31])
set(gca,'XTick',[1 7 13 19 25 31])
xticks = str2double(get(gca,'XTickLabel'));
set(gca,'XTickLabel',{months_extra})


figure('Position',[125   398   439   290])

plot(norm_mn03_15z,'k--','LineWidth', line_wid1)
hold on
plot(norm_mn03_21z,'k-', 'LineWidth', line_wid1)
plot(nonlin_error_mn03,'k--', 'LineWidth', line_wid)

plot(linear_error_dp_mn03,'LineWidth', line_wid,'Color',[0.5 0.5 0.5])
plot(linear_error_mp_mn03,'k', 'LineWidth', line_wid)

hline = refline([0 0]);
set(hline,'Color','k', 'LineWidth', line_wid)

title('24 Hour Forecast Error','FontSize',fontsize,'FontName','TimesNewRoman')
ylabel('Energy (Jkg^{-1})','FontSize',fontsize,'FontName','TimesNewRoman')
xlabel('Date','FontSize',fontsize,'FontName','TimesNewRoman')
legend('15z','21z','Nonlinear Error','Dry Lin Error','Moist Lin Error','Location','East')

xlim([1 31])
set(gca,'XTick',[1 7 13 19 25 31])
xticks = str2double(get(gca,'XTickLabel'));
set(gca,'XTickLabel',{months_extra})





figure('Position',[125   398   439   290])

plot(norm_mn10_15z,'k--','LineWidth', line_wid1)
hold on
plot(norm_mn10_21z,'k-', 'LineWidth', line_wid1)
plot(nonlin_error_mn10,'k--', 'LineWidth', line_wid)

plot(linear_error_dp_mn10,'LineWidth', line_wid,'Color',[0.5 0.5 0.5])
plot(linear_error_mp_mn10,'k', 'LineWidth', line_wid)

hline = refline([0 0]);
set(hline,'Color','k', 'LineWidth', line_wid)

title('24 Hour Forecast Error','FontSize',fontsize,'FontName','TimesNewRoman')
ylabel('Energy (Jkg^{-1})','FontSize',fontsize,'FontName','TimesNewRoman')
xlabel('Date','FontSize',fontsize,'FontName','TimesNewRoman')
legend('15z','21z','Nonlinear Error','Dry Lin Error','Moist Lin Error','Location','East')

xlim([1 31])
set(gca,'XTick',[1 7 13 19 25 31])
xticks = str2double(get(gca,'XTickLabel'));
set(gca,'XTickLabel',{months_extra})






ttest_dry = ttest2(linear_error_dp_dn, linear_error_mp_dn)
ttest_wet = ttest2(linear_error_dp_mn, linear_error_mp_mn)
