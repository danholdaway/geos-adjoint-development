close all
clear
clc

%Dates are for the starting point (take one away)
startdate = '20121226';
enddate   = '20130130';
% enddate   = '20130123';

expid1  = 'x0011dh_a';
expid2  = 'x0011dh_b';
expid3  = 'x0011dh_c';
expid4  = 'x0011dh_d';
expid5  = 'x0011dh_e';

daterange = (datenum(startdate, 'yyyymmdd') : 1 : datenum(enddate, 'yyyymmdd') )';

numdays = length(daterange);

%Add two days to get end point of 24h forecast.
dates = [datestr(daterange, 'yyyymmdd') ; datestr(daterange(end)+1, 'yyyymmdd') ; datestr(daterange(end)+2, 'yyyymmdd')];

%Dry Norm
norm_e00_15z = zeros(numdays,1);
norm_e00_21z = zeros(numdays,1);

%Moist Norm (epsilon = 0.3)
norm_e03_15z = zeros(numdays,1);
norm_e03_21z = zeros(numdays,1);

%Moist Norm (epsilon = 1.0)
norm_e10_15z = zeros(numdays,1);
norm_e10_21z = zeros(numdays,1);

%Norm for PPM9.
norm_def_ppm9_15z = zeros(numdays,1);
norm_def_ppm9_21z = zeros(numdays,1);

for i = 1:numdays
       
    %READ IN DRY AND MOIST (epsilon = 0.3) NORMS - x0011a/x0011b
    path = ['/archive/u/drholdaw/',expid2,'/prog/Y',dates(i,1:4),'/M',dates(i,5:6),'/D',dates(i,7:8),'/H15'];
    cd (path)
    filename = [expid2,'.Jnormf.',dates(i,1:4),dates(i,5:6),dates(i,7:8),'_15z+',dates(i+2,1:4),dates(i+2,5:6),dates(i+2,7:8),'_00z.txt'];
    fid = fopen(filename,'r');
    D_15 = textscan(fid,'%s %f %s %f %s %s %s \n %s %s %s %s %s %s \n %f %f %f %f %f %f \n %s %f %s %f %s %s %s \n %s %s %s %s %s %s \n %f %f %f %f %f %f');
    fclose(fid);

    norm_e00_15z(i) = D_15{1,14};
    norm_e03_15z(i) = D_15{1,33};
    
    path = ['/archive/u/drholdaw/',expid2,'/prog/Y',dates(i,1:4),'/M',dates(i,5:6),'/D',dates(i,7:8),'/H21'];
    cd (path)
    filename = [expid2,'.Jnormf.',dates(i,1:4),dates(i,5:6),dates(i,7:8),'_21z+',dates(i+2,1:4),dates(i+2,5:6),dates(i+2,7:8),'_00z.txt'];
    fid = fopen(filename,'r');
    D_21 = textscan(fid,'%s %f %s %f %s %s %s \n %s %s %s %s %s %s \n %f %f %f %f %f %f \n %s %f %s %f %s %s %s \n %s %s %s %s %s %s \n %f %f %f %f %f %f');
    fclose(fid);

    norm_e00_21z(i) = D_21{1,14};
    norm_e03_21z(i) = D_21{1,33};

    clear D_21 D_15
                          
    %READ IN MOIST (epsilon = 1.0) NORM - x0011c/x0011d - fcst and NL norms are the same.
    path = ['/archive/u/drholdaw/',expid3,'/prog/Y',dates(i,1:4),'/M',dates(i,5:6),'/D',dates(i,7:8),'/H15'];
    cd (path)
    filename = [expid3,'.Jnormf.',dates(i,1:4),dates(i,5:6),dates(i,7:8),'_15z+',dates(i+2,1:4),dates(i+2,5:6),dates(i+2,7:8),'_00z.txt'];
    fid = fopen(filename,'r');
    D_15 = textscan(fid,'%s %f %s %f %s %s %s \n %s %s %s %s %s %s \n %f %f %f %f %f %f \n %s %f %s %f %s %s %s \n %s %s %s %s %s %s \n %f %f %f %f %f %f');
    fclose(fid);

    norm_e10_15z(i) = D_15{1,14};
    
    path = ['/archive/u/drholdaw/',expid3,'/prog/Y',dates(i,1:4),'/M',dates(i,5:6),'/D',dates(i,7:8),'/H21'];
    cd (path)
    filename = [expid3,'.Jnormf.',dates(i,1:4),dates(i,5:6),dates(i,7:8),'_21z+',dates(i+2,1:4),dates(i+2,5:6),dates(i+2,7:8),'_00z.txt'];
    fid = fopen(filename,'r');
    D_21 = textscan(fid,'%s %f %s %f %s %s %s \n %s %s %s %s %s %s \n %f %f %f %f %f %f \n %s %f %s %f %s %s %s \n %s %s %s %s %s %s \n %f %f %f %f %f %f');
    fclose(fid);

    norm_e10_21z(i) = D_21{1,14};

    clear D_21 D_15
                  
    cd /home/drholdaw/Lin_Moist_Physics/Moist_Dev_MWR/
                              
    %READ IN PPM9 - x0011e
    path = ['/archive/u/drholdaw/',expid5,'/prog/Y',dates(i,1:4),'/M',dates(i,5:6),'/D',dates(i,7:8),'/H15'];
    cd (path)
    filename = [expid5,'.Jnormf.',dates(i,1:4),dates(i,5:6),dates(i,7:8),'_15z+',dates(i+2,1:4),dates(i+2,5:6),dates(i+2,7:8),'_00z.txt'];
    fid = fopen(filename,'r');
    D_15 = textscan(fid,'%s %f %s %f %s %s %s \n %s %s %s %s %s %s \n %f %f %f %f %f %f \n %s %f %s %f %s %s %s \n %s %s %s %s %s %s \n %f %f %f %f %f %f');
    fclose(fid);

    norm_def_ppm9_15z(i) = D_15{1,14};
    
    path = ['/archive/u/drholdaw/',expid5,'/prog/Y',dates(i,1:4),'/M',dates(i,5:6),'/D',dates(i,7:8),'/H21'];
    cd (path)
    filename = [expid5,'.Jnormf.',dates(i,1:4),dates(i,5:6),dates(i,7:8),'_21z+',dates(i+2,1:4),dates(i+2,5:6),dates(i+2,7:8),'_00z.txt'];
    fid = fopen(filename,'r');
    D_21 = textscan(fid,'%s %f %s %f %s %s %s \n %s %s %s %s %s %s \n %f %f %f %f %f %f \n %s %f %s %f %s %s %s \n %s %s %s %s %s %s \n %f %f %f %f %f %f');
    fclose(fid);

    norm_def_ppm9_21z(i) = D_21{1,14};

    clear D_21 D_15
                  
    cd /home/drholdaw/Lin_Moist_Physics/Moist_Dev_MWR/
    
end

nonlin_error_e00 = norm_e00_21z - norm_e00_15z;
nonlin_error_e03 = norm_e03_21z - norm_e03_15z;
nonlin_error_e10 = norm_e10_21z - norm_e10_15z;

nonlin_error_def_ppm9 = norm_def_ppm9_21z - norm_def_ppm9_15z;

%NOW READ IN LINEAR MODEL ERORRS
cd /discover/nobackup/drholdaw/odsstats/x0011/

% Dry physics with dry norm
fid = fopen('x0011_e00_pd.txt');
C = textscan(fid,'%f %f %s %f %f');
fclose(fid);
linear_error_dp_e00 = C{5}; clear C

% Dry physics with moist (ep = 0.3) norm
fid = fopen('x0011_e03_pd.txt');
C = textscan(fid,'%f %f %s %f %f');
fclose(fid);
linear_error_dp_e03 = C{5}; clear C

% Dry physics with moist (ep = 1.0) norm
fid = fopen('x0011_e10_pd.txt');
C = textscan(fid,'%f %f %s %f %f');
fclose(fid);
linear_error_dp_e10 = C{5}; clear C

% Moist physics with dry norm
fid = fopen('x0011_e00_pm.txt');
C = textscan(fid,'%f %f %s %f %f');
fclose(fid);
linear_error_mp_e00 = C{5}; clear C

% Moist physics with moist (ep = 0.3) norm
fid = fopen('x0011_e03_pm.txt');
C = textscan(fid,'%f %f %s %f %f');
fclose(fid);
linear_error_mp_e03 = C{5}; clear C

% Moist physics with moist (ep = 1.0) norm
fid = fopen('x0011_e10_pm.txt');
C = textscan(fid,'%f %f %s %f %f');
fclose(fid);
linear_error_mp_e10 = C{5}; clear C

% Default but with PPM9
fid = fopen('x0011_def_ppm9.txt');
C = textscan(fid,'%f %f %s %f %f');
fclose(fid);
linear_error_def_ppm9 = C{5}; clear C

cd /home/drholdaw/Lin_Moist_Physics/Moist_Dev_MWR/

fprintf('Percentage Dry   Pysics Dry (0.0) Norm = %g \n', 100*mean(linear_error_dp_e00./nonlin_error_e00))
fprintf('Percentage Moist Pysics Dry (0.0) Norm = %g \n \n', 100*mean(linear_error_mp_e00./nonlin_error_e00))

fprintf('Percentage Dry   Pysics Moist (0.3) Norm = %g \n', 100*mean(linear_error_dp_e03./nonlin_error_e03))
fprintf('Percentage Moist Pysics Moist (0.3) Norm = %g \n \n', 100*mean(linear_error_mp_e03./nonlin_error_e03))

fprintf('Percentage Dry   Pysics Moist (1.0) Norm = %g \n', 100*mean(linear_error_dp_e10./nonlin_error_e10))
fprintf('Percentage Moist Pysics Moist (1.0) Norm = %g \n \n', 100*mean(linear_error_mp_e10./nonlin_error_e10))

fprintf('Percentage With PPM9 = %g \n \n', 100*mean(linear_error_def_ppm9./nonlin_error_e03))

dates_plot = str2num(dates(2:end-1,end-1:end));


line_wid = 1.0;
line_wid1 = 1.5;
fontsize = 10;

months_extra = ['01-Jan-13'; '07-Jan-13'; '13-Jan-13'; '19-Jan-13'; '25-Jan-13'; '31-Jan-13'];

figure('Position',[125   398   439   290])

plot(norm_e00_15z,'k--','LineWidth', line_wid1)
hold on
plot(norm_e00_21z,'k-', 'LineWidth', line_wid1)
plot(nonlin_error_e00,'k', 'LineWidth', line_wid)

plot(linear_error_dp_e00,'LineWidth', line_wid,'Color',[0.5 0.5 0.5])
plot(linear_error_mp_e00,'r', 'LineWidth', line_wid)

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

plot(norm_e03_15z,'k--','LineWidth', line_wid1)
hold on
plot(norm_e03_21z,'k-', 'LineWidth', line_wid1)
plot(nonlin_error_e03,'k', 'LineWidth', line_wid)

plot(linear_error_dp_e03,'LineWidth', line_wid,'Color',[0.5 0.5 0.5])
plot(linear_error_mp_e03,'r', 'LineWidth', line_wid)

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

plot(norm_e10_15z,'k--','LineWidth', line_wid1)
hold on
plot(norm_e10_21z,'k-', 'LineWidth', line_wid1)
plot(nonlin_error_e10,'k', 'LineWidth', line_wid)

plot(linear_error_dp_e10,'LineWidth', line_wid,'Color',[0.5 0.5 0.5])
plot(linear_error_mp_e10,'r', 'LineWidth', line_wid)

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

plot(norm_def_ppm9_15z,'k--','LineWidth', line_wid1)
hold on
plot(norm_def_ppm9_21z,'k-', 'LineWidth', line_wid1)
plot(nonlin_error_def_ppm9,'k', 'LineWidth', line_wid)

plot(linear_error_mp_e03,'LineWidth', line_wid,'Color',[0.5 0.5 0.5])
plot(linear_error_def_ppm9,'r', 'LineWidth', line_wid)

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
