close all
clear
clc

%Dates are for the starting point (take one away)
startdate = '20121231';
enddate   = '20130130';
% enddate   = '20130123';

expid1  = 'x0011dh_a';
expid2  = 'x0011dh_b';

daterange = (datenum(startdate, 'yyyymmdd') : 1 : datenum(enddate, 'yyyymmdd') )';

numdays = length(daterange);

%Add two days to get end point of 24h forecast.
dates = [datestr(daterange, 'yyyymmdd') ; datestr(daterange(end)+1, 'yyyymmdd') ; datestr(daterange(end)+2, 'yyyymmdd')];

%Dry Norm
norm_15z = zeros(numdays,1);
norm_21z = zeros(numdays,1);

for i = 1:numdays
       
    %READ IN DRY AND MOIST (epsilon = 0.3) NORMS - x0011a/x0011b
    path = ['/archive/u/drholdaw/',expid1,'/prog/Y',dates(i,1:4),'/M',dates(i,5:6),'/D',dates(i,7:8),'/H15'];
    cd (path)
    filename = [expid1,'.Jnormf.',dates(i,1:4),dates(i,5:6),dates(i,7:8),'_15z+',dates(i+2,1:4),dates(i+2,5:6),dates(i+2,7:8),'_00z.txt'];
    fid = fopen(filename,'r');
    D_15 = textscan(fid,'%s %f %s %f %s %s %s \n %s %s %s %s %s %s %s %s \n %f %f %f %f %f %f \n %f %f ');
    fclose(fid);

    norm_15z(i) = D_15{1,16};
    
    path = ['/archive/u/drholdaw/',expid1,'/prog/Y',dates(i,1:4),'/M',dates(i,5:6),'/D',dates(i,7:8),'/H21'];
    cd (path)
    filename = [expid1,'.Jnormf.',dates(i,1:4),dates(i,5:6),dates(i,7:8),'_21z+',dates(i+2,1:4),dates(i+2,5:6),dates(i+2,7:8),'_00z.txt'];
    fid = fopen(filename,'r');
    D_21 = textscan(fid,'%s %f %s %f %s %s %s \n %s %s %s %s %s %s %s %s \n %f %f %f %f %f %f \n %f %f ');
    fclose(fid);

    norm_21z(i) = D_21{1,16};

    clear D_21 D_15
                    
    cd /home/drholdaw/Lin_Moist_Physics/DAS_plotting/Cloud' Fraction'/
    
end

nonlin_error = norm_21z - norm_15z


%NOW READ IN LINEAR MODEL ERORRS
cd /discover/nobackup/drholdaw/odsstats/x0011dh_a/
fid = fopen('accum.sum_odsstats.txt');
C = textscan(fid,'%f %f %s %f %f');
fclose(fid);
linear_error_moi = C{5}; clear C

cd /discover/nobackup/drholdaw/odsstats/x0011dh_b/
fid = fopen('accum.sum_odsstats.txt');
C = textscan(fid,'%f %f %s %f %f');
fclose(fid);
linear_error_dry = C{5}; clear C


cd /home/drholdaw/Lin_Moist_Physics/DAS_plotting/Cloud' Fraction'/

fprintf('Percentage Dry   Pysics Dry (0.0) Norm = %g \n',    100*mean(linear_error_dry(1:numdays)./nonlin_error(1:numdays)))
fprintf('Percentage Moist Pysics Dry (0.0) Norm = %g \n \n', 100*mean(linear_error_moi(1:numdays)./nonlin_error(1:numdays)))


dates_plot = str2num(dates(2:end-1,end-1:end));


line_wid = 1.0;
line_wid1 = 1.5;
fontsize = 16;

months_extra = ['01-Jan-13'; '07-Jan-13'; '13-Jan-13'; '19-Jan-13'; '25-Jan-13'; '31-Jan-13'];

figure('Position',[125   398   439   290])

plot(norm_15z,'k--','LineWidth', line_wid1)
hold on
plot(norm_21z,'k-', 'LineWidth', line_wid1)
plot(nonlin_error,'k', 'LineWidth', line_wid)

plot(linear_error_dry(1:numdays),'b', 'LineWidth', line_wid)
plot(linear_error_moi(1:numdays),'r', 'LineWidth', line_wid)

hline = refline([0 0]);
set(hline,'Color','k', 'LineWidth', line_wid)

title('24 Hour Forecast Error','FontSize',fontsize,'FontName','TimesNewRoman')
ylabel('Energy (Jkg^{-1})','FontSize',fontsize,'FontName','TimesNewRoman')
xlabel('Date','FontSize',fontsize,'FontName','TimesNewRoman')
legend('15z','21z','Nonlinear Error','Dry Lin Error','Moist Lin Error','Location','East')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')


xlim([1 31])
set(gca,'XTick',[1 7 13 19 25 31])
xticks = str2double(get(gca,'XTickLabel'));
set(gca,'XTickLabel',{months_extra})






ttest_dry = ttest2(linear_error_dry, linear_error_moi)

