close all
clear
clc

load coast
coast_lat = lat; clear lat
coast_lon = long; clear long

startdate = '20120409';
enddate   = '20120419';
daterange = (datenum(startdate, 'yyyymmdd') : 1 : datenum(enddate, 'yyyymmdd') )';

numdays = length(daterange);

%Add two days to get end point of 24h forecast.
dates = [datestr(daterange, 'yyyymmdd') ; datestr(daterange(end)+1, 'yyyymmdd') ; datestr(daterange(end)+2, 'yyyymmdd')];

%PICK DAY
for i = 1:length(dates-2);
% for i = 2;

    
    fprintf('Current Day = %s-%s-%s \n', dates(i,1:4), dates(i,5:6), dates(i,7:8));
    
    path = ['/archive/u/drholdaw/iau590a_mp/prog_ppm6/Y',dates(i,1:4),'/M',dates(i,5:6),'/D',dates(i,7:8),'/H15'];
    cd (path)

    filename = ['iau590a_mp.fsens_txe.eta.',dates(i,1:4),dates(i,5:6),dates(i,7:8),'_15z+',dates(i+2,1:4),dates(i+2,5:6),dates(i+2,7:8),'_00z-',dates(i+1,1:4),dates(i+1,5:6),dates(i+1,7:8),'_00z.nc4'];

    lat = ncread(filename,'lat');
    lon = ncread(filename,'lon');

    u_ppm6 = ncread(filename,'u');

    path = ['/archive/u/drholdaw/iau590a/prog/Y',dates(i,1:4),'/M',dates(i,5:6),'/D',dates(i,7:8),'/H15'];
    cd (path)

    filename = ['iau590a.fsens_txe.eta.',dates(i,1:4),dates(i,5:6),dates(i,7:8),'_15z+',dates(i+2,1:4),dates(i+2,5:6),dates(i+2,7:8),'_00z-',dates(i+1,1:4),dates(i+1,5:6),dates(i+1,7:8),'_00z.nc4'];

    u_ppm9 = ncread(filename,'u');

    path = ['/archive/u/drholdaw/iau590a_mp/prog/Y',dates(i,1:4),'/M',dates(i,5:6),'/D',dates(i,7:8),'/H15'];
    cd (path)
    
    filename = ['iau590a_mp.fsens_txe.eta.',dates(i,1:4),dates(i,5:6),dates(i,7:8),'_15z+',dates(i+2,1:4),dates(i+2,5:6),dates(i+2,7:8),'_00z-',dates(i+1,1:4),dates(i+1,5:6),dates(i+1,7:8),'_00z.nc4'];
    
    u_ppm1 = ncread(filename,'u');

    cd /home/drholdaw/Lin_Moist_Physics/DAS_plotting/


    k = 50;
%     while k < 100

        max_u = max(max([u_ppm6(:,:,k);u_ppm9(:,:,k)]));
        min_u = min(min([u_ppm6(:,:,k);u_ppm9(:,:,k)]));
%         max_u = max(max([u_ppm6(:,:,k);u_ppm9(:,:,k);u_ppm1(:,:,k)]));
%         min_u = min(min([u_ppm6(:,:,k);u_ppm9(:,:,k);u_ppm1(:,:,k)]));

%         max_u = max(max([u_ppm1(:,:,k)]));
%         min_u = min(min([u_ppm1(:,:,k)]));

        max_abs_u = max([max_u min_u]);

        cint_u = 2*max_u/40;
 
        fprintf('Current level = %g \n', k);

        figure('position',[509    30   770   889])
        subplot(3,1,1)
        [C,h] = contourf(lon,lat,u_ppm1(:,:,k)','LineStyle','none');
%         set(h,'LevelStep',cint_u);
        hold on; 
        plot(coast_lon,coast_lat,'Color','k')
%         caxis([-max_abs_u max_abs_u])
        title('PPM 1')
        colorbar
        subplot(3,1,2)
        [C,h] = contourf(lon,lat,u_ppm6(:,:,k)','LineStyle','none');
%         set(h,'LevelStep',cint_u);
        hold on; 
        plot(coast_lon,coast_lat,'Color','k')
%         caxis([-max_abs_u max_abs_u])
        title('PPM 6')
        colorbar
        subplot(3,1,3)
        [C,h] = contourf(lon,lat,u_ppm9(:,:,k)','LineStyle','none');
%         set(h,'LevelStep',cint_u);
        hold on; 
        plot(coast_lon,coast_lat,'Color','k')
%         caxis([-max_abs_u max_abs_u])
        title('PPM 9')
        colorbar


%         k = input('Next level? \n')
%         close
% 
%     end

    pause
    close

end
