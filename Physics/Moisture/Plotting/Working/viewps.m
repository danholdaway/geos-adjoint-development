close all
clear
clc

startdate = '20120316';
enddate   = '20120331';

daterange = (datenum(startdate, 'yyyymmdd') : 1 : datenum(enddate, 'yyyymmdd') )';

numdays = length(daterange);

%Add two days to get end point of 24h forecast.
dates = [datestr(daterange, 'yyyymmdd') ; datestr(daterange(end)+1, 'yyyymmdd') ; datestr(daterange(end)+2, 'yyyymmdd')];

boxhor = -37:-34;
boxver = 45:48;

line_wid = 1;
fontsize = 14;

for i = 5%:numdays

    disp([dates(i,1:4) dates(i,5:6) dates(i,7:8)])
    
    path = ['/archive/u/drholdaw/iau590/prog/Y',dates(i,1:4),'/M',dates(i,5:6),'/D',dates(i,7:8),'/H21'];
    cd (path)


    filename = ['iau590.prog.eta.',dates(i,1:4),dates(i,5:6),dates(i,7:8),'_21z+',dates(i+2,1:4),dates(i+2,5:6),dates(i+2,7:8),'_00z.nc4'];
    Ps = ncread(filename,'ps');
    lon = ncread(filename,'lon');
    lat = ncread(filename,'lat');
    
    cd /discover/nobackup/drholdaw/

    % contourf(Ps')

    [C,h] = contour(Ps'./100,960:4:1020);shading flat
    % set(h,'LevelStep',get(h,'LevelStep')/10)
    
    hold on
    plot(boxhor,boxver(1)*ones(1,4),'k')
    plot(boxhor,boxver(end)*ones(1,4),'k')
    plot(boxhor(1)*ones(1,4),boxver,'k')
    plot(boxhor(end)*ones(1,4),boxver,'k')
    
%     hold on; load topo.mat;
%     topo_grenwich(1:180,181:360) = topo(1:180,1:180);
%     topo_grenwich(1:180,1:180) = topo(1:180,181:360);
%     contour(-180:179,-89:90,topo_grenwich,[0 0],'k','LineWidth',1.0)
%     axis equal; box on; xlim([-180 180]); ylim([-90 90])
%     set(gca,'XTick',[-180 -120 -60 0 60 120 180],'Ytick',[-90 -60 -30 0 30 60 90]);
    set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
    
%     pause

end