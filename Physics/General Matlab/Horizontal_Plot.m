%Pre-am
close all
clear
clc
mydir = pwd;

%Matlab coast lines
load coast
coast_lat = lat; clear lat
coast_lon = long; clear long

%Plotting things
fontsize = 16;
line_wid = 1.1;

%Switch to color zeros as white
whiteout = 0;

%Lon/lat/level you want to plot
lon_min = -105;
lon_max = -70;
lat_min = 35;
lat_max = 50;
lev_plo = 68;

%Adjoint initial conditions box
center = [-87.75 42.0];
eastwest = 1;
nortsout = 1;
NW = [center(1)-eastwest/2 center(2)+nortsout/2];
NE = [center(1)+eastwest/2 center(2)+nortsout/2];
SE = [center(1)+eastwest/2 center(2)-nortsout/2];
SW = [center(1)-eastwest/2 center(2)-nortsout/2];
lineN = NW(1):0.1:NE(1);
lineE = SE(2):0.1:NE(2);
lineS = SW(1):0.1:SE(1);
lineW = SW(2):0.1:NW(2);

InitAdj = [359-abs(NW(1))+1 359-abs(NE(1))+1 SE(2) NE(2)]

%File to be read
filename = 'GreatLakes_Hurricane_C180_e34.prog.eta.20101025_00z.nc4';

%Variable to be read
% var = 'u';
% var = 'v';
% var = 'tv';
var = 'sphu';
% var = 'delp';

%Location of file
dir = '/gpfsm/dnb32/bthoover/wrk.greatlakes';

%Move to directory and read grid
cd(dir)
lon = ncread(filename,'lon'); im = length(lon);
lat = ncread(filename,'lat'); jm = length(lat);
lev = ncread(filename,'lev'); lm = length(lev);

%Compute nearest index for chosen plotting region
tmp = abs(lon-lon_min);
lon_min = find(tmp == min(tmp(:)));
tmp = abs(lon-lon_max);
lon_max = find(tmp == min(tmp(:)));
tmp = abs(lat-lat_min);
lat_min = find(tmp == min(tmp(:)));
tmp = abs(lat-lat_max);
lat_max = find(tmp == min(tmp(:)));

%Read fields for the region and level
f1 = ncread(filename,var,[lon_min lat_min lev_plo 1],[lon_max-lon_min+1 lat_max-lat_min+1 1 1]);

%Get maxmin for the field
maxf1 = max(f1(:));
minf1 = min(f1(:));

%Square up coloring if negative values
if minf1 < 0
    maxabs = max(abs([maxf1 maxf1]));
    maxf1 = maxabs;
    minf1 = -maxabs;
    grey = 0.1;
    
    % Option to make 0 of colormap white
    cmap = colormap;
    if whiteout == 1
        cmap(32,:) = [1 1 1];
        cmap(33,:) = [1 1 1];
    end
    close
else
    grey = 0.25;
    % Option to make 0 of colormap white
    cmap = colormap;
    if whiteout == 1
        cmap(1,:) = [1 1 1];
    end
    close
end

%Number of contour intervals
numints = 20;
step = (maxf1-minf1)/numints;

%Plot the field
figure
set(gcf,'position',[152 434 1021 420])
contourf(lon(lon_min:lon_max),lat(lat_min:lat_max),f1','LineStyle','none','LevelStep',step)
hold on
plot(coast_lon,coast_lat,'Color',[grey grey grey],'LineWidth',line_wid)
plot(lineN,NW(2)*ones(1,length(lineN)),'k','LineWidth',line_wid)
plot(lineS,SW(2)*ones(1,length(lineS)),'k','LineWidth',line_wid)
plot(NE(1)*ones(1,length(lineE)),lineE,'k','LineWidth',line_wid)
plot(NW(1)*ones(1,length(lineW)),lineW,'k','LineWidth',line_wid)
caxis([minf1 maxf1])
box on
colormap(cmap)
colorbar
title([var,' at level ',num2str(lev_plo)],'FontSize',fontsize,'FontName','TimesNewRoman')
ylabel('Latitude (\circ)','FontSize',fontsize,'FontName','TimesNewRoman')
xlabel('Longitude (\circ)','FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

cd(mydir)