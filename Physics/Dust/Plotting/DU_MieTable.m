close all
clear
clc
mydir = pwd;


cd /discover/nobackup/drholdaw/btmp.25956/sens.20140202.000000/ExtData/g5chem/x/

% file = 'opticsBands_DU.v14_2.nc';
% file = 'opticsBands_SS.v3_3.nc';
% file = 'opticsBands_SU.v1_3.nc';
% file = 'opticsBands_OC.v1_3.nc';
file = 'opticsBands_BC.v1_3.nc';


rh = ncread(file,'rh');
lambda = ncread(file,'lambda');
radius = ncread(file,'radius');

bext = ncread(file,'bext');
bsca = ncread(file,'bsca');

cd(mydir)

figure
set(gcf,'position',[354    30   925   889])
figure
set(gcf,'position',[354    30   925   889])
for i = 1:size(bext,3)
    
    figure(1)
    subplot(3,2,i)
    contourf(rh,lambda,bext(:,:,i))
    title(['Radius - ',num2str(radius(i))])
    colorbar

    figure(2)
    subplot(3,2,i)
    hold on
    for j = 1:length(lambda)
        grey = j/length(lambda);
        plot(rh,bext(j,:,i),'Color',[grey grey grey])
    end
    hold off
end

figure
set(gcf,'position',[354    30   925   889])
figure
set(gcf,'position',[354    30   925   889])
for i = 1:size(bext,3)
    
    figure(3)
    subplot(3,2,i)
    contourf(rh,lambda,bsca(:,:,i))
    title(['Radius - ',num2str(radius(i))])
    colorbar

    figure(4)
    subplot(3,2,i)
    hold on
    for j = 1:length(lambda)
        grey = j/length(lambda);
        plot(rh,bsca(j,:,i),'Color',[grey grey grey])
    end
    hold off
end