close all
clear
clc

cd /discover/nobackup/drholdaw/ExperimentData/jacobianVSheatingrate/

jacobian_container
cd /home/drholdaw/Lin_Moist_Physics/DAS_plotting/

% J = J4;

fontsize = 14;

k = length(J)/2;

figure
subplot(2,2,1)
contour(J(1:k,1:k))
colorbar
set(gca,'YDir','reverse')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')


subplot(2,2,2)
contour(J(1:k,k+1:2*k))
colorbar
set(gca,'YDir','reverse')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')


subplot(2,2,3)
contour(J(k+1:2*k,1:k))
colorbar
set(gca,'YDir','reverse')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')


subplot(2,2,4)
contour(J(k+1:2*k,k+1:2*k))
colorbar
set(gca,'YDir','reverse')

set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')


