close all
clear
clc

cd /discover/nobackup/drholdaw/ExperimentData/jacobianVSheatingrate
% asd
jacobian_container2
cd /home/drholdaw/Lin_Moist_Physics/DAS_plotting/

fontsize = 10;

k = length(J)/2;

% col_plot_t = 68;
% col_plot_q = 68;
% J(:,1:col_plot_t-1) = 0.0;
% J(:,col_plot_t+1:k) = 0.0;
% J(:,k+1:k+col_plot_q-1) = 0.0;
% J(:,k+col_plot_q+1:k+k) = 0.0;

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

fprintf('Max J1 = %g \n', max(max(J(1:k,1:k))))
fprintf('Max J2 = %g \n', max(max(J(1:k,k+1:2*k))))
fprintf('Max J3 = %g \n', max(max(J(k+1:2*k,1:k))))
fprintf('Max J4 = %g \n\n', max(max(J(k+1:2*k,k+1:2*k))))

figure
subplot(1,2,1)
plot(H,1:72)
set(gca,'YDir','reverse')
ylim([30 72])

subplot(1,2,2)
plot(M,1:72)
set(gca,'YDir','reverse')
ylim([30 72])


[X,L] = eig(eye(2*k) + 1200*J);

l = diag(L);

rl = real(l);
il = imag(l);

% figure
% scatter(rl,il)

disp(max(sort(rl)))
