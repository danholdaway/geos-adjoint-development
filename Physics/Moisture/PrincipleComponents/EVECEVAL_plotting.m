close all
clear
clc

T_evecs1 = ncread('EVECSEVALS.nc4','T_EVECTS');
Q_evecs1 = ncread('EVECSEVALS.nc4','Q_EVECTS');
U_evecs1 = ncread('EVECSEVALS.nc4','U_EVECTS');
V_evecs1 = ncread('EVECSEVALS.nc4','V_EVECTS');

T_evals1 = ncread('EVECSEVALS.nc4','T_EVALUS');
Q_evals1 = ncread('EVECSEVALS.nc4','Q_EVALUS');
U_evals1 = ncread('EVECSEVALS.nc4','U_EVALUS');
V_evals1 = ncread('EVECSEVALS.nc4','V_EVALUS');

dsig = [0.1339988 0.1453457 0.1576430 0.1709903 0.1854701 0.2011475 0.2181856 0.2333368 0.2032636 0.2032642 ...
  0.2032637 0.2032639 0.2032640 0.2032639 0.2032639 0.2032639 0.2032639 0.1659643 0.1659643 0.1659643    ...
  0.1659634 0.1659651 0.1484429 0.1285559 0.1285543 0.1285559 0.1285548 0.1285554 0.1285548 0.1285559    ...
  0.1285559 0.1285538 0.1285554 0.1285569 0.1285548];

ktop = 72 - length(T_evecs1) + 1;

Nk = length(T_evecs1);

Nk_vec = [1:length(T_evecs1)] + ktop - 1;

%Check Normalisation
sum((T_evecs1(:,1).*dsig').^2)

asd

for mode = 1:Nk;

    plot(sqrt(T_evals1),'o')
    hold on
%     plot(mode,T_evals1(mode),'ro')
    hold off
    
    figure
    subplot(1,4,1)
    plot(T_evecs1(:,mode),Nk_vec)
    set(gca,'YDir','reverse')
    ylim([Nk_vec(1) Nk_vec(end)])
    
    subplot(1,4,2)
    plot(Q_evecs1(:,mode),Nk_vec)
    set(gca,'YDir','reverse')
    ylim([Nk_vec(1) Nk_vec(end)])
    
    subplot(1,4,3)
    plot(U_evecs1(:,mode),Nk_vec)
    set(gca,'YDir','reverse')
    ylim([Nk_vec(1) Nk_vec(end)])
    
    subplot(1,4,4)
    plot(V_evecs1(:,mode),Nk_vec)
    set(gca,'YDir','reverse')
    ylim([Nk_vec(1) Nk_vec(end)])

    pause
    close

end

x = -10:0.1:10;
sd = 1;
mu = 0;
y = 1/(sd*sqrt(2*pi))*exp(-(x-mu).^2./(2*sd^2));

% figure
% plot(x,y)
