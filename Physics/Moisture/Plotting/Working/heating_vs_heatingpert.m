close all
clc

clear
heating_vs_heatingpert_container_10e2

scrsz = get(0,'ScreenSize');

J_thresh_1 = zeros(24,3);
J_thresh_2 = zeros(24,3);
J_thresh_5 = zeros(24,3);
J_thresh_10 = zeros(24,3);

HMp_max = zeros(24,4);

for l = 1:12
    i = HMvsTRAJ(1,1,l);
    HMp_max(l,1) = HMvsTRAJ(i,1,l);
    HMp_max(l,2) = HMvsTRAJ(i-1,1,l);
    HMp_max(l,3) = HMvsTRAJ(i,2,l);
    HMp_max(l,4) = HMvsTRAJ(i-1,2,l);
end

disp(HMp_max)

asd

for l = 1:12;

    k = input('Choose Time Interval: ');
%     k = l;
    
    i = HMvsTRAJ(1,1,k);
    
    HMvsTRAJ_plot = HMvsTRAJ(2:i,:,k);
     
    i = i-1;
    
    i90 = round(90*i/100);
    i95 = round(95*i/100);
    i98 = round(98*i/100);
    i99 = round(99*i/100);
    
    close all
    
    %%%%%%%%%%%%
    figure('visible','on','Position',[0 scrsz(4) 0.5*scrsz(3) scrsz(4)])
    
    [HMvsTRAJ_sort,ind] = sort(HMvsTRAJ_plot(:,3));
    HMvsTRAJ_sort_plot = HMvsTRAJ_plot(ind,:);
    subplot(4,4,1)
    scatter(HMvsTRAJ_sort_plot(:,1), HMvsTRAJ_sort_plot(:,3),'b')
    hold on
    xaxi = [min(HMvsTRAJ_sort_plot(:,1)) max(HMvsTRAJ_sort_plot(:,1))];
    plot(xaxi,HMvsTRAJ_sort_plot(i99,3)*ones(1,2),'y','LineWidth',2.0)
    plot(xaxi,HMvsTRAJ_sort_plot(i98,3)*ones(1,2),'r','LineWidth',2.0)
    plot(xaxi,HMvsTRAJ_sort_plot(i95,3)*ones(1,2),'g','LineWidth',2.0)
    plot(xaxi,HMvsTRAJ_sort_plot(i90,3)*ones(1,2),'k','LineWidth',2.0)
    set(gca,'YScale','log')
    set(gca,'XScale','log')
    
    J_thresh_1(k,1) = HMvsTRAJ_sort_plot(i99,3);
    J_thresh_2(k,1) = HMvsTRAJ_sort_plot(i98,3);
    J_thresh_5(k,1) = HMvsTRAJ_sort_plot(i95,3);
    J_thresh_10(k,1) = HMvsTRAJ_sort_plot(i90,3);
    
    [HMvsTRAJ_sort,ind] = sort(HMvsTRAJ_plot(:,4));
    HMvsTRAJ_sort_plot = HMvsTRAJ_plot(ind,:);
    subplot(4,4,2)
    scatter(HMvsTRAJ_sort_plot(:,1), HMvsTRAJ_sort_plot(:,4),'b')
    hold on
    xaxi = [min(HMvsTRAJ_sort_plot(:,1)) max(HMvsTRAJ_sort_plot(:,1))];
    plot(xaxi,HMvsTRAJ_sort_plot(i99,4)*ones(1,2),'y','LineWidth',2.0)
    plot(xaxi,HMvsTRAJ_sort_plot(i98,4)*ones(1,2),'r','LineWidth',2.0)
    plot(xaxi,HMvsTRAJ_sort_plot(i95,4)*ones(1,2),'g','LineWidth',2.0)
    plot(xaxi,HMvsTRAJ_sort_plot(i90,4)*ones(1,2),'k','LineWidth',2.0)
    set(gca,'YScale','log')
    set(gca,'XScale','log')
    
    J_thresh_1(k,2) = HMvsTRAJ_sort_plot(i99,4);
    J_thresh_2(k,2) = HMvsTRAJ_sort_plot(i98,4);
    J_thresh_5(k,2) = HMvsTRAJ_sort_plot(i95,4);
    J_thresh_10(k,2) = HMvsTRAJ_sort_plot(i90,4);

    [HMvsTRAJ_sort,ind] = sort(HMvsTRAJ_plot(:,5));
    HMvsTRAJ_sort_plot = HMvsTRAJ_plot(ind,:);
    subplot(4,4,3)
    scatter(HMvsTRAJ_sort_plot(:,1), HMvsTRAJ_sort_plot(:,5),'b')
    hold on
    xaxi = [min(HMvsTRAJ_sort_plot(:,1)) max(HMvsTRAJ_sort_plot(:,1))];
    plot(xaxi,HMvsTRAJ_sort_plot(i99,5)*ones(1,2),'y','LineWidth',2.0)
    plot(xaxi,HMvsTRAJ_sort_plot(i98,5)*ones(1,2),'r','LineWidth',2.0)
    plot(xaxi,HMvsTRAJ_sort_plot(i95,5)*ones(1,2),'g','LineWidth',2.0)
    plot(xaxi,HMvsTRAJ_sort_plot(i90,5)*ones(1,2),'k','LineWidth',2.0)
    set(gca,'YScale','log')
    set(gca,'XScale','log')
    
    J_thresh_1(k,3) = HMvsTRAJ_sort_plot(i99,5);
    J_thresh_2(k,3) = HMvsTRAJ_sort_plot(i98,5);
    J_thresh_5(k,3) = HMvsTRAJ_sort_plot(i95,5);
    J_thresh_10(k,3) = HMvsTRAJ_sort_plot(i90,5);

    [HMvsTRAJ_sort,ind] = sort(HMvsTRAJ_plot(:,6));
    HMvsTRAJ_sort_plot = HMvsTRAJ_plot(ind,:);
    subplot(4,4,4)
    scatter(HMvsTRAJ_sort_plot(:,1), HMvsTRAJ_sort_plot(:,6),'b')
    hold on
    xaxi = [min(HMvsTRAJ_sort_plot(:,1)) max(HMvsTRAJ_sort_plot(:,1))];
    plot(xaxi,HMvsTRAJ_sort_plot(i99,6)*ones(1,2),'y','LineWidth',2.0)
    plot(xaxi,HMvsTRAJ_sort_plot(i98,6)*ones(1,2),'r','LineWidth',2.0)
    plot(xaxi,HMvsTRAJ_sort_plot(i95,6)*ones(1,2),'g','LineWidth',2.0)
    plot(xaxi,HMvsTRAJ_sort_plot(i90,6)*ones(1,2),'k','LineWidth',2.0)
    set(gca,'YScale','log')
    set(gca,'XScale','log')
    
    J_thresh_1(k,4) = HMvsTRAJ_sort_plot(i99,6);
    J_thresh_2(k,4) = HMvsTRAJ_sort_plot(i98,6);
    J_thresh_5(k,4) = HMvsTRAJ_sort_plot(i95,6);
    J_thresh_10(k,4) = HMvsTRAJ_sort_plot(i90,6);

    [HMvsTRAJ_sort,ind] = sort(HMvsTRAJ_plot(:,3));
    HMvsTRAJ_sort_plot = HMvsTRAJ_plot(ind,:);
    subplot(4,4,5)
    scatter(HMvsTRAJ_sort_plot(:,2), HMvsTRAJ_sort_plot(:,3),'b')
    hold on
    xaxi = [min(HMvsTRAJ_sort_plot(:,2)) max(HMvsTRAJ_sort_plot(:,2))];
    plot(xaxi,HMvsTRAJ_sort_plot(i99,3)*ones(1,2),'y','LineWidth',2.0)
    plot(xaxi,HMvsTRAJ_sort_plot(i98,3)*ones(1,2),'r','LineWidth',2.0)
    plot(xaxi,HMvsTRAJ_sort_plot(i95,3)*ones(1,2),'g','LineWidth',2.0)
    plot(xaxi,HMvsTRAJ_sort_plot(i90,3)*ones(1,2),'k','LineWidth',2.0)
    set(gca,'YScale','log')
    set(gca,'XScale','log')
    xlim([10e-15 10e-4])

    [HMvsTRAJ_sort,ind] = sort(HMvsTRAJ_plot(:,4));
    HMvsTRAJ_sort_plot = HMvsTRAJ_plot(ind,:);
    subplot(4,4,6)
    scatter(HMvsTRAJ_sort_plot(:,2), HMvsTRAJ_sort_plot(:,4),'b')
    hold on
    xaxi = [min(HMvsTRAJ_sort_plot(:,2)) max(HMvsTRAJ_sort_plot(:,2))];
    plot(xaxi,HMvsTRAJ_sort_plot(i99,4)*ones(1,2),'y','LineWidth',2.0)
    plot(xaxi,HMvsTRAJ_sort_plot(i98,4)*ones(1,2),'r','LineWidth',2.0)
    plot(xaxi,HMvsTRAJ_sort_plot(i95,4)*ones(1,2),'g','LineWidth',2.0)
    plot(xaxi,HMvsTRAJ_sort_plot(i90,4)*ones(1,2),'k','LineWidth',2.0)
    set(gca,'YScale','log')
    set(gca,'XScale','log')
    xlim([10e-15 10e-4])

    [HMvsTRAJ_sort,ind] = sort(HMvsTRAJ_plot(:,5));
    HMvsTRAJ_sort_plot = HMvsTRAJ_plot(ind,:);
    subplot(4,4,7)
    scatter(HMvsTRAJ_sort_plot(:,2), HMvsTRAJ_sort_plot(:,5),'b')
    hold on
    xaxi = [min(HMvsTRAJ_sort_plot(:,2)) max(HMvsTRAJ_sort_plot(:,2))];
    plot(xaxi,HMvsTRAJ_sort_plot(i99,5)*ones(1,2),'y','LineWidth',2.0)
    plot(xaxi,HMvsTRAJ_sort_plot(i98,5)*ones(1,2),'r','LineWidth',2.0)
    plot(xaxi,HMvsTRAJ_sort_plot(i95,5)*ones(1,2),'g','LineWidth',2.0)
    plot(xaxi,HMvsTRAJ_sort_plot(i90,5)*ones(1,2),'k','LineWidth',2.0)
    set(gca,'YScale','log')
    set(gca,'XScale','log')
    xlim([10e-15 10e-4])

    [HMvsTRAJ_sort,ind] = sort(HMvsTRAJ_plot(:,6));
    HMvsTRAJ_sort_plot = HMvsTRAJ_plot(ind,:);
    subplot(4,4,8)
    scatter(HMvsTRAJ_sort_plot(:,2), HMvsTRAJ_sort_plot(:,6),'b')
    hold on
    xaxi = [min(HMvsTRAJ_sort_plot(:,2)) max(HMvsTRAJ_sort_plot(:,2))];
    plot(xaxi,HMvsTRAJ_sort_plot(i99,6)*ones(1,2),'y','LineWidth',2.0)
    plot(xaxi,HMvsTRAJ_sort_plot(i98,6)*ones(1,2),'r','LineWidth',2.0)
    plot(xaxi,HMvsTRAJ_sort_plot(i95,6)*ones(1,2),'g','LineWidth',2.0)
    plot(xaxi,HMvsTRAJ_sort_plot(i90,6)*ones(1,2),'k','LineWidth',2.0)
    set(gca,'YScale','log')
    set(gca,'XScale','log')
    xlim([10e-15 10e-4])

    [HMvsTRAJ_sort,ind] = sort(HMvsTRAJ_plot(:,7));
    HMvsTRAJ_sort_plot = HMvsTRAJ_plot(ind,:);
    subplot(4,4,9)
    scatter(HMvsTRAJ_sort_plot(:,1), HMvsTRAJ_sort_plot(:,7),'b')
    hold on
    xaxi = [min(HMvsTRAJ_sort_plot(:,1)) max(HMvsTRAJ_sort_plot(:,1))];
    plot(xaxi,HMvsTRAJ_sort_plot(i99,7)*ones(1,2),'y','LineWidth',2.0)
    plot(xaxi,HMvsTRAJ_sort_plot(i98,7)*ones(1,2),'r','LineWidth',2.0)
    plot(xaxi,HMvsTRAJ_sort_plot(i95,7)*ones(1,2),'g','LineWidth',2.0)
    plot(xaxi,HMvsTRAJ_sort_plot(i90,7)*ones(1,2),'k','LineWidth',2.0)
    set(gca,'YScale','log')
    set(gca,'XScale','log')

    [HMvsTRAJ_sort,ind] = sort(HMvsTRAJ_plot(:,8));
    HMvsTRAJ_sort_plot = HMvsTRAJ_plot(ind,:);
    subplot(4,4,10)
    scatter(HMvsTRAJ_sort_plot(:,1), HMvsTRAJ_sort_plot(:,8),'b')
    hold on
    xaxi = [min(HMvsTRAJ_sort_plot(:,1)) max(HMvsTRAJ_sort_plot(:,1))];
    plot(xaxi,HMvsTRAJ_sort_plot(i99,8)*ones(1,2),'y','LineWidth',2.0)
    plot(xaxi,HMvsTRAJ_sort_plot(i98,8)*ones(1,2),'r','LineWidth',2.0)
    plot(xaxi,HMvsTRAJ_sort_plot(i95,8)*ones(1,2),'g','LineWidth',2.0)
    plot(xaxi,HMvsTRAJ_sort_plot(i90,8)*ones(1,2),'k','LineWidth',2.0)
    set(gca,'YScale','log')
    set(gca,'XScale','log')

    [HMvsTRAJ_sort,ind] = sort(HMvsTRAJ_plot(:,9));
    HMvsTRAJ_sort_plot = HMvsTRAJ_plot(ind,:);
    subplot(4,4,11)
    scatter(HMvsTRAJ_sort_plot(:,1), HMvsTRAJ_sort_plot(:,9),'b')
    hold on
    xaxi = [min(HMvsTRAJ_sort_plot(:,1)) max(HMvsTRAJ_sort_plot(:,1))];
    plot(xaxi,HMvsTRAJ_sort_plot(i99,9)*ones(1,2),'y','LineWidth',2.0)
    plot(xaxi,HMvsTRAJ_sort_plot(i98,9)*ones(1,2),'r','LineWidth',2.0)
    plot(xaxi,HMvsTRAJ_sort_plot(i95,9)*ones(1,2),'g','LineWidth',2.0)
    plot(xaxi,HMvsTRAJ_sort_plot(i90,9)*ones(1,2),'k','LineWidth',2.0)
%     set(gca,'YScale','log')
    set(gca,'XScale','log')

    [HMvsTRAJ_sort,ind] = sort(HMvsTRAJ_plot(:,10));
    HMvsTRAJ_sort_plot = HMvsTRAJ_plot(ind,:);
    subplot(4,4,12)
    scatter(HMvsTRAJ_sort_plot(:,1), HMvsTRAJ_sort_plot(:,10),'b')
    hold on
    xaxi = [min(HMvsTRAJ_sort_plot(:,1)) max(HMvsTRAJ_sort_plot(:,1))];
    plot(xaxi,HMvsTRAJ_sort_plot(i99,10)*ones(1,2),'y','LineWidth',2.0)
    plot(xaxi,HMvsTRAJ_sort_plot(i98,10)*ones(1,2),'r','LineWidth',2.0)
    plot(xaxi,HMvsTRAJ_sort_plot(i95,10)*ones(1,2),'g','LineWidth',2.0)
    plot(xaxi,HMvsTRAJ_sort_plot(i90,10)*ones(1,2),'k','LineWidth',2.0)
    set(gca,'YScale','log')
    set(gca,'XScale','log')

    [HMvsTRAJ_sort,ind] = sort(HMvsTRAJ_plot(:,7));
    HMvsTRAJ_sort_plot = HMvsTRAJ_plot(ind,:);
    subplot(4,4,13)
    scatter(HMvsTRAJ_sort_plot(:,2), HMvsTRAJ_sort_plot(:,7),'b')
    hold on
    xaxi = [min(HMvsTRAJ_sort_plot(:,2)) max(HMvsTRAJ_sort_plot(:,2))];
    plot(xaxi,HMvsTRAJ_sort_plot(i99,7)*ones(1,2),'y','LineWidth',2.0)
    plot(xaxi,HMvsTRAJ_sort_plot(i98,7)*ones(1,2),'r','LineWidth',2.0)
    plot(xaxi,HMvsTRAJ_sort_plot(i95,7)*ones(1,2),'g','LineWidth',2.0)
    plot(xaxi,HMvsTRAJ_sort_plot(i90,7)*ones(1,2),'k','LineWidth',2.0)
    set(gca,'YScale','log')
    set(gca,'XScale','log')
    xlim([10e-15 10e-4])

    [HMvsTRAJ_sort,ind] = sort(HMvsTRAJ_plot(:,8));
    HMvsTRAJ_sort_plot = HMvsTRAJ_plot(ind,:);
    subplot(4,4,14)
    scatter(HMvsTRAJ_sort_plot(:,2), HMvsTRAJ_sort_plot(:,8),'b')
    hold on
    xaxi = [min(HMvsTRAJ_sort_plot(:,2)) max(HMvsTRAJ_sort_plot(:,2))];
    plot(xaxi,HMvsTRAJ_sort_plot(i99,8)*ones(1,2),'y','LineWidth',2.0)
    plot(xaxi,HMvsTRAJ_sort_plot(i98,8)*ones(1,2),'r','LineWidth',2.0)
    plot(xaxi,HMvsTRAJ_sort_plot(i95,8)*ones(1,2),'g','LineWidth',2.0)
    plot(xaxi,HMvsTRAJ_sort_plot(i90,8)*ones(1,2),'k','LineWidth',2.0)
    set(gca,'YScale','log')
    set(gca,'XScale','log')
    xlim([10e-15 10e-4])

    [HMvsTRAJ_sort,ind] = sort(HMvsTRAJ_plot(:,9));
    HMvsTRAJ_sort_plot = HMvsTRAJ_plot(ind,:);
    subplot(4,4,15)
    scatter(HMvsTRAJ_sort_plot(:,2), HMvsTRAJ_sort_plot(:,9),'b')
    hold on
    xaxi = [min(HMvsTRAJ_sort_plot(:,2)) max(HMvsTRAJ_sort_plot(:,2))];
    plot(xaxi,HMvsTRAJ_sort_plot(i99,9)*ones(1,2),'y','LineWidth',2.0)
    plot(xaxi,HMvsTRAJ_sort_plot(i98,9)*ones(1,2),'r','LineWidth',2.0)
    plot(xaxi,HMvsTRAJ_sort_plot(i95,9)*ones(1,2),'g','LineWidth',2.0)
    plot(xaxi,HMvsTRAJ_sort_plot(i90,9)*ones(1,2),'k','LineWidth',2.0)
%     set(gca,'YScale','log')
    set(gca,'XScale','log')
    xlim([10e-15 10e-4])

    [HMvsTRAJ_sort,ind] = sort(HMvsTRAJ_plot(:,10));
    HMvsTRAJ_sort_plot = HMvsTRAJ_plot(ind,:);
    subplot(4,4,16)
    scatter(HMvsTRAJ_sort_plot(:,2), HMvsTRAJ_sort_plot(:,10),'b')
    hold on
    xaxi = [min(HMvsTRAJ_sort_plot(:,2)) max(HMvsTRAJ_sort_plot(:,2))];
    plot(xaxi,HMvsTRAJ_sort_plot(i99,10)*ones(1,2),'y','LineWidth',2.0)
    plot(xaxi,HMvsTRAJ_sort_plot(i98,10)*ones(1,2),'r','LineWidth',2.0)
    plot(xaxi,HMvsTRAJ_sort_plot(i95,10)*ones(1,2),'g','LineWidth',2.0)
    plot(xaxi,HMvsTRAJ_sort_plot(i90,10)*ones(1,2),'k','LineWidth',2.0)
    set(gca,'YScale','log')
    set(gca,'XScale','log')
    xlim([10e-15 10e-4])
    
end
