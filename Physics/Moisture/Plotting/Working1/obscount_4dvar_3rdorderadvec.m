close all
clear
clc



i = 1;
x(i) = 20;    a(i) = 3219025   ; b(i) = 3219025;    c(i) = 3219025; i = i+1;
x(i) = 20.25; a(i) = 2961314   ; b(i) = 2961343;    c(i) = 2961343; i = i+1;
x(i) = 20.5;  a(i) = 3085336   ; b(i) = 3085618;    c(i) = 3095092; i = i+1;
x(i) = 20.75; a(i) = 3096561   ; b(i) = 3098262;    c(i) = 3111183; i = i+1;
x(i) = 21;    a(i) = 3174215   ; b(i) = 3175656;    c(i) = 3191605; i = i+1;
x(i) = 21.25; a(i) = 2939520   ; b(i) = 2939946;    c(i) = 2956459; i = i+1;
x(i) = 21.5;  a(i) = 3086038   ; b(i) = 3086793;    c(i) = 3110653; i = i+1;
x(i) = 21.75; a(i) = 3110970   ; b(i) = 3114131;    c(i) = 3136787; i = i+1;
x(i) = 22;    a(i) = 3229432   ; b(i) = 3230086;    c(i) = 3255491; i = i+1;
x(i) = 22.25; a(i) = 3065649   ; b(i) = 3064769;    c(i) = 3081535; i = i+1;
x(i) = 22.5;  a(i) = 3129013   ; b(i) = 3129170;    c(i) = 3147447; i = i+1;
x(i) = 22.75; a(i) = 3149570   ; b(i) = 3148351;    c(i) = 3168014; i = i+1;
x(i) = 23;    a(i) = 3265296   ; b(i) = 3264072;    c(i) = 3281872; i = i+1;
x(i) = 23.25; a(i) = 2999913   ; b(i) = 3001616;    c(i) = 3017540; i = i+1;
x(i) = 23.5;  a(i) = 3134767   ; b(i) = 3136830;    c(i) = 3156004; i = i+1;
x(i) = 23.75; a(i) = 3095369   ; b(i) = 3096566;    c(i) = 3115248; i = i+1;
x(i) = 24;    a(i) = 3225846   ; b(i) = 3227230;    c(i) = 3248791; i = i+1;
x(i) = 24.25; a(i) = 2934684   ; b(i) = 2935356;    c(i) = 2955187; i = i+1;
x(i) = 24.5;  a(i) = 3153354   ; b(i) = 3153853;    c(i) = 3178486; i = i+1;
x(i) = 24.75; a(i) = 2849774   ; b(i) = 2848251;    c(i) = 2867570; i = i+1;
x(i) = 25;    a(i) = 3237650   ; b(i) = 3235551;    c(i) = 3256428; i = i+1;
x(i) = 25.25; a(i) = 2986586   ; b(i) = 2986929;    c(i) = 3009360; i = i+1;
x(i) = 25.5;  a(i) = 3130392   ; b(i) = 3129343;    c(i) = 3154272; i = i+1;
% x(i) = 25.75; a(i) = 3118338   ; b(i) = 3123176;    c(i) = ; i = i+1;
% x(i) = 26;    a(i) = 3146712   ; b(i) = 3147632;    c(i) = ; i = i+1;
% x(i) = 26.25; a(i) = 2866421   ; b(i) = 2870830;    c(i) = ; i = i+1;
% x(i) = 26.5;  a(i) = 3103156   ; b(i) = 3104044;    c(i) = ; i = i+1;
% x(i) = 26.75; a(i) = 3083895   ; b(i) = 3086251;    c(i) = ; i = i+1;
% x(i) = 27;    a(i) = 3131647   ; b(i) = 3132470;    c(i) = ; i = i+1;
% x(i) = 27.25; a(i) = 2934911   ; b(i) = 2933747;    c(i) = ; i = i+1;
% x(i) = 27.5;  a(i) = 3066691   ; b(i) = 3068197; i = i+1;
% x(i) = 27.75; a(i) = 3044235   ; b(i) = 3042200; i = i+1;
% x(i) = 28;    a(i) = 3146481   ; b(i) = 3148380; i = i+1;
% x(i) = 28.25; a(i) = 2940620   ; b(i) = 2942022; i = i+1;
% x(i) = 28.5;  a(i) = 3065610   ; b(i) = 3064646; i = i+1;
% x(i) = 28.75; a(i) = 3006212   ; b(i) = 3009858; i = i+1;
% x(i) = 29;    a(i) = 3203042   ; b(i) = 3204058; i = i+1;
% x(i) = 29.25; a(i) = 3026142   ; b(i) = 3028210; i = i+1;
% x(i) = 29.5;  a(i) = 3104615   ; b(i) = ; i = i+1;
% x(i) = 29.75; a(i) = 3045131   ; b(i) = ; i = i+1;
% x(i) = 30;    a(i) = 3205132   ; b(i) = ; i = i+1;
% x(i) = 30.25; a(i) = 2957357   ; b(i) = ; i = i+1;
% x(i) = 30.5;  a(i) = 3119260   ; b(i) = ; i = i+1;
% x(i) = 30.75; a(i) = 3026689   ; b(i) = ; i = i+1;
% x(i) = 31;    a(i) = 3111242   ; b(i) = ; i = i+1;
% x(i) = 31.25; a(i) = 2928331   ; b(i) = ; i = i+1;
% x(i) = 31.5;  a(i) = 3105011   ; b(i) = ; i = i+1;
% x(i) = 31.75; a(i) = 3086901   ; b(i) = ; i = i+1;
% x(i) = 32;    a(i) = 3226349   ; b(i) = ; i = i+1;
% x(i) = 32.25; a(i) = 2326280   ; b(i) = ; i = i+1;
% x(i) = 32.5;  a(i) = 3102735   ; b(i) = ; i = i+1;
% x(i) = 32.75; a(i) = 3107563   ; b(i) = ; i = i+1;
% x(i) = 33;    a(i) = 3165872   ; b(i) = ; i = i+1;
% x(i) = 33.25; a(i) = 2939291   ; b(i) = ; i = i+1;
% x(i) = 33.5;  a(i) = 3084509   ; b(i) = ; i = i+1;
% x(i) = 33.75; a(i) = 3048784   ; b(i) = ; i = i+1;
% x(i) = 34;    a(i) = 3124925   ; b(i) = ; i = i+1;


fontsize = 13;
lin_wid = 1.5;
j = i-1;

figure
set(gcf,'position',[528   377   751   542])
plot(x(1:j),b(1:j)-a(1:j),'b','LineWidth',lin_wid)
hold on
plot(x(1:j),c(1:j)-b(1:j),'r','LineWidth',lin_wid)
hold on
% plot(x(1:4:j),b(1:4:j)-a(1:4:j),'or','LineWidth',lin_wid)
box on
xlabel('Assimilation Time (2014)','FontSize',fontsize,'FontName','TimesNewRoman')
ylabel('Moist - Dry','FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
xlim([20 34])
% ylim([0 26000])
% set(gca,'xTick',20:2:34)
% set(gca,'xTickLabel',{'20-Apr','22-Apr','24-Apr','26-Apr','28-Apr','30-Apr','2-May','4-May'})
% set(gca,'YTick',[0,2000,4000,6000,8000,10000,12000,14000,16000,18000,20000,22000,24000,26000'])
% set(gca,'yTickLabel',{'0','2000','4000','6000','8000','10000','12000','14000','16000','18000','20000','22000','24000','26000'})
legend('All times','0000 UTC times','Location','NorthWest')

legend boxoff