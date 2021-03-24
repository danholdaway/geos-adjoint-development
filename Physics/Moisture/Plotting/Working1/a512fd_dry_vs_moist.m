close all
clear
clc
fontsize = 12;
lin_wid = 1.6;

i = 1;
x(i) = 20;    a(i) = 3219025   ; b(i) = 3219025; i = i+1;
x(i) = 20.25; a(i) = 2961314   ; b(i) = 2964298; i = i+1;
x(i) = 20.5;  a(i) = 3085336   ; b(i) = 3092398; i = i+1;
x(i) = 20.75; a(i) = 3096561   ; b(i) = 3103614; i = i+1;
x(i) = 21;    a(i) = 3174215   ; b(i) = 3184752; i = i+1;
x(i) = 21.25; a(i) = 2939520   ; b(i) = 2946366; i = i+1;
x(i) = 21.5;  a(i) = 3086038   ; b(i) = 3096465; i = i+1;
x(i) = 21.75; a(i) = 3110970   ; b(i) = 3124840; i = i+1;
x(i) = 22;    a(i) = 3229432   ; b(i) = 3242457; i = i+1;
x(i) = 22.25; a(i) = 3065649   ; b(i) = 3075305; i = i+1;
x(i) = 22.5;  a(i) = 3129013   ; b(i) = 3145261; i = i+1;
x(i) = 22.75; a(i) = 3149570   ; b(i) = 3163991; i = i+1;
x(i) = 23;    a(i) = 3265296   ; b(i) = 3278374; i = i+1;
x(i) = 23.25; a(i) = 2999913   ; b(i) = 3012751; i = i+1;
x(i) = 23.5;  a(i) = 3134767   ; b(i) = 3150575; i = i+1;
x(i) = 23.75; a(i) = 3095369   ; b(i) = 3109790; i = i+1;
x(i) = 24;    a(i) = 3225846   ; b(i) = 3242051; i = i+1;
x(i) = 24.25; a(i) = 2934684   ; b(i) = 2949058; i = i+1;
x(i) = 24.5;  a(i) = 3153354   ; b(i) = 3171400; i = i+1;
x(i) = 24.75; a(i) = 2849774   ; b(i) = 2864300; i = i+1;
x(i) = 25;    a(i) = 3237650   ; b(i) = 3252215; i = i+1;
x(i) = 25.25; a(i) = 2986586   ; b(i) = 3005237; i = i+1;
x(i) = 25.5;  a(i) = 3130392   ; b(i) = 3147968; i = i+1;
x(i) = 25.75; a(i) = 3118338   ; b(i) = 3143402; i = i+1;
x(i) = 26;    a(i) = 3146712   ; b(i) = 3163021; i = i+1;
x(i) = 26.25; a(i) = 2866421   ; b(i) = 2880404; i = i+1;
x(i) = 26.5;  a(i) = 3103156   ; b(i) = 3113429; i = i+1;
x(i) = 26.75; a(i) = 3083895   ; b(i) = 3101567; i = i+1;
x(i) = 27;    a(i) = 3131647   ; b(i) = 3148050; i = i+1;
x(i) = 27.25; a(i) = 2934911   ; b(i) = 2947996; i = i+1;
x(i) = 27.5;  a(i) = 3066691   ; b(i) = 3083487; i = i+1;
x(i) = 27.75; a(i) = 3044235   ; b(i) = 3063696; i = i+1;
x(i) = 28;    a(i) = 3146481   ; b(i) = 3163194; i = i+1;
x(i) = 28.25; a(i) = 2940620   ; b(i) = 2952317; i = i+1;
x(i) = 28.5;  a(i) = 3065610   ; b(i) = 3080296; i = i+1;
x(i) = 28.75; a(i) = 3006212   ; b(i) = 3024286; i = i+1;
x(i) = 29;    a(i) = 3203042   ; b(i) = 3219275; i = i+1;
% x(i) = 29.25; a(i) =    ; b(i) = 3040218; i = i+1;
% x(i) = 29.5;  a(i) =    ; b(i) = 3119212; i = i+1;
% x(i) = 29.75; a(i) =    ; b(i) = 3055639; i = i+1;
% x(i) = 30;    a(i) =    ; b(i) = 3223967; i = i+1;
% x(i) = 30.25; a(i) =    ; b(i) = 2966397; i = i+1;
% x(i) = 30.5;  a(i) =    ; b(i) = 3133390; i = i+1;
% x(i) = 30.75; a(i) =    ; b(i) = 3039555; i = i+1;


figure
set(gcf,'position',[409   499   870   420])
plot(x,a,'-xb')
hold on
plot(x,b,'-xr')
box on
xlabel('Assimilation Time','FontSize',fontsize,'FontName','TimesNewRoman')
ylabel('Number of observations','FontSize',fontsize,'FontName','TimesNewRoman')
title('Number of observations assimilated using 4DVAR','FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
xlim([20 26])
ylim([2.8e6 3.3e6])
set(gca,'xTickLabel',{'20-Apr-14','21-Apr-14','22-Apr-14','23-Apr-14','24-Apr-14','25-Apr-14','26-Apr-14'})
legend('Dry Physics','Moist Physics','Location','SouthWest')

j = i-1;

figure
set(gcf,'position',[409   499   870   420])
plot(x(1:j),b(1:j)-a(1:j),'-x')
hold on
plot(x(1:4:j),b(1:4:j)-a(1:4:j),'or')
box on
xlabel('Assimilation Time','FontSize',fontsize,'FontName','TimesNewRoman')
ylabel('Moist - Dry','FontSize',fontsize,'FontName','TimesNewRoman')
title('Number of observations gained using moist physics','FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
% xlim([20 26])
% ylim([0 26000])
% set(gca,'xTickLabel',{'20-Apr-14','21-Apr-14','22-Apr-14','23-Apr-14','24-Apr-14','25-Apr-14','26-Apr-14'})
% set(gca,'YTick',[0,2000,4000,6000,8000,10000,12000,14000,16000,18000,20000,22000,24000,26000'])
% set(gca,'yTickLabel',{'0','2000','4000','6000','8000','10000','12000','14000','16000','18000','20000','22000','24000','26000'})
legend('All times','0000 UTC times','Location','NorthWest')

figure
set(gcf,'position',[409   499   870   420])
plot(x(1:j),(b(1:j)-a(1:j))./a(1:j),'-x')
hold on
plot(x(1:4:j),(b(1:4:j)-a(1:4:j))./a(1:4:j),'or')
box on
xlabel('Assimilation Time','FontSize',fontsize,'FontName','TimesNewRoman')
ylabel('Moist - Dry','FontSize',fontsize,'FontName','TimesNewRoman')
title('Number of observations gained using moist physics','FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
xlim([20 26])
set(gca,'xTickLabel',{'20-Apr-14','21-Apr-14','22-Apr-14','23-Apr-14','24-Apr-14','25-Apr-14','26-Apr-14'})
legend('All times','0000 UTC times','Location','NorthWest')