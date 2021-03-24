close all
clear
clc

t = -100:0.1:100;

f = t.^2;


a1 = 1.0;

f_tl = a1^2 + (t-a1)*2*a1;

t1 = 1.0;
t2 = 2.0;
t3 = 3.0;

f_tl_t1 = a1^2 + (t1-a1)*2*a1

f_tl_t2 = a1^2 + (t2-t1)*2*a1

a2 = 2.0;

f_tl_2 = a2^2 + (t-a2)*2*a2;

f_tl_t3 = a2^2 + (t3-t2)*2*a2

line_wid = 1.5;
mark_siz = 10.0;
fontsize = 14;


plot(t,f, 'LineWidth', line_wid)
hold on
plot(t,f_tl,'r', 'LineWidth', line_wid)
plot(t,f_tl_2,'g', 'LineWidth', line_wid)

plot(t1,f_tl_t1,'*k','MarkerSize',mark_siz, 'LineWidth', line_wid);
plot(t2,f_tl_t2,'*k','MarkerSize',mark_siz, 'LineWidth', line_wid);
plot(t3,f_tl_t3,'*k','MarkerSize',mark_siz, 'LineWidth', line_wid);

xlim([-2.0 5])
ylim([-1.5 12])



ylabel('f(x)','FontSize',fontsize,'FontName','TimesNewRoman')
xlabel('x','FontSize',fontsize,'FontName','TimesNewRoman')

legend('f = x^2','f_{tl} = a^2+(x-a)2a','f_{tl} = a^2+(x-a)2a','Linear Solution','Location','NorthWest')


set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
grid on
box on