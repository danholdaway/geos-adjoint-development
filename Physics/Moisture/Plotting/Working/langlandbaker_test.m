close all
clear
clc


xt = [2.567 3.456];
xf = [2.532 3.352];
xg = [2.347 3.532];

C = [1.23 4.67; 2.43 0.97];

xft = xf - xt;
xgt = xg - xt;

xfg = xf - xg;

Cxft = C*xft';
Cxgt = C*xgt';

ef = xft*C*xft';
eg = xgt*C*xgt';

defg = ef - eg


dJfdxf = C*xft';
dJgdxg = C*xgt';

dJdx_sum = dJfdxf + dJfdxf;

defg_lin = xfg*dJdx_sum

