close all
clear
clc

steps = 33;

divs = steps/3;

x = 1:steps;

b = zeros(1,length(x));
g = zeros(1,length(x));
r = zeros(1,length(x));

b(0*divs+1:1*divs) = 1.0;
g(0*divs+1:1*divs) = fliplr(1-x(1:divs)/divs);
r(0*divs+1:1*divs) = 0.0;

b(1*divs+1:2*divs) = 1-x(1:divs)/divs;
g(1*divs+1:2*divs) = 1.0;
r(1*divs+1:2*divs) = x(1:divs)/divs;

b(2*divs+1:3*divs) = 0.0;
g(2*divs+1:3*divs) = 1-x(1:divs)/divs;
r(2*divs+1:3*divs) = 1.0;

r(end) = 0.90;
b = fliplr(r);

figure
hold on

plot(x,r,'r')
plot(x,g,'g')
plot(x,b,'b')


cmap = zeros(length(x),3);

cmap(:,1) = r;
cmap(:,2) = g;
cmap(:,3) = b;


cmap1 = colormap;
r1 = cmap1(:,1);
g1 = cmap1(:,2);
b1 = cmap1(:,3);

x1 = 1:length(r1);


figure
hold on

plot(x1,r1,'r')
plot(x1,g1,'g')
plot(x1,b1,'b')
