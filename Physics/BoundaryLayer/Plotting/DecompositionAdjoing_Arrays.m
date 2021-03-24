close all
clear
clc
mydir = pwd;

scale = 100;

% Create random arrays
im = 5;
jm = 5;
lm = 10;

A = rand(im,jm,lm)/scale;
B = rand(im,jm,lm)/scale;
C = rand(im,jm,lm)/scale;
Xp = rand(im,jm,lm)/scale;
Yhat = rand(im,jm,lm)/scale;

A_LU = zeros(im,jm,lm);
B_LU = zeros(im,jm,lm);
C_LU = zeros(im,jm,lm);

Yp1 = zeros(im,jm,lm);
Xhat1 = zeros(im,jm,lm);

for i = 1:im
    for j = 1:jm
        
        [A_LU(i,j,:),B_LU(i,j,:),C_LU(i,j,:)] = LUTriDiag(lm,A(i,j,:),B(i,j,:),C(i,j,:));

        Yp1(i,j,:) = LUTriDiagSolve(lm,A_LU(i,j,:),B_LU(i,j,:),C_LU(i,j,:),Xp(i,j,:),1);

        Xhat1(i,j,:) = LUTriDiagSolve(lm,A_LU(i,j,:),B_LU(i,j,:),C_LU(i,j,:),Yhat(i,j,:),2);
        
        
    end
end

dot1 = dot(Yp1(:),Yhat(:));
dot2 = dot(Xp(:),Xhat1(:));
test1 = abs(dot1-dot2)/dot2;

disp([dot1 dot2])
disp(test1)



[A_LU1,B_LU1] = LUTriDiagArray(lm,A,B,C);
C_LU1 = C;

Yp2   = LUTriDiagSolve(lm,A_LU1,B_LU1,C_LU1,Xp,1);
Xhat2 = LUTriDiagSolve(lm,A_LU1,B_LU1,C_LU1,Yhat,2);


dot1 = dot(Yp2(:),Yhat(:));
dot2 = dot(Xp(:),Xhat2(:));
test1 = abs(dot1-dot2)/dot2;

disp([dot1 dot2])
disp(test1)


a(i,:) = [   539786.704904748        539786.704904750    4.744720566329065E-015 ]; i = i + 1;











