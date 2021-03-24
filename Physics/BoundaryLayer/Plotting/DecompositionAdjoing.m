close all
clear
clc
mydir = pwd;


% Create random A
im = 1000;
A = zeros(im,im);
for i = 1:im
    A(i,i) = rand;
    if i > 1
        A(i-1,i) = rand;
        A(i,i-1) = rand;
    end
end
Xp = rand(im,1);
Yhat = rand(im,1);

A1 = [0; diag(A,-1)];
B1 = diag(A);
C1 = [diag(A,1); 0];


%LU
[A1,B1,C1] = LUTriDiag(im,A1,B1,C1);

%Tangent linear
Yp = LUTriDiagSolve(im,A1,B1,C1,Xp,1);

%Adjoint
Xhat = LUTriDiagSolve(im,A1,B1,C1,Yhat,2);


dot1 = Yp'*Yhat
dot2 = Xp'*Xhat

abs(dot1-dot2)/dot2














