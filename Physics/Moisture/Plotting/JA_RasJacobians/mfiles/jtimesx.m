function [J] = jtimesx(J,x)


Nk = length(J);

for i = 1:Nk
    
    J(i,:) = J(i,:).*x';
    
end