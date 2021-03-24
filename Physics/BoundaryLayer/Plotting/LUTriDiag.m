function [A,B,C] = LUTriDiag(LM,A,B,C)

B(1) = 1.0 / B(1);
for L = 2:LM
   A(L) = A(L) * B(L-1);
   B(L) = 1.0 ./ ( B(L) - C(L-1) * A(L) );
end

