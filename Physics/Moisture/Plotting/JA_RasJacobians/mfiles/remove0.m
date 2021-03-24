function [J] = remove0(J,percent_cut)

max_J = max(max(J));
min_J = min(min(J));

Nk = length(J);

for i = 1:Nk
    for j = 1:Nk
        
        if J(i,j) > 0  && J(i,j) < percent_cut*max_J
            J(i,j) = 0.0;
        elseif J(i,j) < 0  && J(i,j) > percent_cut*min_J
            J(i,j) = 0.0;
        end
    end
end