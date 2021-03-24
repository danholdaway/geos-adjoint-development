function Y = LUTriDiagSolve(im,A,B,C,Y,Phase)


if Phase == 1

    % Sweep down, modifying rhs with multiplier A
    for l=2:im
        Y(l) = Y(l) - A(l)*Y(l-1);
    end

    % Sweep up, solving for updated value. Note B has the inverse of the main diagonal

    Y(im) = B(im-1)*Y(im)/(B(im-1)-A(im)*(1.0+C(im-1)*B(im-1)));
    for l=im-1:-1:1
        Y(l) = B(l)*(Y(l)-C(l)*Y(l+1) );
    end
    
elseif Phase == 2

    % Adjoint of sweep up, solving for updated value. Note B has the inverse of the main diagonal
    for l=1:im-1
        TMP = B(l)*Y(l);
        Y(l+1) = Y(l+1) - C(l)*TMP;
        Y(l) = TMP;
    end
    Y(im) = B(im-1)*Y(im)/(B(im-1)-A(im)*(C(im-1)*B(im-1)+1.0));

    % Adjoint of sweep down, modifying rhs with multiplier A
    for l=im:-1:2
        Y(l-1) = Y(l-1) - A(l)*Y(l);
    end
    
    
end