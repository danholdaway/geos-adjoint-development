function Y = LUTriDiagSolveArray(LM,A,B,C,Y,PHASE)

if (PHASE == 1)

    %Sweep down, modifying rhs with multiplier A
    for l=2:LM
      Y(:,:,l) = Y(:,:,l) - A(:,:,l)*Y(:,:,l-1);
    end

    %Sweep up, solving for updated value. Note B has the inverse of the main diagonal
    Y(:,:,LM) = B(:,:,LM-1)*Y(:,:,LM)/(B(:,:,LM-1)-A(:,:,LM)*(1.0+C(:,:,LM-1)*B(:,:,LM-1)));
    for l=LM-1:-1:1
      Y(:,:,l) = B(:,:,l)*(Y(:,:,l)-C(:,:,l)*Y(:,:,l+1) );
    end

elseif (PHASE == 2)

    %Adjoint of sweep up, solving for updated value. Note B has the inverse of the main diagonal
    for l=1:LM-1
      TMP = B(:,:,l)*Y(:,:,l);
      Y(:,:,l+1) = Y(:,:,l+1) - C(:,:,l)*TMP;
      Y(:,:,l) = TMP;
    end
    Y(:,:,LM) = B(:,:,LM-1)*Y(:,:,LM)/(B(:,:,LM-1)-A(:,:,LM)*(C(:,:,LM-1)*B(:,:,LM-1)+1.0));

    %Adjoint of sweep down, modifying rhs with multiplier A
    for l=LM:-1:2
      Y(:,:,l-1) = Y(:,:,l-1) - A(:,:,l)*Y(:,:,l);
    end

end




    
    
end