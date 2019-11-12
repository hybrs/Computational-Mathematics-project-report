function [x0, lambda] = feasible_sp(A, Q, q)

    sz = size(A);
    n = sz(2);
    m = sz(1);
    
    lambda = struct('eqlin', rand(m, 1), 'lower', rand(n, 1));
    x0 = (1/n)*ones(n, 1);
end
    %initializing x0
    
    %tmpA = (1./sum(A, 2)).*A;
    
    %for c=1:n
    %    cel = tmpA(1:end, c);
    %    for jj=1:size(cel)
    %    if not(cel(jj) == 0) && cel(jj) < x0(c)
    %        x0(c) = cel(jj);
    %    end
    %    end
    %end
    
    %Check for x_i  that should be 0 (0-column)
    %ind = sum(A, 1) == 0;
    %if any( ind )
    %    x0(ind) = 1e-6*rand;
    %end
    
    %initializing lambda_s (solve gradient lagragian  = 0)
    %lambda.lower = 2*Q*x0 + q +A'*lambda.eqlin;
    
    
    
%end