function [x0, lambda] = feasible_sp(A, Q, q)

    sz = size(A);
    n = sz(2);
    m = sz(1);
    
    lambda = struct('eqlin', rand(m, 1), 'lower', ones(n, 1));
   
	sz = size(A);
    n = sz(2);
    m = sz(1);
    
    % We init lambda.eqlin as rand positive
    lambda = struct('eqlin', zeros(m, 1), 'lower', rand(n, 1));
    
    x0 = ones(n, 1);
    % To initialize x0 we exploit the structure of the problem 
    % to compute a feasible starting point. First we compute a new matrix
    % 
    % A_new(i,j) = A(i,j) * 1/sum(A(i,1:end))
    % 
    % The we choose x0 as the maximum along all columns. If x0 has a 0
    % component we set it to a small positive number to avoid future
    % divisions by zero.
    tmpA = (1./sum(A, 2)).*A;
    x0 = sum(tmpA, 1)';
    idx = (x0 == 0);
    if any(idx)
    x0(idx) = 1e-6;
    end
    
    % Initializig lambda.eqlin such that solves gradient lagragian = 0
    % and then we adjust the negative components.
    % nocedal - wrigth
    lambda.lower = 2*Q*x0 + q + A'*lambda.eqlin;
    
    ind = lambda.lower < 0;
    if any( ind )
        dl = max(0,-(3/2)*min(lambda.lower));
    end
    
    lambda.lower = dl + lambda.lower;    
end