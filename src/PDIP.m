function [x, fval, lambda, exit_code, time, it, gaps] =  PDIP(problem, maxitt, eps, solver)
%function [x, fval, lambda, exit_code, time, it, gaps] = PDIP(problem, maxitt, eps, solver)
%
%  Solves a CQP problem P with a primal-dual interior-point method
%   
% (P) = { min x'Qx + q'x | Ax = b, x >= 0}
% 
% Input:
%
% - problem (struct): the CQP problem
%
% - maxit (integer): maximum number of iteration allowed
%
% - eps (real, default 1e-15): complementary gap termination threshold
%
% Output: the solution, with the following fields:
%
% - x(n \times 1  real vector): the final coordinate of x
%
% - fval(scalar): the value of the function in x, 
%
% - lambda(struct): contains the lagrangian multiplicators for equalities 
%                   and lower bounds, respectively in lambda.eqlin and lambda.lower
% 
% - exit_code(integer): 1 if the gap reduced under the user-def threshold eps, 
%                       2 if the number of iteration exceeded maxit
%
% - time(real): time in seconds
%
% - it(integer) :total number of iterations
% 
% - gaps(it \times 1 vector of real): vector containing the gap for each iteration
%
warning("off", "all");

Q = problem.Q;
q = problem.q;

n = size(Q);
n = n(1);

m = problem.m;

A = problem.A;
b = problem.b;

k = size(A);
k = k(1);


x = []; fval = inf; lambda = []; exit_code = -1; time = 0; it = 0; gap = inf;

% Input checks
%
% Check the number of constraints in A   
if m > round(n/2)
    fprintf("Too many constraints in A, please retry.");
    return;
end
% Check for the correctnes of the structure of A
if any(sum(A, 2) < 2) || any(not(sum(A, 1) == 1)) || not(rank(A) == k)
    fprintf("A is not a set of k-disjoint simplices, please retry")
    exit_code = -2;
end
% Check for positive-semidefinite Q
eigQ = eig(Q);
if any(eigQ < 0)
    fprintf( 'Q is not positive-semidef, please retry.\n');
    exit_code = -3;
end

tic;

%init the starting point (x0, lambda_eq0, lambda_s0)
[x0, lambda] = feasible_sp(A, Q, q);
lambda_eq = lambda.eqlin;
lambda_s = lambda.lower;

% init mu and sigma
alpha = 0.1;
mu = (x0'*lambda_s)/n;
sigma = 0.3;

fprintf( 'Primal-Dual Interior Point method\n');
fprintf( '\titer\t\tfval\t\tgap\t\t\t  dualf\t\tprimalf\t\t  sTx\n' );
fprintf('-----------------------------------------------------------------------------\n');

x = x0;
it = 0;
gap = inf;
rd = inf;

gaps = []; residuals = struct('rd', [], 'rp', []);

while true
    it=it+1;
    Qx = Q * x;
    qx = q'*x;
    fval = x'*Qx + qx;
    gradL = 2*Qx + q +A'*lambda_eq - lambda_s;
    
    primal = fval;
    dual = fval + lambda_eq' *(A*x -b) - lambda_s'*x ;
    
    gapprec = gap;
    gap = ( primal - dual ) / max( abs( primal ) , 1 );
    
    rd = norm(gradL, 2);
    rp = norm(A*x - b, 2);
    rxs = norm(lambda_s'*x, 2);
    
    gaps = [gaps; gap]; 
    %residuals.rd = [residuals.rd; rd]; 
    %residuals.rp = [residuals.rp; rp];
    
    fprintf( '%4d\t\t%1.3e\t%1.3e\t\t%1.3e\t%1.3e\t%1.3e\n' , it , fval, gap, rd, rp, rxs);
    
    if it == maxitt
        fprintf( 'Execution terminated because reached iteration limit\n');
        lambda.lower = lambda_s;
        lambda.eqlin = lambda_eq;
        exit_code = 2;
        break
    end
    
    if gap <= eps
        fprintf( 'Execution terminated because duality gap reduced under the threshold\n');
        lambda.lower = lambda_s;
        lambda.eqlin = lambda_eq;
        exit_code = 1;
        break
    end


    X = diag(x);
    S = diag(lambda_s);
    M=2*Q+diag(lambda_s.*(1./x));
    
    J = [M, A';A, zeros(k, k)];
    r = -[ gradL - diag(1./x)*(sigma*mu*ones(n,1)) + lambda_s ;A*x-b];
    
    tol = 1e-14;
    maxit = size(J);
    maxit = maxit(1);
    
    if solver=="minres"
        [d, tmp] = minres(J, r, tol, maxit);
    end
    if solver=="gmres"
        [d, tmp] = gmres(J, r, [], tol, maxit);
    end
    if solver=="ldl"
        d = J\r;
    end
    
    dx    = d(1:n);
    dl_eq = d(n+1:k+n);
    dl_s  = diag(1./x)*(sigma*mu*ones(n, 1) - S*dx) - lambda_s; 
    

    alphax = alpha;
    alphas = alpha;

    % We check that negative direction entries do not violate non-neg constarints
    ind = dx < 0;  
    if any( ind )
       alphax = min(1 ,min( - x( ind ) ./ dx( ind ) ));
    end

    ind = dl_s < 0; 
    if any( ind )
       alphas = min(1, min( - lambda_s( ind ) ./ dl_s( ind ) ));
    end

    xprec = x; 
    fprec = fval;
    % Updating current point
    x = x + (0.995*alphax)*dx;
    lambda_s = lambda_s + (0.995*alphas)*dl_s;
    lambda_eq = lambda_eq + dl_eq;
    
    mu = sigma*(x'*lambda_s)/n;
end
time = toc;
fprintf( 'Primal-Dual Interior Point method terminated in %d iterations, elapsed time is %1.3e\nfval = %1.3e and complementary gap = %1.3e\n', it, time, fval, gap);
end