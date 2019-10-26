function [x, fval, lambda, exit_code] =  solvePDIPCQP(problem, eps)

Q = problem.Q;
q = problem.q;

n = size(Q);
n = n(1);

A = problem.A;
b = problem.b;

k = size(A);
k = k(1);

x0 = problem.x0;
maxIt = problem.maxit;

lambda_eq = rand(k, 1);%1e-6 * ones( k , 1 );
lambda_s = rand(n, 1);%1e-6 * ones( n , 1 );
lambda = struct('leq', lambda_eq, 'ls', lambda_s);

xprec =x0;
x = x0;
it = 0;

alpha = 0.1;
mu = (x'*lambda_s)/n;
sigma = 0.5;

fprintf( 'Primal-Dual Interior Point method\n');
fprintf( 'iter\tgap\tdualf\tprimalf\tsTx\tmu\tdeltax\tK\n\n' );
K = 0.0;


while true
    it=it+1;
    Qx = Q * x;
    fval = x'*Qx + q'*x;
    primal = x'*Qx + q' * x;
    dual = - lambda_eq' * b - x'*Qx ;
    gap = ( primal - dual ) / max( abs( primal ) , 1 );
    condd = norm(2*Qx + q +A'*lambda_eq - lambda_s, 2);
    condp = norm(A*x - b, 2);
    condxs = abs(lambda_s'*x);
    
    %local minima <=> x^k+1-x^k < eps
    %if it > 1 && norm(x-xprec, 2) <= eps
    %    fprintf( 'Execution terminated because of local minimum\n');
    %    lambda = struct('leq', lambda_eq, 'ls', lambda_s);
    %    exit_code = 1;
    %    break
    %end
    if abs(gap) <= eps
        fprintf( '%4d\t%1.3e |\t%1.3e\t%1.3e\t%1.3e\t%1.3e |\t%1.3e\t%1.3e\n' , it , gap, condd, condp, condxs, mu, norm(x-xprec, 2),K );
        fprintf( 'iter\tgap\tdualf\tprimalf\tsTx\tmu\tdeltax\tK\n\n' );
        fprintf( 'Execution terminated because stopping criterion was matched\n');
        lambda = struct('leq', lambda_eq, 'ls', lambda_s);
        exit_code = 2;
        break
    end
    if condp  <= eps && condd <= eps && condxs <= eps
        fprintf( '%4d\t%1.3e |\t%1.3e\t%1.3e\t%1.3e\t%1.3e |\t%1.3e\t%1.3e\n' , it , gap, condd, condp, condxs, mu, norm(x-xprec, 2),K );
        fprintf( 'iter\tgap\tdualf\tprimalf\tsTx\tmu\tdeltax\tK\n\n' );
        fprintf( 'Execution terminated because prima/dual feasibility and sTx conditions were matched\n');
        lambda = struct('leq', lambda_eq, 'ls', lambda_s);
        exit_code = 3;
        break
    end
    if it >= maxIt
        fprintf( '%4d\t%1.3e |\t%1.3e\t%1.3e\t%1.3e\t%1.3e |\t%1.3e\t%1.3e\n' , it , gap, condd, condp, condxs, mu, norm(x-xprec, 2),K );
        fprintf( 'iter\tgap\tdualf\tprimalf\tsTx\tmu\tdeltax\tK\n\n' );
        fprintf( 'Execution terminated because iteration reached threshold limit\n');
        lambda = struct('leq', lambda_eq, 'ls', lambda_s);
        exit_code = -1;
        break
    end
    
    if mu > eps
        mu = (x'*lambda_s)/(n*10);
    end
    
    X = diag(x);
    S = diag(lambda_s);
    M=2*Q+diag(lambda_s.*(1./x));
    
    J = [M, A';A, zeros(k, k)];
    r = -[ 2*Q*x+ q + A'*lambda_eq - lambda_s - diag(1./x)*(sigma*mu*ones(n,1));A*x-b];
    
    K = cond(J,2);
    
    fprintf( '%4d\t%1.3e |\t%1.3e\t%1.3e\t%1.3e\t%1.3e |\t%1.3e\t%1.3e\n' , it , gap, condd, condp, condxs, mu, norm(x-xprec, 2),K );
    
    tol = 1e-8;
    maxit  = size(J);
    maxit = maxit(1);
    
    d = gmres(J, r, [], tol, maxit);
    %d = linsolve(J, r);
    
    dx    = d(1:n);
    dl_eq = d(n+1:k+n);
    dl_s  = diag(1./x)*(sigma*mu*ones(n, 1) - S*dx) - lambda_s; 
    
    
    %min(dl);
    alphax = alpha;
    alphas = alpha;
    
    ind = dx < 0;  % negative direction entries
    if any( ind )
       alphax = min(1 ,min( - x( ind ) ./ dx( ind ) ));
    end
    ind = dl_s < 0;  % negative direction entries
    if any( ind )
       alphas = min(1, min( - lambda_s( ind ) ./ dl_s( ind ) ));
    end
    %maxt = 0.995*maxt;
    xprec = x; 
    x = x + (0.995*alphax)*dx;
    lambda_s = lambda_s + (0.995*alphas)*dl_s;
    lambda_eq = lambda_eq + (0.995*alphas)*dl_eq;
 %   if lambda < 0
 %       fprintf("Execution terminated because lambda < 0\n");
 %       exit_code = -2;
 %       break
 %   end
end