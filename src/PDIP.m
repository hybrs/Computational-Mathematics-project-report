function [x, fval, lambda, exit_code] =  PDIP(problem, maxit, eps, epsd, epsp)
warning("off", "all");
Q = problem.Q;
q = problem.q;

n = size(Q);
n = n(1);

A = problem.A;
b = problem.b;

k = size(A);
k = k(1);

%init the starting point (x0, lambda_eq0, lambda_s0)
[x0, lambda] = feasible_sp(A, Q, q);
lambda_eq = lambda.eqlin;
lambda_s = lambda.lower;

alpha = 0.1;
mu = (x0'*lambda_s)/n;
sigma = 0.3;

fprintf( 'Primal-Dual Interior Point method\n');
fprintf( 'iter\t\tgap\t\tdualf\tprimalf\tsTx\n' );
fprintf('-------------------------------------\n');

x = x0;
it = 0;

while true
    it=it+1;
    Qx = Q * x;
    qx = q'*x;
    fval = x'*Qx + qx;
    gradL = 2*Qx + q +A'*lambda_eq - lambda_s;
    
    primal = fval;
    dual = - lambda_eq' * b - x'*Qx ;
    
    gap = ( primal - dual ) / max( abs( primal ) , 1 );
    
    
    rd = norm(gradL, 2);
    rp = norm(A*x - b, 2);
    rxs = norm(lambda_s'*x, 2);
    
    fprintf( '%4d\t\t%1.3e\t\t%1.3e\t%1.3e\t%1.3e\n' , it , gap, rd, rp, rxs);
    
    if it == maxit
        fprintf( 'iter\t\tgap\t\tdualf\tprimalf\tsTx\n\n' );
        fprintf( 'Execution terminated because reached iteration limit\n');
        lambda.lower = lambda_s;
        lambda.eqlin = lambda_eq;
        exit_code = 3;
        break
    end
    
    if abs(gap) <= eps
        fprintf( 'iter\t\tgap\t\tdualf\tprimalf\tsTx\n\n' );
        fprintf( 'Execution terminated because duality gap reduced under the threshold\n');
        lambda.lower = lambda_s;
        lambda.eqlin = lambda_eq;
        exit_code = 1;
        break
    end
    
    if rd <= epsd && rp <= epsp && rxs <= eps
        fprintf( 'iter\t\tgap\t\tdualf\tprimalf\tsTx\n\n' );
        fprintf( 'Execution terminated because stopping criteria all matched\n');
        lambda.lower = lambda_s;
        lambda.eqlin = lambda_eq;
        exit_code = 2;
        break
    end
    
    
    
    if it > 1 && norm(x-xprec, 2) <= epsd && abs(fval-fprec) <= epsd
        fprintf( 'iter\t\tgap\t\t\tdualf\t\tprimalf\t\tsTx\n\n' );
        fprintf( 'Execution terminated because saddle point\n');
        lambda.lower = lambda_s;
        lambda.eqlin = lambda_eq;
        exit_code = 4;
        break
    end
    
    X = diag(x);
    S = diag(lambda_s);
    M=2*Q+diag(lambda_s.*(1./x));
    
    J = [M, A';A, zeros(k, k)];
    r = -[ gradL - diag(1./x)*(sigma*mu*ones(n,1)) + lambda_s ;A*x-b];
    
    tol = 1e-15;
    maxit  = size(J);
    maxit = maxit(1);
    
    
    [d, tmp] = gmres(J, r, [], tol, maxit);
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
    fprec = fval;
    x = x + (0.995*alphax)*dx;
    lambda_s = lambda_s + (0.995*alphas)*dl_s;
    lambda_eq = lambda_eq + dl_eq;
    
    mu = sigma*(x'*lambda_s)/n;
 %   if lambda < 0
 %       fprintf("Execution terminated because lambda < 0\n");
 %       exit_code = -2;
 %       break
 %   end
end