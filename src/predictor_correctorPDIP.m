function [x, fval, lambda, exit_code] =  predictor_correctorPDIP(problem, maxit, eps1, esp2)
warning("off", "all");
Q = problem.Q;
q = problem.q;

n = size(Q);
n = n(1);

A = problem.A;
b = problem.b;

k = size(A);
k = k(1);
alpha = 0.5; 


%init the starting point (x0, lambda_eq0, lambda_s0)
x0 = ones(n,1);
lambda_eq = ones(k, 1);%1e-6 * ones( k , 1 );
lambda_s = ones(n, 1);%1e-6 * ones( n , 1 );

fprintf( 'Primal-Dual Interior Point method\n');
fprintf( 'iter\t\tgap\t\tdualf\tprimalf\tsTx\n' );
fprintf('-------------------------------------\n');

x = x0;
it = 0;

r = [1;norm(Q,2);norm(A,2);norm(q,2);norm(b,2)];
ro = max(r);

while true
    
    sigma = 0;
    
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
    rxs = min(x); rxs = min(lambda_s); rxs = min(x.*lambda_s);%lambda_s'*x;
    
    fprintf( '%4d\t\t%1.3e\t\t%1.3e\t%1.3e\t%1.3e\n' , it , gap, rd, rp, rxs);
    
    if it == 50
        fprintf( 'iter\t\tgap\t\t\tdualf\t\tprimalf\t\tsTx\n\n' );
        fprintf( 'Execution terminated because reached iteration limit\n');
        lambda = struct('leq', lambda_eq, 'ls', lambda_s);
        exit_code = 3;
        break
    end
    
    if abs(gap) <= eps
        fprintf( 'iter\t\tgap\t\t\tdualf\t\tprimalf\t\tsTx\n\n' );
        fprintf( 'Execution terminated because duality gap reduced under the threshold\n');
        lambda = struct('leq', lambda_eq, 'ls', lambda_s);
        exit_code = 1;
        break
    end
    
    if rd <= ro*eps1 && rp <= ro*eps2 && rxs <= eps2
        fprintf( 'iter\t\tgap\t\t\tdualf\t\tprimalf\t\tsTx\n\n' );
        fprintf( 'Execution terminated because stopping criteria all matched\n');
        lambda = struct('leq', lambda_eq, 'ls', lambda_s);
        exit_code = 2;
        break
    end
    
    if it > 1 && norm(x-xprec, 2) <= 1e-10 && abs(fval-fvalp) <= 1e-10
        fprintf( 'iter\t\tgap\t\t\tdualf\t\tprimalf\t\tsTx\n\n' );
        fprintf( 'Execution terminated because saddle point\n');
        lambda = struct('leq', lambda_eq, 'ls', lambda_s);
        exit_code = 4;
        break
    end
    
    X = diag(x);
    S = diag(lambda_s);
    M = 2*Q+diag(lambda_s.*(1./x));
    
    J = [M, A';A, zeros(k, k)];
    r_aff = -[ gradL ;A*x-b];
    
    tol = 1e-15;
    maxit  = size(J);
    maxit = maxit(1);
    
    
    [d_aff, tmp] = gmres(J, r_aff, [], tol, maxit);
    %d = linsolve(J, r);
    
    mu = (x'*lambda_s)/n;
    
    
    dx_aff    = d_aff(1:n);
    dl_eq_aff = d_aff(n+1:k+n);
    dl_s_aff  = diag(1./x)*(- S*dx_aff) - lambda_s; 
    
    
    %min(dl);
    alphax_aff = alpha;
    alphas_aff = alpha;
    
    ind = dx_aff < 0;  % negative direction entries
    if any( ind )
       alphax_aff = min(1 ,min( - x( ind ) ./ dx_aff( ind ) ));
    end
    ind = dl_s_aff < 0;  % negative direction entries
    if any( ind )
       alphas_aff = min(1, min( - lambda_s( ind ) ./ dl_s_aff( ind ) ));
    end
    %maxt = 0.995*maxt;
    
    mu_aff = ((x + alphax_aff*dx_aff)'*(lambda_s+alphas_aff*dl_s_aff))/n;
    sigma = (mu_aff/mu)^3;
    tau = 0.995;
    
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %       
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

    
    r = -[ 2*Q*x+ q + A'*lambda_eq - lambda_s - diag(1./x)*(sigma*mu*ones(n,1));A*x-b];
    
    
    [d, tmp]= gmres(J, r, [], tol, maxit);
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
    xprec = x;
    fvalp = fval;
    x = x + (tau*alphax)*dx;
    lambda_s = lambda_s + (tau*alphas)*dl_s;
    lambda_eq = lambda_eq + (tau)*dl_eq;
end