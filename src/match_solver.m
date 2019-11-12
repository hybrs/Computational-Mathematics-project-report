function [sol] = match_solver(p, solver)

    maxit = 100; eps = 1e-14; epsd = 1e-11; epsp = 1e-13;

    n = p.n;

    if solver == 1
        solver = 'ldl';
        [x, fval, lambda, exit_code, time] =  PDIP(p, maxit, eps, epsd, epsp, solver);
    else
        if solver == 2
            solver = 'gmres';
            [x, fval, lambda, exit_code, time] =  PDIP(p, maxit, eps, epsd, epsp, solver);
        else
            if solver == 3
                solver = 'minres';
                [x, fval, lambda, exit_code, time] =  PDIP(p, maxit, eps, epsd, epsp, solver);
            else
                solver = 'quadprog';
                opts = optimoptions(@quadprog,'Algorithm','interior-point-convex');
                tic;
                [x, fval, exit_code, out, lambda] = quadprog(2*p.Q, p.q, [], [], p.A, p.b, zeros(n,1), [],[], opts);
                time = toc;
            end
        end
    end
    grad = 2*p.Q*x + p.q +p.A'*lambda.eqlin - lambda.lower;
    primal = x'*p.Q*x + p.q'*x;
    dual = - lambda.eqlin' * p.b - x'*p.Q*x ;
    gap = ( primal - dual ) / max( abs( primal ) , 1 );
    rp_vec =p.A*x-p.b;
    
    sol = struct('x',x, 'fval', fval, 'exit_code', exit_code, 'lambda',lambda, 'solver', solver, 'time', time, 'n', p.n, 'm', p.m, 'density',p.density, 'gap', gap, 'rd', grad, 'rp', rp_vec);
end
