function [sol] = match_solver(p, solver)

    maxit = 100; eps = 1e-14; epsd = 1e-11; epsp = 1e-13;

    n = p.n;

    if solver == "quadprog"
        %opts = optimoptions(@quadprog,'Display', 'iter-detailed','Algorithm', 'interior-point-convex');
        opts = optimoptions(@quadprog,'Algorithm', 'interior-point-convex');
        tic;
        [x, fval, exit_code, out, lambda] = quadprog(2*p.Q, p.q, [], [], p.A, p.b, zeros(n,1), [],[], opts);
        its = out.iterations;
        time = toc;
        
        grad = norm(2*p.Q*x + p.q +p.A'*lambda.eqlin - lambda.lower, 2);
        primal = x'*p.Q*x + p.q'*x;
        dual = - lambda.eqlin' * p.b - x'*p.Q*x ;
        gap = ( primal - dual ) / max( abs( primal ) , 1 );
        rp = norm(p.A*x-p.b, 2);
        gaps = [gap];
        res = struct('rd', [norm(grad, 2)], 'rp', [rp]);
        
    else
        [x, fval, lambda, exit_code, time, its, gaps, res] =  PDIP(p, maxit, eps, epsd, epsp, solver);
    end
    
    
    
    if solver == "ldl"
        solver = 1;
    else if solver == "gmres"
            solver = 2;
        else if solver == "minres"
                solver = 3;
            else 
            solver = 4;
            end
        end
   end    
 
    
    sol = struct('x',x, 'fval', fval, 'exit_code', exit_code, 'lambda',lambda, 'solver', solver, 'time', time, 'n', p.n, 'm', p.m, 'density',p.density, 'gap', gaps, 'res', res, 'iterations', its);
end
