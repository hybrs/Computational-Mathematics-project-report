function [sol, qpsol] = compareWquadprog(p)

maxit  = 100; eps = 1e-14, epsd = 1e-11; epsp = 1e-13;

n = size(p.q);
n = n(1);

tic;
[x, fval, lambda, exit_code] =  PDIP(p, maxit, eps, epsd, epsp);
fprintf("PDIP took:\n");
toc;
fprintf("\n_______________________________________________\n")
sol = struct('x',x, 'fval', fval, 'exit_code', exit_code, 'lambda',lambda);

tic;
[x1, fval1, exit_code1, out, lambda1] = quadprog(2*p.Q, p.q, [], [], p.A, p.b, zeros(n,1), []);
fprintf("quadprog took:\n")
toc;
fprintf("_______________________________________________\n\n")
qpsol = struct('x',x1, 'fval', fval1, 'exit_code', exit_code1, 'out', out, 'lambda',lambda1);

grad = 2*p.Q*x + p.q +p.A'*lambda.eqlin - lambda.lower;
qpgrad =  2*p.Q*x1 + p.q +p.A'*lambda1.eqlin - lambda1.lower;

primal = x'*p.Q*x + p.q'*x;
primal_qp = x1'*p.Q*x1 + p.q'*x1;
dual = - lambda.eqlin' * p.b - x'*p.Q*x ;
dual_qp =  - lambda1.eqlin' * p.b - x1'*p.Q*x1 ;
gap = ( primal - dual ) / max( abs( primal ) , 1 );
qpgap = ( primal_qp - dual_qp ) / max( abs( primal_qp ) , 1 );

fprintf("fval = '%4d\t\t\tfval_qp = '%4d\n", fval, fval1);
fprintf("gap = '%4d\t\t\tgap_qp = '%4d\n", gap, qpgap);
fprintf("|Ax-b| = '%4d\t\t\t|Ax_qp-b| = '%4d\n", norm(p.A*x-p.b, 2), norm(p.A*x1-p.b, 2));
fprintf("|grad_L| = '%4d\t\t\t|grad_L_qp| = '%4d\n\n-------------------------------------------------------------------\n", norm(grad, 2), norm(qpgrad, 2));
fprintf("|delta_x| = '%4d\n", norm(x-x1, 2));
fprintf("|delta_lambda.eqlin| = '%4d\n",norm(lambda.eqlin- lambda1.eqlin, 2));
fprintf("|delta_lambda.lower| = '%4d\n",norm(lambda.lower- lambda1.lower, 2));

end