function [gap, rd, rp] = print_stats(sols)

sz = size(sols); sz = sz(2);
fprintf("**************************** PRINTING %d SOLUTIONS STATS***********************************************\n", sz);

if sz == 1 || sz > 2
    for idx=1:sz
        sol = sols(idx);
        fval = sol.fval;
        fprintf("\n[%s]\nQ in R %dx%d with density %.2d, A in {0,1} %dx%d\n", sol.solver, sol.n, sol.n, sol.density, sol.m, sol.n);
        fprintf("took: %4d seconds\n", sol.time);
        fprintf("fval = %4d\n", fval);
        fprintf("gap = %4d\n", sol.gap(length(sol.gap)));
        fprintf("|Ax-b| = %4d\n", norm(sol.res.rp(length(sol.res.rp)), 2));
        fprintf("|grad_L| = %4d\n", norm(sol.res.rp(length(sol.res.rp)), 2));
        if idx < sz
            fprintf("----------------------------------\n")
        end
    end
    fprintf("\n**************************** %d STATS PRINTED ***********************************************\n", sz); 
else
    sol = sols(1);
    sol1 = sols(2);
        
    x = sol.x;
    lambda = struct(); lambda1 = struct();
    lambda.eqlin = sol.lambda.eqlin;
    lambda.lower = sol.lambda.lower;
    fval = sol.fval;
    x1 = sol1.x;
    lambda1.eqlin = sol1.lambda.eqlin;
    lambda1.lower = sol1.lambda.lower;
    fval1 = sol1.fval;

    fprintf("\nQ in R %dx%d with density %.2d\nA in {0,1} %dx%d\n\n\t\t\t\t%d\t\t\t\t\t%d\n", sol.n, sol.n, sol.density, sol.m, sol.n, sol.solver, sol1.solver);
    fprintf("-------------------------------------------------------------------\n");
    fprintf("fval \t\t= %.4d\t\t\t%.4d\n", fval, fval1);
    fprintf("gap \t\t= %.4d\t\t\t%.4d\n", sol.gap(length(sol.gap)), sol1.gap(length(sol1.gap)));
    fprintf("|Ax-b| \t\t= %.4d\t\t\t%.4d\n", norm(sol.res.rp(length(sol.res.rp)), 2), norm(sol1.res.rp(length(sol1.res.rp)), 2));
    fprintf("|grad_L| \t= %.4d\t\t\t%.4d\n\nDeltas between final solutions:\n", norm(sol.res.rd(length(sol.res.rd)), 2), norm(sol1.res.rd(length(sol1.res.rd)), 2));
    fprintf("|delta_x| \t\t\t\t= %.4d\n", norm(x-x1, 2));
    fprintf("|delta_lambda.eqlin| \t= %.4d\n",norm(lambda.eqlin- lambda1.eqlin, 2));
    fprintf("|delta_lambda.lower| \t= %.4d\n",norm(lambda.lower- lambda1.lower, 2));
    fprintf("-------------------------------------------------------------------\n");
    fprintf("%d took: %.4d seconds\n", sol.solver, sol.time);
    fprintf("%d took: %.4d seconds\n", sol1.solver, sol1.time);fprintf("\n***************************************************************************\n"); 
end
end