function [sols] = solveKDSQP(p, solvers)
sz = size(solvers);
sz = sz(2);

[sol] = match_solver(p, solvers(1));
sols = [sol];

for idx=2:sz
    solver = solvers(idx);
    [sol] = match_solver(p, solver);
    sols(idx) = sol;
end