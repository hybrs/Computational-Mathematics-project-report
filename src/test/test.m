ns = [10, 15, 20, 35, 50, 70, 100];%, 120, 140, 160, 200, 250, 400, 500, 1000];
ms = [0.1, 0.3, 0.5, 0.7, 0.8];

sols = [];

for idxn=1:length(ns)
    for idxm=1:length(ms)
        n = idxn; m = round(idxm*n); density = 0.4; p = genProblem(n, m, density);
        [sol] = solveKDSQP(p, ["ldl", "minres", "gmres", "quadprog"]);
        sols = [sols; sol];
    end
end
