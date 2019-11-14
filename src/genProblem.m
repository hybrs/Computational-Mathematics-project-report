function [problem] = genProblem(n,m, density)

a = -1e-6;
b = 1;
eigs = (b-a).*rand(n,1) + a;

eigs(eigs<0) = 0;

[Q,q]=genQF(n, density, eigs);
[A] = generateAdisjointed(m,n);
b = ones(m, 1);
if rank(A) == m
    problem=struct('Q',Q,'q',q,'A',A,'b',b, 'density', density, 'n', n, 'm', m);
else
    fprintf("A has rank %d < %d", rank(A), m);
end
end