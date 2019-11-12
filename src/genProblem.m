function [problem] = genProblem(n,m, density)

a = 0;
b = 1;
r = (b-a).*rand(n,1) + a;

r(r<0) = 0;

eigs = r;

[Q,q]=genQF(n, density, eigs);
[A] = generateAdisjointed(m,n);
b = ones(m, 1);
if rank(A) == m
    problem=struct('Q',Q,'q',q,'A',A,'b',b, 'density', density, 'n', n, 'm', m);
else
    fprintf("A has rank %d < %d", rank(A), m);
end
end