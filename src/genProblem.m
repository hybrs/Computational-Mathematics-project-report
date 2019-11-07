function [problem] = genProblem(n,m)

[Q,q]=genQF(n);
[A,b]=genKDS(m, n);
if rank(A) == m
    problem=struct('Q',Q,'q',q,'A',A,'b',b);
end
end