function [problem] = genProblem(n,m)

[Q,q]=genQF(n);
[A,b]=genKDS(m, n);
 x0 = rand(n,1);
 maxit = 100;

 problem=struct('Q',Q,'q',q,'A',A,'b',b,'x0',x0,'maxit',maxit);
 
end