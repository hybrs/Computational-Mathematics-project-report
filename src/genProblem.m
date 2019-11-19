function [problem] = genProblem(n,m, density)
%function [problem] = genProblem(m, density, eigs)
%
%  Generate an instance of a QCP problem with k-disjoint simplices as constraints.
%				
%					 (P) = { min x'Qx + q'x | Ax = b, x >= 0}
% 
% Input:
%
% - n (integer): number of variables
%
% - m (integer): number of constraints
%
% - density (real \in (0, 1]): density*n*n will be approximately the number of non-zero el in Q 
%
% Output: a struct representing the problem, with the following fields:
%
% - problem.Q (n \times n  real matrix): positive-semidefinite matrix
%
% - problem.q (n \times 1 real vector): the q vector
%
% - problem.A (m \times n matrix): the simplices
%
% - problem.b (m \times 1 vector): vector of ones
%
% - problem.density (scalar): density of Q
%


% We generate a vector of eigenvalues >= 0
a = -1e-6;
b = 1;
eigs = (b-a).*rand(n,1) + a;
eigs(eigs<0) = 0;
% Generate Q and q
[Q,q]=genQF(n, density, eigs);
% Generate A and b
[A] = generateAdisjointed(m,n);
b = ones(m, 1);

problem=struct('Q',Q,'q',q,'A',A,'b',b, 'density', density, 'n', n, 'm', m);

end