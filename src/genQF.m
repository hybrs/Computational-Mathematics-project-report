function [Q, q] = genQF(n, density, eigs)
%function [Q, q] = genQF(m, density, eigs)
%
%  Generate a positive-semidefinite matrix Q \in R^{n \times n} and a vector q \in R^n
%  that will express the generic quadratic fuction 
%
%    			f(x) = x'Qx + q'x 
%
% Input:
%
% - n (integer): dimension of Q
%
% - density (real \in (0, 1]): density*n*n will be approximately the number of non-zero el in Q 
%
% - eigs (n \times 1 vector of real): eigenvectors of Q
%
% Output: the solution, with the following fields:
%
% - Q (n \times n  real matrix): positive-semidefinite matrix
%
% - q (n \times 1 real vector): the q vector in f(x)ing the gap for each iteration
%
Q = sprandsym(n,density,eigs);
a = -1;
b = 1;
q = (b-a).*rand(n,1) + a;

end