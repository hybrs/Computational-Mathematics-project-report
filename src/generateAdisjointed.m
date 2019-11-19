function [A] = generateAdisjointed(m,n)
%function [A] = generateAdisjointed(m,n)
%
%  Generate m-disjointed simplices (MDS) from n variables. All simplex have this form:
%
%								\sum_{i \in I^h} x_i = 1 
% Input:
%
% - n (integer): number of variables
%
% - m (integer): number of constraints
%
% Output: the matrix A representing MDS
%
% - A (m \times n matrix): the constraints matrix


%generate a matrix of all zeros
A = zeros(m,n);

%step 1- 
%random permutation of 2*m (this ensure that each column has at
%least one element = 1 and each row at least 2 elements
randRows1 = randperm(m);
randRows2 = randperm(m);

%step 2- 
%generate a vector of (n-2m) random rows
randRows3 = randi(m,[1,n-2*m]); 

%concat the vectors
randRows = [randRows1,randRows2,randRows3];

%shuffling the elements
randRows=randRows(randperm(length(randRows)));

%update the matrix
for c=1:n
    A(randRows(1,c),c)=1;
end