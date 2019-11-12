function [A] = generateAdisjointed(m,n)

%generate a matrix of all zeros
A = zeros(m,n);

%step 1- 
%random permutation of m (this ensure that each column has at
%least one element = 1

randRows1 = randperm(m);

%step 2- 
%generate a vector of (n-m) random rows
randRows2 = randi(m,[1,n-m]); 

%concat the two vectors
randRows = [randRows1,randRows2];

%shuffling the elements
randRows=randRows(randperm(length(randRows)));

%update the matrix
for c=1:n
    A(randRows(1,c),c)=1;
end