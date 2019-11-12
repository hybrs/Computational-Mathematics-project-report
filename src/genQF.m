function [Q, q] = genQF(n, density, eigs)
Q = sprandsym(n,density,eigs);
q = rand(n, 1);
end