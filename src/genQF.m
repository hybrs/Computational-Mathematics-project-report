function [Q, q] = genQF(n, density, eigs)
Q = sprandsym(n,density,eigs);
a = -1;
b = 1;
q = (b-a).*rand(n,1) + a;
end