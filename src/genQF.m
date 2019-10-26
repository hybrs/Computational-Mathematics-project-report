function [Q, q] = genQF(n)

G = rand(n, n);
Q = G'*G;
q = rand(n, 1);

end