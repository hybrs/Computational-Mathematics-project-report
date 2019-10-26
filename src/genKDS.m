function [A, b] = genKDS(k, n)

found = 0;
it = 0;
b = ones(k, 1);

while found == 0
    
    it=it+1;
    A = randi([0 1],k,n);
    if rank(A) == k
        found = 1;
        break;
    end
    
    if it > 100
        mat = zeros(k);
        break
    end

end


end