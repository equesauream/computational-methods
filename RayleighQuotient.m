% based on an initial vector v0, find the "closest" eigenvector using an
% iterative approach
function [v, lambda, iter] = RayleighQuotient(A, v0, maxiter, tol)

    [m, n] = size(A);
    
    vold = v0/norm(v0);
    lambda = vold' * (A * vold);
    for i = 1:maxiter
        w = (A - lambda * eye(n)) \ vold;
        vnew = w/norm(w);
    
        lambda = vnew' * (A * vnew);
    
        vold = vnew;
        v = vnew;
        iter = i;
    
        if norm(A * vnew - lambda * vnew) < tol
            break
        end
    end
end