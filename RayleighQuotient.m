% based on an initial vector v0, find the "closest" eigenvector using an
% iterative approach

% use maxiter for the max number of iterations, and tol for the error bound
% i.e. return [V, Lambda] such that ||A * V[i] - Lambda[i] * V[i]|| < tol
% for each i
function [v, lambda, iter] = RayleighQuotient(A, v0, maxiter, tol)

    [m, n] = size(A);
    
    vold = v0/norm(v0);
    lambda = vold' * (A * vold);
    for i = 1:maxiter
        % solve the shifted system by removing the current closest
        % eigenvalue
        w = (A - lambda * eye(n)) \ vold;
        vnew = w/norm(w);
    
        % compute the rayleigh quotient
        lambda = vnew' * (A * vnew);
    
        vold = vnew;
        v = vnew;
        iter = i;
    
        if norm(A * vnew - lambda * vnew) < tol
            break
        end
    end
end