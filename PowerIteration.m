% finds an eigenvector of A whose associated eigenvalue is the largest magnitude
% by using an iterative approach from a starting vector v0

% use maxiter for the max number of iterations, and tol for the error bound
% i.e. return [v, lambda] such that ||A * v - lambda * v|| < tol

% Power Iteration uses the Rayleigh Quotient to approximate the eigenvalue
function [v, lambda, iter] = PowerIteration(A, v0, maxiter, tol)

    vold = v0;
    for i = 1:maxiter
        vnew = A * vold / norm(A * vold);
        vold = vnew;
        lambda = vnew' * (A * vnew);
    
        v = vnew;
        iter = i;
        if norm(A * vnew - lambda * vnew) < tol
            break
        end
    end

end