% iterative approach to solving a system Ax = b based on initial guess x_initial

% Jacobi uses update matrix M = D, where D is the diagonal of A

% use maxiter for the max number of iterations, and tol for the error bound
% i.e. return x such that ||Ax - b|| <= tol * ||b||
function [x, iter] = Jacobi(A, b, x_initial, maxiter, tol)

    xold = x_initial;
    [m, n] = size(x_initial);
    
    Dinv = eye(n)./diag(A);
    
    for k = 1:maxiter
        iter = k;
        
        % solving the update system is easier done by directly computing the
        % inverse of D, which is I./D
        xnew = xold + Dinv * (b - A * xold);
        x = xnew;
        xold = xnew;
        r = A * xnew - b;
        if (norm(r) <= tol * norm(b))
            break
        end
    end

end