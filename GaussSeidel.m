% iterative approach to solving a system Ax = b
function [x, iter] = GaussSeidel(A, b, x_initial, maxiter, tol)

    xold = x_initial;
    [m, n] = size(x_initial);
    
    D = diag(A);
    L = -tril(A);
    
    for k = 1:maxiter
        iter = k;
        
        % this direct solve is relatively cheap because D - L is triangular
        xnew = xold + (D - L) \ (b - A * xold);
        x = xnew;
        xold = xnew;
        r = A * xnew - b;
        if (norm(r) <= tol * norm(b))
            break
        end
    end

end