% iterative approach to solving a system Ax = b
function [x, iter] = Jacobi(A, b, x_initial, maxiter, tol)

    xold = x_initial;
    [m, n] = size(x_initial);
    
    D = zeros(m, m);
    for i = 1:m
        D(i, i) = A(i, i);
    end
    Dinv = D;
    for i = 1:m
        Dinv(i, i) = 1/D(i, i);
    end
    
    for k = 1:maxiter
        iter = k;
        
        xnew = xold + Dinv * (b - A * xold);
        x = xnew;
        xold = xnew;
        r = A * xnew - b;
        if (norm(r) <= tol * norm(b))
            break
        end
    end

end