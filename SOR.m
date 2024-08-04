% iterative approach to solving a system Ax = b with relaxation parameter
% omega (omega < 1 is under-relaxation, omega > 1 is over-relaxation)
function [x, iter] = SOR(omega, A, b, x_initial, maxiter, tol)

xold = x_initial;
D = diag(A);
L = -tril(A);

for k = 1:maxiter
    iter = k;
    
    xnew = xold + (D./omega - L) \ (b - A * xold);
    x = xnew;
    xold = xnew;
    r = A * xnew - b;
    if (norm(r) <= tol * norm(b))
        break
    end
end

end