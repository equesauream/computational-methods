% iterative approach to solving a system Ax = b based on initial guess
% x_initial and relaxation parameter omega (omega < 1 is under-relaxation, omega > 1 is over-relaxation)

% SOR uses update matrix M = D/omega - L, where D is the diagonal of A, and L
% is the lower triangle of A

% use maxiter for the max number of iterations, and tol for the error bound
% i.e. return x such that ||Ax - b|| <= tol * ||b||
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