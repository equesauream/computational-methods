% find the spectrum (set of eigenvalues) of A and their associated
% eigenvectors V using the QR iteration algorithm
function [V, Lambda, iter] = QRIteration(A, maxiter, tol)
    [m, n] = size(A);
    
    Aprev = A;
    Vprev = eye(n);
    Lambda = zeros(n, 1);
    for i = 1:maxiter
        iter = i;
        
        [Q, R] = qr(Aprev);
        Vprev = Vprev * Q;
        Aprev = R * Q;
    
        V = Vprev;
        for k = 1:n
            Lambda(k) = Aprev(k, k);
        end
    
        exitflag = true;
        for k = 1:n
            if norm(A * Vprev(1:n, k) - Aprev(k, k) * Vprev(1:n, k)) >= tol
                exitflag = false;
            end
        end
    
        if exitflag == true
            break
        end
    end
end