% find the spectrum (set of eigenvalues) of A and their associated
% eigenvectors V using the QR iteration algorithm

% this algorithm computes the QR decomposition QR = A, then uses RQ as the
% next A. Eventually, this process will converge to a diagonal matrix, and
% since each intermediate matrix is similar to the last, the eigenvalues
% and eigenvectors are the same.

% use maxiter for the max number of iterations, and tol for the error bound
% i.e. return [V, Lambda] such that ||A * V[i] - Lambda[i] * V[i]|| < tol
% for each i
function [V, Lambda, iter] = QRIteration(A, maxiter, tol)
    [m, n] = size(A);
    
    Aprev = A;
    Vprev = eye(n);
    Lambda = zeros(n, 1);
    for i = 1:maxiter
        iter = i;
        
        % find the QR decomposition of A
        [Q, R] = qr(Aprev);
        Vprev = Vprev * Q;

        % reverse the decomposition
        Aprev = R * Q;
    
        V = Vprev;
        Lambda(k) = diag(Aprev);
    
        % check if eigenvalues and eigenvectors are within tolerance
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