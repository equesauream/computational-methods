% finds the Cholesky decomposition of a symmetric PSD matrix A
% i.e. finds G so that A = GG^T where G is lower triangular
function G = Cholesky(A)

    [n, m] = size(A);
    G = zeros(n, n);
    
    % find the Cholesky factor
    for i = 1:n
        G(i, i) = sqrt(A(i, i));
    
        % set v/sqrt(alpha)
        G((i + 1):n, i) = A((i + 1):n, i) / sqrt(A(i, i));
    
        % update G_22 by taking B <- B - vv^T/alpha
        op = ( A(i + 1:n, i) * A(i + 1:n, i).' ) / A(i, i);
        A(i + 1:n, i + 1:n) = A(i + 1:n, i + 1:n) - op;
    end
end