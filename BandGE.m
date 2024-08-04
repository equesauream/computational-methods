% perform Gaussian Elimination on a matrix A (size n x n) with 
% upper bandwidth p and lower bandwidth q

% e.g. a diagonal matrix has upper bandwidth 0 and lower bandwidth 0
% whereas a tridiagonal matrix has upper bandwidth 1 and lower bandwidth 0

% runtime: O(npq)
function x = BandGE(A, b, p, q)
    [n, m] = size(A);
    
    x = zeros(n, 1);
    
    % iterate over columns
    for i = 1:n
        % iterate over rows after the i-th row
        for j = (i + 1):min(i + p + 1, n)
            leading_val = A(j, i);
            if (leading_val ~= 0)
                % assume the entries A(j, [1:i-1]) have already been zeroed out
                A(j, i:n) = A(j, i:n) - ( leading_val / A(i, i) ) * A(i, i:n);
    
                b(j) = b(j) - leading_val * b(i) / A(i, i);
            end
        end
    end
    
    % now do the same, but on the top diagonal entries
    for i = n:-1:1
        for j = max(i - q, 1):(i - 1)
            b(j) = b(j) - A(j, i) * b(i) / A(i, i);
            A(j, i) = 0;
        end
    end
    
    for i = 1:n
        x(i) = b(i) / A(i, i);
    end

end