% solves Ax = b by using Gaussian Elimination
% assume that A is invertible with size n x n, and that no row swaps are needed

% runtime: O(n^3)
function x = GaussElim(A, b)

    [n, m] = size(A);
    
    % iterate over columns
    for i = 1:n
        % iterate over rows after the i-th row
        for j = (i + 1):n
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
        for j = (i - 1):-1:1
            b(j) = b(j) - A(j, i) * b(i) / A(i, i);
            A(j, i) = 0;
        end
    end
    
    % now, A is a diagonal matrix
    x = zeros(n, 1);
    for i = 1:n
        x(i) = b(i) / A(i, i);
    end

end