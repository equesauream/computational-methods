% modified QR factorization method to give a QL factorization:
% A = QL where Q is orthogonal (i.e. Q^T = Q^-1), and L is lower triangular

% similar algorithm as QR, except the computed Hausholder transformations
% act on the first k - 1 entries instead of the last k - 1 entries
function [Q, L] = QL_Factor(A)

    [m, n] = size(A);
    % A is square, so m = n

    v = zeros(n,n);
    
    % get triangular part
    for k = n:-1:1
        x = A(1:k, k);
        vk = x;
        vk(k) = vk(k) + sign(x(1)) * norm(x);
        vk = vk/norm(vk);
        v(1:k, k) = vk;
        for j = 1:k
            A(1:k, j) = A(1:k, j) - 2 * vk * (vk' * A(1:k, j));
        end
    end

    L = A;
    Q = zeros(n, n);

    % y = Qx
    function y = Qx(x)
        for p = 1:n
            x(1:p) = x(1:p) - 2 * v(1:p, p) * (v(1:p, p)' * x(1:p));
        end
        y = x;
    end

    I = eye(n);

    % recover Q by computing Qe1, Qe2, ...
    for q = 1:n
        Q(1:n, q) = Qx(I(:, q));
    end
end

% get the sign of a real number x
function sgn = sign(x)
    if x > 0
        sgn = 1;
    else
        sgn = -1;
    end
end