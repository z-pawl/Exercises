% Tridiagonal matrix algorithm
%     [ ⋱ ⋱       ]
%     [ ⋱ ⋱ c     ]
% A = [   ⋱ b ⋱   ]
%     [     a ⋱ ⋱ ]
%     [       ⋱ ⋱ ]
% b,d - nx1 vectors
% a,c - nx1 vectors - program will ignore the first and the last element of the vector respectively
function sol = tdma(a,b,c,d)

    % Dimension of the matrix A
    n = length(b);

    % Forward sweep
    for i = 2:n
        w = a(i)/b(i-1);
        b(i) = b(i) - w*c(i-1);
        d(i) = d(i) - w*d(i-1);
    end

    % Preallocating the solution vector
    sol = zeros(size(b));

    % Backward substitution
    sol(n) = d(n)/b(n);
    for i = (n-1):(-1):1
        sol(i) = (d(i)-c(i)*sol(i+1))/b(i);
    end
end

