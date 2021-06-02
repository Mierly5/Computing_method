% Here I implement Lagrange-interpolation algorithm
% to approximate a periodic function
%
% PB18111679 fanweneddie (from USTC)

clear, clc

% The periodic function to be solved
syms x;
F = @(x) sin(2 * pi .* x) .* exp( cos(2 * pi .* x));
% the base function when n is odd
L_k_odd = @(x,k,n) (-1)^k / n ...
            * sin(n*pi .* x) * csc( pi*(x - k/n) );
% the base function when n is even
L_k_even = @(x,k,n) (-1)^k / n ...
            * sin(n*pi .* x) * cot( pi*(x - k/n) );


% number of interpolation points
n = 2^6;

% number of testing points to test error 
k = 1000;

Lagrange_interpolation(F, L_k_odd, L_k_even, n, k);

% use Lagrange interpolation to approximate a function
% @F: the function to approximate in [-1,1]
% @L_k_odd: the base function when n is odd
% @L_k_even: the base function when n is even
% @n: number of interpolation points
% @k: number of testing points in [0,1]
% plot the semilogy of errors
function Lagrange_interpolation(F,L_k_odd, L_k_even, n, k)

    % get the sequence of interpolation points
    x = linspace(0, 1, n);
    
    % select base function
    if( mod(n,2) == 1)
        l = L_k_odd;
    else
        l = L_k_even;
    end
    
    % get the testing points
    test_x = linspace(0,1,k);
    % the actual function value of testing points
    f = F(test_x);
    % the evaluated function value of testing points
    L = zeros(1,k);
    
    % get the error of each testing point
    for i = 1 : k
        tmp = 0;
        % sum on each interpolation point
        for j = 1 : n
            tmp = tmp + l(test_x(i), j, n) * f(i);
        end
        L(i) = tmp;
    end
    % get the errors on each testing point and plot them
    error = abs(L-f);
    semilogy(test_x,error);
    grid on;
end
