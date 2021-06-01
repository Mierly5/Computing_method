% Here I implement Newton-interpolation algorithm
% 
% PB18111679 fanweneddie (from USTC)

clear, clc

% The Runge function to be solved
syms x;
F = @(x) 1 ./ (1 + 25 .* x .^ 2);

% at lease, we use 2^2 testing points
% at most we use 2^7 testing points
exp_down = 2;
exp_up = 7;
% stores the number of points to be inserted
% here we only use 2^2, 2^3, ... ,2^7 testing points 
n_list = zeros(1,exp_up - exp_down + 1);
for i = 1 : exp_up - exp_down + 1
    n_list(i) = 2^(i + exp_down - 1);
end

% number of testing points to test error 
k = 2000;

Newton_interpolation(F, n_list, k,0);

% use Newton interpolation to approximate a function
% @F: the function to approximate in [-1,1]
% @n_list: stores the number of points to be inserted
% @k: number of testing points in [-1,1]
% @rand: whether to randomly choose the sequence of
%   interpolation points(1) or not(0).
% plot the semilogy of max error on each n in n_list
function Newton_interpolation(F, n_list, k,rand)

    % max_errors stores the max error 
    % of testing points on each n
    [~,n_list_col] = size(n_list);
    max_errors = zeros(1,n_list_col);
    
    % get the max_error on each n
    for index = 1 : n_list_col
        
        n = n_list(index);
        % get the sequence of interpolation points
        % { ( x(i),f(x(i)) ) }
        if(rand == 1)
            rng(22);
            x0 = randperm(n + 1);
            x0 = x0 - 1;
            x0 = x0 *(pi/n);
        else
            x0 = linspace(0, pi, n + 1);
        end
        x = cos(x0);
        f = F(x);
    
        % u is the list of testing points
        u = linspace(-1, 1, k);
        % fu is the exact function value
        % of each testing points
        fu = F(u);
        % Nu is the approximation function value
        % of each testing points
        Nu = u;
         
        % g is difference of each order
        g = f;
        for i = 2 : n + 1
            for j = n + 1 : -1: i
                g(j) = (g(j) - g(j-1)) ...
                        / (x(j) - x(j-i+1)); 
            end
        end
        
        % get the approximate value of each testing point
        for i = 1 : k
            % init the parameter t and interpolation value newton
            t = 1;
            newton = g(1);
            for j = 2 : n + 1
                t = t * ( u(i) - x(j-1) );
                newton = newton + t * g(j); 
            end
            Nu(i) = newton;
        end
        max_errors(index) = max( abs(fu - Nu) );
    end
    
    % show the max error for each case
    fprintf('max error for each case is\n');
    for i = 1 : n_list_col
        fprintf('%20.12f\n',max_errors(i));
    end    
end
