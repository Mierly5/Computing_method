% We use spline to approximate a function in [-1,1]
% Here, funtion's derivative on edge is known,
% which satisfies boundary condition 2
% 
% PB18111679 fanweneddie (from USTC)

clear, clc

% The function to be solved
syms x;
F = @(x)sin(4 .* x .^ 2) + (sin(4 .* x)) .^ 2;

% at lease, we use 2^4 testing points
% at most we use 2^10 testing points
exp_down = 4;
exp_up = 10;
% stores the log number of points to be inserted
% here we only use 2^4, 2^5, ... ,2^10 testing points 
n_list = zeros(1,exp_up - exp_down + 1);
for i = 1 : exp_up - exp_down + 1
    n_list(i) = 2^(i + exp_down - 1);
end


% number of testing points to get the error on each point
k = 2000;

%option: 0 for (2) and 1 for (3)
opt = 1;

% call our spline function
spline_boundary_condition_2(F,n_list,k,opt);

% implementing spline with boundary condition 2
% @F: the function to approximate in [-1,1]
% @n_list: stores the number of intervals in [-1,1]
% @k: number of testing points in [-1,1]
% @opt: 0: plot the semilogy of errors on each point in the end
%       1: plot the semilogy of max error in each case 
function spline_boundary_condition_2(F,n_list,k,opt)
    
    % max_errors stores the max error 
    % of testing points on each n
    [~,n_list_col] = size(n_list);
    max_errors = zeros(1,n_list_col);
    
    % get the max_error on each n
    for index = 1 : n_list_col
    
        n = n_list(index);
        % get the info of intervals
        % x marks interval end points
        % h marks interval length
        % df marks the diff of the function value 
        % at the start and the end of the interval
        x = linspace(-1, 1, n + 1);
        f = F(x);
        h = diff(x);
        df = diff(f);
    
        % lambda and mu are parameters in A
        lambda = h(2:n) ./ ( h(2:n) + h(1:n - 1) );
        mu = 1 - lambda;
    
        % d is on LHS of the linear equation
        % get d(1) ~ d(n-1)
        d = 6 * ( df(2:n) ./ h(2:n) - df(1:n - 1) ./ h(1:n-1) )...
            ./ ( h(1:n-1) + h(2:n) );
    
        % for boundary condition 2,
        % get d(0) and d(n)
        m_0 = 0;
        m_n = 0;
    
        d = d.';
        d_0 = 6 * ( df(1) / h(1) - m_0 ) / h(1);
        d_n = 6 * (m_n - df(n) / h(n) ) / h(n);
        d = [d_0; d; d_n];
    
        % A is the parameters of this linear equation
        % get A
        A = 2 * eye(n+1);
        A(1,2) = 1;
        A(n+1,n) = 1;
        for i = 2 : n - 1
            A(i,i-1) = mu(i);
            A(i,i+1) = lambda(i);
        end
    
        % get M where AM = d
        M = A \ d;
        errors = get_errors(F,x,h,M,k);
        % for opt = 0, now k = 2^4
        % plot the errors of approximation on each testing point
        if(index == 1 && opt == 0)
            test_x = linspace(-1,1,k);
            semilogy(test_x,errors);
        end
        max_errors(index) = max(abs(errors));
    end
    
    % for opt = 1,
    % plot the max error of each interpolation
    if(opt == 1)
        loglog(n_list,max_errors);
        grid on;
    end
end

% get the error on each testing point
% @F: the function in [-1,1]
% @x: the points of each interval in [-1,1]
% @h: the length of each interval
% @M: the second order derivative
%   on the point of each interval
% @k: number of testing points
function errors = get_errors(F,x,h,M,k)
    test_x = linspace(-1,1,k);
    errors = test_x;
    i = 1;
    for j = 1 : k
        % get the interval that the testing point is in
        x_0 = test_x(j);
        while( x_0 > x(i+1) )
            i = i + 1;
        end
        gap_up = x(i+1) - x_0;
        gap_down = x_0 - x(i);
        % S is the approximate value
        S = ( gap_up^3 * M(i) + gap_down^3 * M(i+1) ) / (6*h(i)) + ...
            ( gap_up * F(x(i)) + gap_down * F(x(i+1)) ) / h(i) - ...
            h(i) * (gap_up * M(i) + gap_down * M(i+1)) / 6;
        errors(j) = abs( F(x_0) - S );
    end
end


