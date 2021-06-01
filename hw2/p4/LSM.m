% Here I implement Least Square Method(LSM)
% to make approximation based on f(x) = x / (a + bx)
%
% PB18111679 fanweneddie (from USTC)

clear, clc

% x-value of four original input points 
X = [2.1, 2.5, 2.8, 3.2];
% y-value of four original input points
Y = [0.6087, 0.6849, 0.7368, 0.8111];

% preprocess the data to fit in linear LSM
% now y = x / (a + bx) <=> y_inv = a*x_inv + b
X_inv = 1 ./ X;
Y_inv = 1 ./ Y;

alpha = Least_square_method(X_inv,Y_inv);

a = alpha(2);
b = alpha(1);

% The approximate function
syms x;
F = @(x) x ./ (a + b .* x);

% get the 2-norm of error on each points
Y_appro = F(X);
errors = Y - Y_appro;
err = norm(errors,2);
% print a,b,err 
fprintf('Approximate f(x) = x / (%10.6f +%10.6f x )\n',a,b);
fprintf('The 2-norm of errors is %10.6f\n',err);

% plot the approximation function
x= 2 : 0.01 :4;
y = F(x);
scatter(X,Y,'k*');
hold on;
plot(x,y);

% implement Least_square_method
% to get a linear approximate function Y = aX + b
% @X: x-value of input points 
% @Y: y-value of input points
% return the coeffient vector alpha
function alpha = Least_square_method(X,Y)
    % coefficients of linear function
    alpha = zeros(2,1);
    % number of input points
    n = length(X);
    % A stores different order of xi
    A = ones(n,2);
    for i = 1 : n
        A(i,2) = X(i);
    end
    
    % A^T * A * alpha = A^T * Y^T
    % that is, L * alpha = R
    L = A.' * A;
    R = A.' * Y.';
    alpha = L \ R;
end
