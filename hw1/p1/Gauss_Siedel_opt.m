% implement Jacobi iteration
% to solve linear equation set
% PB18111679 fanweneddie

%clear, clc

% A is the coefficent matrix on LHS
A = [ 2,-1, 0, 0, 0, 0, 0, 0, 0, 0;
     -1, 2,-1, 0, 0, 0, 0, 0, 0, 0;
      0,-1, 2,-1, 0, 0, 0, 0, 0, 0;
      0, 0,-1, 2,-1, 0, 0, 0, 0, 0;
      0, 0, 0,-1, 2,-1, 0, 0, 0, 0;
      0, 0, 0, 0,-1, 2,-1, 0, 0, 0;
      0, 0, 0, 0, 0,-1, 2,-1, 0, 0;
      0, 0, 0, 0, 0, 0,-1, 2,-1, 0;
      0, 0, 0, 0, 0, 0, 0,-1, 2,-1;
      0, 0, 0, 0, 0, 0, 0, 0,-1, 2;];
% b is the constant matrix on RHS
b = [ 2;-2; 2;-1; 0; 0; 1;-2; 2;-2];

% the exact solution to this equation set
x_exact = [1;0;1;0;0;0;0;-1;0;-1];
% the error bound in iteration
epsilon = 10^-15;
% maximal number of loops
max_loop = 10000;

% call function to implement gauss siedel iteration
gauss_siedel_opt_solution(A,b,x_exact,epsilon,max_loop,1);

sum_t = 0;
for i = 1:10
    tic;
    gauss_siedel_opt_solution(A,b,x_exact,epsilon,max_loop,0);
    t = toc;
    sum_t = sum_t + t;
end

fprintf('calling optimize Gauss-Siedel for 10 times, sum time = %6fs\n',sum_t);


% implementing optimized Gauss-Siedel iteration
% @A: the coefficent matrix on LHS
% @b: the constant matrix on RHS
% @x_exact: the exact solution to this equation set
% @epsilon: the error bound
% @max_loop: the max number of loops
% @print: whether to print the final result(1) or not(0)
% print the info during each iteration
function gauss_siedel_opt_solution(A,b,x_exact,epsilon,max_loop,print)
    [A_row,A_col] = size(A);
    [b_row,~] = size(b);
    % check whether the dimension of A and b matches.
    if(A_row ~= b_row)
        fprintf('the dimension of A and b does not match.\n');
        return;
    end
    
    % -------------------------------------------------
    % originally, it should be
    % X(k+1) = SX(k) + f
    % where S = -(D+L)_inv U and f = -(D+L)_inv b
    % however, I don't use matrix multiplication
    % since the cost of getting inversion is unacceptible 
    % -------------------------------------------------
    
    % x_cur is the solution of the current step of iteration
    % x_cur is initialized as a vector full of 0
    x_cur = zeros(A_col,1);
    
    % x_next is the solution of the next step of iteration
    % x_cur is initialized as a vector full of 1
    x_next = ones(A_col,1);
    
    
    % show the items for the output info
    % x_size means the number of variables
    [x_size,~] = size(x_exact);
    % the row of items(such as loop time,option,variables)
    %fprintf('  loop    option    ');
    %for i = 1:x_size
     %   fprintf('  x%d      ',i);
    %end
    %fprintf('\n');
    
    % loop time
    loop = 0;
    % main iteration
    while( norm(x_cur-x_next,inf) > epsilon ...
            && loop < max_loop )
        loop = loop + 1;
        % update the current solution
        x_cur = x_next;
        % get the next solution 
        for i = 1:x_size
            sum = 0;
            for j = max(i-1,1):min(i+1,x_size)
                sum = sum + A(i,j)*x_next(j,1);
            end
            x_next(i,1) = (b(i,1) - sum) / A(i,i) + x_next(i,1);
        end
        % show the info during each loop
        %print_info(loop,x_cur,x_exact);
    end
    if(print == 1)
        print_info(loop,x_cur,x_exact);
    end
end

% show the info of solution and error in each iteration
% @loop: loop time
% @x_cur: the current solution
% @x_exact: the exact solution(calculated manually beforehand)
function print_info(loop,x_cur,x_exact)
    % x_size means the number of variables
    [x_size,~] = size(x_cur);
    
    fprintf('--------------------------------------------------\n');
    % show the solution in this loop
    fprintf('%5d     value  ',loop);
    for i = 1:x_size
       fprintf('%20.15f',x_cur(i,1)); 
    end
    fprintf('\n');
    
    % show the error in this loop
    fprintf('          error  ');
    for i = 1:x_size
       fprintf('%10f',abs(x_exact(i,1) - x_cur(i,1))); 
    end
    fprintf('\n');
end
    
