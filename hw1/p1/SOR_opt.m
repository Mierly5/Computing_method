% implement optimized SOR iteration
% to solve linear equation set
% PB18111679 fanweneddie

% clear, clc

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

% w is relaxation factor
w = 0.8;
% the error bound in iteration
epsilon = 10^-15;
% maximal number of loops
max_loop = 10000;


% call function to implement SOR iteration
SOR_opt_solution(A,b,w,x_exact,epsilon,max_loop,1);

sum_t = 0;
for i = 1:10
    tic;
    SOR_opt_solution(A,b,w,x_exact,epsilon,max_loop,0);
    t = toc;
    sum_t = sum_t + t;
end

fprintf('calling optimized SOR with w = %5f for 10 times, sum time = %6fs\n'...
                    ,w,sum_t);


% implementing optimized SOR iteration
% @A: the coefficent matrix on LHS
% @b: the constant matrix on RHS
% @w: the relaxation factor
% @x_exact: the exact solution to this equation set
% @epsilon: the error bound
% @max_loop: the maximal number of loops
% @print: whether to print the result
function SOR_opt_solution(A,b,w,x_exact,epsilon,max_loop,print)
    [A_row,A_col] = size(A);
    [b_row,~] = size(b);
    % check whether the dimension of A and b matches.
    if(A_row ~= b_row)
        fprintf('the dimension of A and b does not match.\n');
        return;
    end
    % check whether w is out of valid bound
    if(w < 0)
        fprintf('relaxation factor w should be positive.\n');
        return;
    end
    
    % -------------------------------------------------
    % X(k+1) = Sw*X(k) + f
    % but we use loop rather than matrix multiplication
    % -------------------------------------------------
    
    % x_cur is the solution of the current step of iteration
    % x_cur is initialized as a vector full of 0
    x_cur = zeros(A_col,1);
    
    % x_next is the solution of the next step of iteration
    % x_cur is initialized as a vector full of 1
    x_next = ones(A_col,1);
    
    % the temporary vector in each loop
    x_temp = ones(A_col,1);
    
    % the error between current solution 
    % and exact solution in each loop
    error = zeros(A_col,1);
    % an vector to store the infinite norm of errors in each iteration
    error_list = zeros(max_loop,1);
    
    % show the items for the output info
    % x_size means the number of variables
    [x_size,~] = size(x_exact);
    if(print == 1)
        % the row of items(such as loop time,option,variables)
        fprintf('  loop    option    ');
        for i = 1:x_size
            fprintf('  x%d      ',i);
        end
        fprintf('\n');
    end
    
    % loop time
    loop = 0;
    % main iteration
    while( norm(x_cur-x_next,inf) > epsilon ...
            && loop < max_loop)
        loop = loop + 1;
        x_cur = x_next;
        % update x_next with optimization to sparse matrix A
        for i = 1:A_row
            sum_1 = 0;
            sum_2 = 0;
            if(i > 1)
                sum_1 = A(i,i-1)*x_next(i-1,1);
            end
            if(i < A_row)
                sum_2 = A(i,i+1)*x_cur(i+1,1);
            end
            x_temp(i,1) = (b(i,1) - sum_1 - sum_2)/A(i,i);
            x_next(i,1) = w*x_temp(i,1) + (1-w)*x_cur(i,1);
        end
        % record the error
        error = x_exact - x_cur;
        error_list(loop,1) = norm(error,inf);
    end
    if(print == 1)
        print_info(loop,x_cur,error);
    end
end

% show the info of solution and error in each iteration
% @loop: loop time
% @x_cur: the current solution
% @error: the error between current solution 
% and exact solution in each loop
function print_info(loop,x_cur,error)
    % x_size means the number of variables
    [x_size,~] = size(x_cur);
    
    fprintf('--------------------------------------------------\n');
    % show the solution in this loop
    fprintf('%5d     value  ',loop);
    for i = 1:x_size
       fprintf('%20.15f',x_cur(i,1)); 
    end
    fprintf('\n');
    
    % show the error of each element in this loop
    fprintf('          error  ');
    for i = 1:x_size
       fprintf('%10f',abs(error(i,1)) ); 
    end
    fprintf('\n');
end
    
