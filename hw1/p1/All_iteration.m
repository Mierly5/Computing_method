% implement Jacobi iteration
% to solve linear equation set
% PB18111679 fanweneddie

clear, clc

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
% the error bound
epsilon = 10^-15;
% maximal number of loops
max_loop = 10000;

% call function to implement jacobi iteration
jacobi_solution(A,b,x_exact,epsilon,max_loop);
gauss_siedel_solution(A,b,x_exact,epsilon,max_loop);
% call function to implement SOR iteration
SOR_solution(A,b,0.9,x_exact,epsilon,max_loop);
SOR_solution(A,b,0.8,x_exact,epsilon,max_loop);
SOR_solution(A,b,0.6,x_exact,epsilon,max_loop);
SOR_solution(A,b,0.5,x_exact,epsilon,max_loop);
SOR_solution(A,b,1.1,x_exact,epsilon,max_loop);
SOR_solution(A,b,1.3,x_exact,epsilon,max_loop);
SOR_solution(A,b,1.5,x_exact,epsilon,max_loop);
SOR_solution(A,b,1.6,x_exact,epsilon,max_loop);
SOR_solution(A,b,1.7,x_exact,epsilon,max_loop);
SOR_solution(A,b,1.8,x_exact,epsilon,max_loop);
% implementing jacobi iteration
% @A: the coefficent matrix on LHS
% @b: the constant matrix on RHS
% @x_exact: the exact solution to this equation set
% @epsilon: the error bound
% @max_loop: the maximal number of loops
% print the info during each iteration
% and plot the error in each iteration
function jacobi_solution(A,b,x_exact,epsilon,max_loop)
    [A_row,A_col] = size(A);
    [b_row,~] = size(b);
    % check whether the dimension of A and b matches.
    if(A_row ~= b_row)
        fprintf('the dimension of A and b does not match.\n');
        return;
    end
    
    % -------------------------------------------------
    % X(k+1) = RX(k) + g
    % -------------------------------------------------
    
    % x_cur is the solution of the current step of iteration
    % x_cur is initialized as a vector full of 0
    x_cur = zeros(A_col,1);
    
    % x_next is the solution of the next step of iteration
    % x_cur is initialized as a vector full of 1
    x_next = ones(A_col,1);
    
    % the error between current solution 
    % and exact solution in each loop
    error = zeros(A_col,1);
    % an vector to store the infinite norm of errors in each iteration
    error_list = zeros(max_loop,1);
    
    % R is the factor matrix in the loop
    % when i != j, R(i,j) = -A(i,j)/A(i,i)
    R = zeros(size(A));
    for i = 1:A_row
        for j = 1:A_col
            if i ~= j
                R(i,j) = -A(i,j)/A(i,i);
            end
        end
    end
    
    % g is the matrix to be added in each iteration
    % gi = bi/aii
    g = b;
    for i = 1:b_row
        g(i,1) = b(i,1)/A(i,i);
    end
    
    % show the items for the output info
    % x_size means the number of variables
    [x_size,~] = size(x_exact);
    % the row of items(such as loop time,option,variables)
    fprintf('  loop    option    ');
    for i = 1:x_size
        fprintf('  x%d      ',i);
    end
    fprintf('\n');
    
    % loop time
    loop = 0;
    % main iteration
    while( norm(x_cur-x_next,inf) > epsilon && loop <= max_loop )
        loop = loop + 1;
        % update x_cur
        x_cur = x_next;
        % update x_next
        x_next = R*x_cur + g;
        % record the error
        error = x_exact - x_cur;
        error_list(loop,1) = norm(error,inf);
        % show the info during each loop
        print_info(loop,x_cur,error);
    end
    % plot the max errors in each iteration
    title('Bisection')
    semilogy(error_list);
    hold on;
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
       fprintf('%10f',x_cur(i,1)); 
    end
    fprintf('\n');
    
    % show the error of each element in this loop
    fprintf('          error  ');
    for i = 1:x_size
       fprintf('%10f',abs(error(i,1)) ); 
    end
    fprintf('\n');
end

% implementing jacobi iteration
% @A: the coefficent matrix on LHS
% @b: the constant matrix on RHS
% @x_exact: the exact solution to this equation set
% @epsilon: the error bound
% @max_loop: the maximal number of loops
% print the info during each iteration
function gauss_siedel_solution(A,b,x_exact,epsilon,max_loop)
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
    % -------------------------------------------------
    
    % x_cur is the solution of the current step of iteration
    % x_cur is initialized as a vector full of 0
    x_cur = zeros(A_col,1);
    
    % x_next is the solution of the next step of iteration
    % x_cur is initialized as a vector full of 1
    x_next = ones(A_col,1);
    
    
    % D is the diagnal matrix of A
    D = diag(diag(A));
    
    % L is the lower triangle matrix of A
    L = tril(A,-1);
    
    % U is the upper triangle matrix of A
    U = triu(A,1);
    
    inv_D_plus_L = inv(D+L);
    S = -inv_D_plus_L * U;
    f = inv_D_plus_L * b;
    % show the items for the output info
    % x_size means the number of variables
    [x_size,~] = size(x_exact);
    
    % the error between current solution 
    % and exact solution in each loop
    error = zeros(A_col,1);
    % an vector to store the infinite norm of errors in each iteration
    error_list = zeros(max_loop,1);
    
    % the row of items(such as loop time,option,variables)
    fprintf('  loop    option    ');
    for i = 1:x_size
       fprintf('  x%d      ',i);
    end
    fprintf('\n');
    
    % loop time
    loop = 0;
    % main iteration
    while( norm(x_cur-x_next,inf) > epsilon ...
                && loop < max_loop )
        loop = loop + 1;
        % update the current solution
        x_cur = x_next;
        % get the next solution
        x_next = S*x_cur + f;
        % record the error
        error = x_exact - x_cur;
        error_list(loop,1) = norm(error,inf);
        % show the info during each loop
        print_info(loop,x_cur,error);
    end
    % plot the max errors in each iteration
    semilogy(error_list);
    hold on;
end

% implementing SOR iteration
% @A: the coefficent matrix on LHS
% @b: the constant matrix on RHS
% @w: the relaxation factor
% @x_exact: the exact solution to this equation set
% @epsilon: the error bound
% @max_loop: the maximal number of loops
% print the info during each iteration
function SOR_solution(A,b,w,x_exact,epsilon,max_loop)
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
    % -------------------------------------------------
    
    % x_cur is the solution of the current step of iteration
    % x_cur is initialized as a vector full of 0
    x_cur = zeros(A_col,1);
    
    % x_next is the solution of the next step of iteration
    % x_cur is initialized as a vector full of 1
    x_next = ones(A_col,1);
    
    % D is the diagnal matrix of A
    D = diag(diag(A));
    % temp result to calculate Sw and f
    w_D_inv = w*inv(D);
    
    % L is the lower triangle matrix of A
    L = tril(A,-1);
    
    % U is the upper triangle matrix of A
    U = triu(A,1);
    
    % I is identity matrix
    I = eye(A_col);
    
    % temp result to calculate Sw and f
    temp = inv(I + w_D_inv * L);
    
    % Sw is the factor matrix in the loop
    Sw = temp *( (1-w) * I - w_D_inv * U );
    
    % f is the matrix to be added in each iteration
    f = temp * w_D_inv * b;
    
    % the error between current solution 
    % and exact solution in each loop
    error = zeros(A_col,1);
    % an vector to store the infinite norm of errors in each iteration
    error_list = zeros(max_loop,1);
    
    % show the items for the output info
    % x_size means the number of variables
    [x_size,~] = size(x_exact);
    % the row of items(such as loop time,option,variables)
    fprintf('  loop    option    ');
    for i = 1:x_size
        fprintf('  x%d      ',i);
    end
    fprintf('\n');
    
    % loop time
    loop = 0;
    % main iteration
    while( norm(x_cur-x_next,inf) > epsilon ...
            && loop < max_loop)
        loop = loop + 1;
        %semilogy(loop,norm(x_cur-x_next,inf));
        x_cur = x_next;
        x_next = Sw*x_cur + f;
        % record the error
        error = x_exact - x_cur;
        error_list(loop,1) = norm(error,inf);
        % show the info during each loop
        print_info(loop,x_cur,error);
    end
    % plot the max errors in each iteration
    semilogy(error_list);
    hold on;
end
