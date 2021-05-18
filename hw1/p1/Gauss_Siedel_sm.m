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
gauss_siedel_solution(A,b,x_exact,epsilon,max_loop,1);

sum_t = 0;
for i = 1:10
    tic;
    gauss_siedel_solution(A,b,x_exact,epsilon,max_loop,0);
    t = toc;
    sum_t = sum_t + t;
end

fprintf('calling optimized Gauss-Siedel for 10 times, sum time = %6fs\n',sum_t);


% implementing jacobi iteration
% @A: the coefficent matrix on LHS
% @b: the constant matrix on RHS
% @x_exact: the exact solution to this equation set
% @epsilon: the error bound
% @max_loop: the maximal number of loops
% @print: whether to print the final result(1) or not(0)
% print the info during each iteration
function gauss_siedel_solution(A,b,x_exact,epsilon,max_loop,print)
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
%     fprintf('  loop    option    ');
%     for i = 1:x_size
%        fprintf('  x%d      ',i);
%     end
%     fprintf('\n');

    % loop time
    loop = 0;
    % main iteration
    while( norm(x_cur-x_next,inf) > epsilon ...
                && loop < max_loop )
        loop = loop + 1;
        % update the current solution
        x_cur = x_next;
        % get the next solution
        x_next = sparse_matrix_mult(S,x_cur) + f;
        % record the error
        error = x_exact - x_cur;
        error_list(loop,1) = norm(error,inf);
        % show the info during each loop
        %print_info(loop,x_cur,error);
    end
    if(print)
       print_info(loop,x_cur,error); 
    end
    % plot the max errors in each iteration
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
    


% function to implement sparse matrix multiplication
% @M1: one sparse matrix operand
% @M2: another sparse matrix operand
% @Mr: the result matrix of the multiplication
% T = O( M1_row*M1_col + M2_row*M2_col
%        + #effective multiplication + M1_row*M2_col)
%   << O(M1_row*M1_col*M2_col) before optimization
function [Mr,state] = sparse_matrix_mult(M1,M2)

    [M1_row,M1_col] = size(M1);
    [M2_row,M2_col] = size(M2);
    
    % check whether the dimension of two operands matches.
    if(M1_col ~= M2_row)
        fprintf('the dimension of two matrix operands does not match\n');
        result = 0;
        state = 1;
        return;
    end
    
    % init Mr
    Mr_row = M1_row;
    Mr_col = M2_col;
    Mr = zeros(Mr_row,Mr_col);
    
    % specific data structure for sparse matrix M1,M2
    % -------------------------------------------------
    % triple for a non-zero element in matrix
    % marks the row and col it is in and its value
    triple.row = 0;
    triple.col = 0;
    triple.value = 0;
    
    % number of non-zero elements in M1
    M1_non_zero_num = 0;
    % a sequential collection of all non-zero elements in M1
    M1_seq = [];
    % the position in M1_triples
    % of each row's first non-zero element in M1
    % init all elements in M1_seq_pos as M1_size + 1
    M1_seq_pos = zeros(M1_row,1);
    
    % same for M2
    M2_non_zero_num = 0;
    M2_seq = [];
    M2_seq_pos = zeros(M2_row,1);
    % --------------------------------------------------
    
    
    % get the triples their position in sequential list for M1
    % T1 = O( M1_row*M1_col )
    for i = 1:M1_row
        % marks that row i's elements are all zero
        all_zero_row = 1;
        % init M1_seq_pos(i,1)
        M1_seq_pos(i,1) = M1_non_zero_num;
        for j = 1:M1_col
            % record the non-zero elements in each row
            if(M1(i,j) ~= 0)
                M1_non_zero_num = M1_non_zero_num + 1;
                triple.row = i;
                triple.col = j;
                triple.value = M1(i,j);
                % append the triple of this element to M1_seq
                M1_seq = [M1_seq,triple];
                % record the first non-zero element in row i
                if(all_zero_row == 1)
                    M1_seq_pos(i,1) = M1_non_zero_num;
                    all_zero_row = 0;
                end
            end
        end
    end
    
    % get the triples their position in sequential list for M2
    % T2 = O( M2_row*M2_col )
    for i = 1:M2_row
        % marks that row i's elements are all zero
        all_zero_row = 1;
        % init M2_seq_pos(i,1)
        M2_seq_pos(i,1) = M2_non_zero_num;
        for j = 1:M2_col
            % record the non-zero elements in each row
            if(M2(i,j) ~= 0)
                M2_non_zero_num = M2_non_zero_num + 1;
                triple.row = i;
                triple.col = j;
                triple.value = M2(i,j);
                % append the triple of this element to M1_seq
                M2_seq = [M2_seq,triple];
                % record the first non-zero element in row i
                if(all_zero_row == 1)
                    M2_seq_pos(i,1) = M2_non_zero_num;
                    all_zero_row = 0;
                end
            end
        end
    end
    
    % iterate on M1's every row, to get Mr's every row
    % Amortized T3 = O( #effective multiplication + M1_row*M2_col) 
    for M1_r = 1:M1_row
        % accumulator for elements of each column in that row
        Mr_sum = zeros(Mr_col,1);
        % the position in M1_seq of the last 
        % non-zero element to be calculated
        M1_last_ele_seq_pos = 0;
        if(M1_r < M1_row)
            M1_last_ele_seq_pos = M1_seq_pos(M1_r + 1,1) - 1;
        else
            M1_last_ele_seq_pos = M1_non_zero_num;
        end
        
        % iterate on every non-zero elements in this M1's row
        for M1_pos = M1_seq_pos(M1_r,1):M1_last_ele_seq_pos
            % a row of elements in M2 will multiply with that element in M1
            % and that row corresponds to the element's column in M1
            M2_r = M1_seq(M1_pos).col;
            % the position in M2_seq of the last 
            % non-zero element to be calculated
            M2_last_ele_seq_pos = 0;
            if(M2_r < M2_row)
                M2_last_ele_seq_pos = M2_seq_pos(M2_r + 1,1) - 1;
            else
                M2_last_ele_seq_pos = M2_non_zero_num;
            end
            
            % iterate on every non-zero elements in this M2's row
            % and update the accumulator Mr_sum
            for M2_pos = M2_seq_pos(M2_r,1):M2_last_ele_seq_pos
                Mr_c = M2_seq(M2_pos).col;
                Mr_sum(Mr_c,1) = Mr_sum(Mr_c,1) + ...
                        M1_seq(M1_pos).value * M2_seq(M2_pos).value;
            end
        end
        
        % get the row in Mr from the final result of accumulator
        for Mr_c = 1: Mr_col
            Mr(M1_r,Mr_c) = Mr_sum(Mr_c,1);
        end
    end
        
    state = 0;
end    
