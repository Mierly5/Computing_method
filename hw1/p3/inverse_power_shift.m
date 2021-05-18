% implement inverse-power and shift 
% to get the eigenvalue closest to 0.8-0.6i
% and its corresponding eigenvector in a stochastic matrix
% PB18111679 fanweneddie

clear,clc

rng(2);
A = rand(100,100);

% the regulated max number of loops
max_loop = 10000;
% the predefined error bound
epsilon = 10^-14;
% the predefined eigenvalue
pre_value = 0.8-0.6i;

% call the function to get the eigenvalue and eigenvector that we want
inverse_power_shift_solution(A,max_loop,epsilon,pre_value);


% implement inverse-power and shift method
% to get eigenvalue closest to parameter p
% also we can get its corresponding eigenvector
% @A: the input square matrix
% @max_loop: the max number of loops
% @epsilon: the predefined error bound 
% @pre_value: the predefined eigenvalue
% print maximal eigenvalue and the corresponding eigenvector
function inverse_power_shift_solution(A,max_loop,epsilon,pre_value)

    % check whether A is a square matrix
    [A_row,A_col] = size(A);
    if(A_row ~= A_col)
        fprintf('the input matrix is not a square matrix.\n');
        return;
    end
    
    % print the items of info
    fprintf('   loop    eigenvalue\n');
    
    % last iteration
    q_last = zeros(A_col,1);
    q_last(1,1) = 1;
    
    % current iteration
    q_cur = q_last;
    
    % the eigenvalue
    numbda = 1;
    
    % the factor to multiply in each loop
    M = A - pre_value * eye(A_row,A_col);
    
    for loop = 1:max_loop
        % iterate on q_cur and normalize it
        % (A - pre_value*I)*X_(k+1) = Y(k)
        % that is, M*q_cur = q_last; 
        q_cur = LU(M,q_last);
        
        % normalize q_cur
        % get the element with max module
        max_module = -1;
        numbda = 0;
        for i = 1:A_col
            if( abs(q_cur(i,1)) > max_module )
                max_module = abs(q_cur(i,1));
                numbda = q_cur(i,1);
            end
        end
        q_cur = q_cur/numbda;
       
        % print the loop time, current eigenvalue and eigenvector
        fprintf('%6d     %6f + %6fi\n',loop,real(numbda),imag(numbda));
        
        % check whether X(k) converges
        % if so, there is only a positive eigenvalue with min_modulo 
        % return numbda as eigenvalue and q_cur as eigenvector
        if( norm(q_cur - q_last,inf) < epsilon )
            break;
        % check whether X(k) and X(k+1) converges to opposite vectors
        % if so, there is only a negative eigenvalue with min_modulo  
        % return -numbda as eigenvalue and q_cur as eigenvector
        elseif( norm(q_cur + q_last,inf) < epsilon )
            numbda = -numbda;
            break;
        else
            % update eigenvector in last iteration
            q_last = q_cur;
        end 
    end
    
    % print the last result
    % the eigenvalue closest to pre_value
    result = pre_value + 1/numbda;
    fprintf('\nat last\n loop =%4d,eigenvalue = %20.13f + %20.13fi\n' ...
            ,loop,real(result),imag(result));
    % the corresponding eigenvector
    fprintf('eigenvector is [\n');
    for i = 1:A_col
        fprintf(' %20.13f + %20.13fi\n ', ...
            real(q_cur(i,1)),imag(q_cur(i,1)));
    end
    fprintf(']\n');
end

% solve the equation set AX = b by LU decomposition
% @A: the factor square matrix on LHS
% @b: the result matrix on RHS
% return the solution matrix X
function X = LU(A,b)
    
    [A_row,A_col] = size(A);
    [b_row,~] = size(b);
    
    % check whether the dimension of A and b is valid
    if(A_row ~= A_col)
        fprintf('input factor matrix A should be a square matrix.\n');
        return;
    elseif(A_row ~= b_row)
        fprintf('the dimension of A and b does not match.\n');
        return;
    end
    
    % the lower triangle matrix whose diagnal elements are 1
    L = eye(A_row,A_col);
    % the upper triangle matrix whose diagnal elements are 0
    U = zeros(A_row,A_col);
    % the temp result when solving this equation
    Y = zeros(A_row,1);
    % the final solution of the equation.
    X = zeros(A_row,1);
    % temporary sum to be used in iteration
    sum = 0;
    
    % Doolittle decomposition
    % set A = LU
    % then AX = b <=> LUX = b.
    % set Y = UX,
    % then after solving LY = b, we can get X quickly
    
    % 1. get L and U
    % iterate with row number
    for k = 1:A_row
        % calculate the elements in row k in U
        for j = k:A_col
            sum = 0;
            for r = 1:k-1
                sum = sum + L(k,r)*U(r,j);
            end
            U(k,j) = A(k,j) - sum;
        end
        % calculate the elements in colomn k in L
        for i = k+1:A_row
            sum = 0;
            for r = 1:k-1
                sum = sum + L(i,r)*U(r,k);
            end
            L(i,k) = (A(i,k) - sum)/U(k,k);
        end
    end
    
    % 2. solve LY = b
    for i = 1:A_row
        sum = 0;
        for j = 1:i-1
            sum = sum + L(i,j)*Y(j,1);
        end
        Y(i,1) = b(i,1) - sum;
    end
    
    % 3. solve UX = Y
    for i = A_row:-1:1
        sum = 0;
        for j = i+1:A_row
            sum = sum + U(i,j)*X(j,1);
        end
        X(i,1) = ( Y(i,1) - sum ) / U(i,i);
    end
end

