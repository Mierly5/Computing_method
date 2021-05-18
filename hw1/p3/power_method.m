% implement power method
% to get the maximal eigenvalue 
% and the corresponding eigenvector of a square matrix
% PB18111679 fanweneddie

% clear, clc

% the square matrix to be solved
% -- A = [4,-1,1;16,-2,-2;16,-3,-1];

A = [ -148, -105,  -83,  -67; ...
       488,  343,  269,  216; ...
      -382, -268, -210, -170; ...
        50,   38,   32,   29];
A = -A;

A = [ 222, 580, 584, 786; ...
      -82,-211,-208,-288; ...
       37,  98, 101, 132; ...
      -30, -82, -88,-109];

% the regulated max number of loops
max_loop = 100;
% the predefined error bound
epsilon = 10^-15;

% call the function to use power method
power_method_solution(A,max_loop,epsilon);

% implement power method
% to get maximal eigenvalue and eigenvector
% @A: the input square matrix
% @max_loop: the max number of loops
% @epsilon: the predefined error bound 
% print maximal eigenvalue and the corresponding eigenvector
function power_method_solution(A,max_loop,epsilon)

    % check whether A is a square matrix
    [A_row,A_col] = size(A);
    if(A_row ~= A_col)
        fprintf('the input matrix is not a square matrix.\n');
        return;
    end
    
    % print the items of info
    fprintf('   loop    eigenvalue     eigenvector\n');
    
    % the eigenvector of "last" iteration
    % last iteration
    q_last = ones(A_col,1);
    % last odd iteration
    q_odd_last = q_last;
    % last even iteration
    q_even_last = q_last;
    
    
    % the eigenvector of "current" iteration
    % current iteration
    q_cur = q_last;
    % current iteration after being normalized
    q_cur_norm = q_last;
    % current odd iteration
    q_odd_cur = q_last;
    % current even iteration
    q_even_cur = q_last;
    
    
    % the eigenvalue
    lambda = 1;
    % the final eigenvector
    q_1 = q_cur;
    q_2 = q_cur;
    % marks whether there is only one max eigenvalue
    one_max_eigen = 1;
    
    
    for loop = 1:max_loop
        % iterate on q_cur and normalize it
        % odd loop
        if(mod(loop,2) == 1)
            q_odd_cur = A*q_last;
            q_cur = q_odd_cur;
            q_last = q_even_last;
        % even loop
        else
            q_even_cur = A*q_last;
            q_cur = q_even_cur;
            q_last = q_odd_last;
        end
        
        % normalize q_cur
        % get the element with max module
        max_module = -1;
        lambda = 0;
        for i = 1:A_col
            if( abs(q_cur(i,1)) > max_module )
                max_module = abs(q_cur(i,1));
                lambda = q_cur(i,1);
            end
        end
        q_cur_norm = q_cur/lambda;
 
        % print the loop time, current eigenvalue and eigenvector
        fprintf('%6d     %6f    (',loop,lambda);
        for i = 1:A_col
            fprintf('%6f  ',q_cur_norm(i,1));
        end
        fprintf(')\n');
        
        % update current eigenvector in the loop
        % odd loop
        if(mod(loop,2) == 1)
            q_odd_cur = q_cur_norm;
        % even loop
        else
            q_even_cur = q_cur_norm;
        end
        
        % check whether X(k) converges
        % if so, there is only a positive eigenvalue with max_modulo 
        % return lambda as eigenvalue and q_cur as eigenvector
        if( norm(q_cur_norm - q_last,inf) < epsilon )
            q_1 = q_cur_norm;
            break;
        % check whether X(k) and X(k+1) converges to opposite vectors
        % if so, there is only a negative eigenvalue with max_modulo  
        % return -lambda as eigenvalue and q_cur as eigenvector
        elseif( norm(q_cur_norm + q_last,inf) < epsilon )
            q_1 = q_cur_norm;
            lambda = -lambda;
            break;
        % check whether X(2k) and X(2k+1) converges to different vectors
        % if so, there is a pair of opposite eigenvalues with max_modulo
        % get lambda, -lambda and the corresponding eigenvectors and return
        elseif(norm(q_odd_cur - q_odd_last,inf) < epsilon ...
                && norm(q_even_cur - q_even_last,inf) < epsilon)
            % q_cur is before normalization
            q_next = A * q_cur;
            % lambda and -lambda are the eigenvalues that we want
            lambda = sqrt( q_next(1,1) / q_last(1,1) );
            % get two normalized eigenvectors
            q_1 = q_next + lambda * q_cur;
            q_1 = q_1 / norm(q_1,inf);
            q_2 = q_next - lambda * q_cur;
            q_2 = q_2 / norm(q_2,inf);
            % show that there are two max eigenvalues
            one_max_eigen = 0;
            break;
        else
            % update eigenvector in last iteration
            % odd loop
            if(mod(loop,2) == 1)
                q_odd_last = q_cur_norm;
                q_last = q_cur_norm;
            % even loop
            else
                q_even_last = q_cur_norm;
                q_last = q_cur_norm;
            end 
        end
    end
    
    % print the last result
    % one max eigenvalue
    if(one_max_eigen == 1)
        fprintf('\nat last\neigenvalue = %20.15f\n',lambda);
        fprintf('eigenvector is (');
        for i = 1:A_col
            fprintf('%20.15f\n',q_1(i,1));
        end
        fprintf(')\n');
    % two max eigenvalue
    else
        % show eigenvalue_1 and eigenvector_1
        fprintf('\nat last\neigenvalue_1 = %20.15f\n',lambda);
        fprintf('eigenvector_1 is (');
        for i = 1:A_col
            fprintf('%20.15f\n',q_1(i,1));
        end
        fprintf(')\n');
        
        % show eigenvalue_2 and eigenvector_2
        fprintf('\neigenvalue_2 = %20.15f\n',-lambda);
        fprintf('eigenvector_2 is (');
        for i = 1:A_col
            fprintf('%20.15f\n',q_2(i,1));
        end
        fprintf(')\n');
    end
end