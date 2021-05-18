% implement Newton iteration
% to get the solution of equation f(x) = 0
% where f(x) = x^3 - 3*x^2 + 2
% PB18111679 fanweneddie

clear, clc

% The function to be solved
syms x;
f = x^3 - 3*x^2 + 2;

% the bound of three ranges
% left range is [-3,0]
left_l = -3;
right_l = 0;
% middle range is [0,2]
left_m = 0;
right_m = 2;
% right range is [2,4]
left_r = 2;
right_r = 4;

% max loop time
max_loop = 100;
% error bound in iteration
epsilon = 10^-6;

% call the function to use Newton iteration
% to calcuate the approximate solution in three ranges
newton_solution(left_l,right_l,f,max_loop,epsilon);
newton_solution(left_m,right_m,f,max_loop,epsilon);
newton_solution(left_r,right_r,f,max_loop,epsilon);

%   use Newton iteration to get the solution
% of f(x) = 0 where left <= x <= right
%   print the info of temporary solution 
% and error during each interation
% @left: left boundary of the solution
% @right: right boundary of the solution
% @f: the consecutive real function
% @max_loop: maximal number of loops
% @epsilon: the error bound in this range
function newton_solution(left,right,f,max_loop,epsilon)
    % x_last is the solution in the last iteration
    x_last = (left + right)/2;
    % x_cur is the solution in the current iteration
    x_cur = 1;
    % f's differential
    diff_f = diff(f);
    % error in this loop
    error = 0;
    % the vector to store errors in each iteration
    errors = zeros(max_loop,1);
    % the evaluated order in iteration
    order = 0;
   
    % print the info
    dash_str = repmat('-', 1, 50);
    fprintf('%s\n',dash_str);
    display(f);
    fprintf('looking for a root in [%f,%f]\n',left,right);
    fprintf('  loop        x         error       order\n');
    
    % main loop
    for loop = 1:max_loop
        % set x_cur = x_last - f(x_last)/f'(x_last)
        x_cur = x_last - ...
               subs(f,symvar(f),x_last)/subs(diff_f,symvar(diff_f),x_last);
        error = abs(x_cur - x_last);
        errors(loop,1) = error;
        % print the current solution, current log error 
        % and approximate order
        fprintf('%5d    %10f    %10f',loop,x_cur,error);
        if(loop >= 3)
            order = log(errors(loop,1)/errors(loop-1,1)) ...
                    / log(errors(loop-1,1)/errors(loop-2,1));
            fprintf('  %10f\n',order);
        else
            fprintf('\n');
        end
        if( abs(x_cur - x_last) < epsilon )
            break;
        else
            x_last = x_cur;
        end
    end
    fprintf('at last, the solution is x = %10f\n',x_next);
end