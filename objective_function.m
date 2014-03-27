% created the total objective function , did not come in any use

% f_l(x) = (1/L) g(x -y_l) 
% g(x) = c^2 [ ( |x|/c ) - log ( 1 + ( |x|/c ) ) ]
function [ cost ] = objective_function ( x  , Y_l , L, c)
cost = 0;
   for i = 1: L 
        cost = cost + function_prototype( x, Y_l(i), L, c);
    end
end

