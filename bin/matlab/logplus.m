function [ z ] = logplus( x, y )
%[ z ] = logplus( x, y )
%   This function performs the addition of two logarithmic values of the
%   form log(exp(x) + exp(y)). It avoids problems of dynamic range when the
%   exponentiated values of x or y are very large or small.

if x > y
    z = x+log(1.+exp(y-x));
else
    z = y+log(1.+exp(x-y));
end

end

