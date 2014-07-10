function Y = iseven(N)
% ISEVEN - Returns 1 if the input is a positive even integer, 0 otherwise.
%
% usage:
%
%   Y = iseven(N)
%
%  N    Numeric array.
%  
%  Y    Logical array of same size as N, with value 1 (0) if the
%       corresponding element of N is (is not) an even positive integer.  
%
% Note that ISEVEN tests the value of N, not the data storage type (i.e., 
% it does not check for integer type).
%
% See also isodd.

% ---- Check for sufficient command line arguments.
error(nargchk(1, 1, nargin));

% ---- Verify that N is a numeric array.
if ~isnumeric(N)
    error('N must be a numeric array.');
end

% ---- Verify that N is even.
Nbase = N/2;
Y = (Nbase == round(Nbase)) & N>0;

return
