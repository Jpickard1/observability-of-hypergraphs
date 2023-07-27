function [Ae] = kronexp(A, e)
%KRONEXP Kronecker exponentiation of x with itself e times
%
%   Ae = A kron A kron ... kron A
%        |  ---  e  terms  ---  |
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: July 27, 2023

Ae = A;
for i=1:e-1
    Ae = kron(Ae, A);
end

end

