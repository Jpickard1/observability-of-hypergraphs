function [Jp] = Jp4(HG, p, x)
%JP4 This function computes the Jp matrices using the LyapProduct method
% from Kronecker Tools.
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: March 31, 2023

% Extract Amat, n, k
A = HG.adjTensor;           % Adjacency tensor as multi-way array (i.e. not tensor toolbox)
A = tensor(A);              % A tensor as tensor object
Amat = tenmat(A,1);         % Unfold A
Amat = Amat(:,:);           % tensor -> matrix
n = size(Amat,1);
k = length(size(A));

v = KroneckerPower(x, n*k-(2*n-1));
for p=n:2
    d = (p-1)*k-(2*p-3);
    v = LyapProduct(Amat, v, d);
end
Jp = v;
end


