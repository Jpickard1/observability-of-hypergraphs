function [Jp] = Jp3(HG, p, S)
%JP This function computes the Jp vectors in the shared overleaf document
%   in a reasonable amount of time. At no point is a matrix larger than
%   Amat stored, and there are never more than k-2 kronecker products taken
%   consecutively.
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

% Recursive base case
if p == k - 2
    xx = S(:,1);
    for i=2:size(S, 2)
        xx = kron(xx, S(:,i));
    end
    Jp = Amat * xx;
    return
end

% Recursive case
b = (p-1)*k-(2*p-3);
Snew = cell(b, 1);
%% Calculate new sets (old)
for j=1:length(Snew)
    ss = sym('x_%d', [n, b]);
    offset = 0;
    for i=1:b
        if i ~= j
            ss(:,i) = S(:,i+offset);
        else
            xx = S(:,i);
            offset = k-2;
            for l=1:k-2
                xx = kron(xx, S(:,i+l));
            end
            ss(:,i) = Amat * xx;
        end
    end
    Snew{j} = ss;
end

%% Calculate new sets (new)
for j=1:length(Snew)
    ss = sym('x_%d', [n, b]);
    offset = 0;
    for i=1:b
        if i ~= j
            ss(:,i) = S(:,i+offset);
        else
            xx = S(:,i);
            offset = k-2;
            for l=1:k-2
                xx = kron(xx, S(:,i+l));
            end
            ss(:,i) = Amat * xx;
        end
    end
    Snew{j} = ss;
end

%% Recursively call for each set
for i=1:length(Snew)
    if i ~= 1
        Jp = Jp + Jp3(HG, p-1, Snew{i});
    else
        Jp = Jp3(HG, p-1, Snew{i});
    end
end

end


