%% Observability of Toy Hypergraphs
%
%   This file computes the minimal set of observable nodes for a set of
%   small hypergraphs
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: July 27, 2023

function toy_hypergraphs(tiArr)
%% Preamble

% Add software dependencies
addpath(genpath('/nfs/turbo/umms-indikar/Joshua/Hypergraph-Analysis-Toolbox/'));
addpath(genpath('/nfs/turbo/umms-indikar/Joshua/Kronecker Tools'));
addpath(genpath('/nfs/turbo/umms-indikar/Joshua/tensor_toolbox/'));

% Set of possible parameters
N=3:8;
colNames = cell(length(N), 1);
for i=1:length(N)
    colNames{i} = char(string(N(i)));
end
K=2:7;
type = ["hyperring", "hyperchain", "hyperstar"];

r = containers.Map;
for ti=1:length(type)
    t = type(ti);

    % rows are number of vertices, columns are order k
    r(t) = cell2table(cell(length(K), length(N)), 'VariableNames', colNames);
end

%% Execute Experiment
tic;
for ki=1:length(K)
    k = K(ki);
    for ni=1:length(N)
        n = N(ni);
        if n < k
            continue
        end
        t = type(tiArr);

        % Print parameters
        disp(string(t) + "(" + string(n) + "," + string(k) + ")");

        % Experiment
        HG = HAT.toyHG(n,k,t);               % Get hypergraph
        O = HGObsvSym(HG);
        [D, ~] = greedyMON(O, n);             % Greedy Node Selection

        T = r(t);
        T{ki, ni} = {D};
        r(t) = T;

        % Save result
        disp(D);

        fileName = "toyHG/" + string(t) + "_sym_2.mat";
        cmd = "save " + fileName + " " + "r -v7.3";
        eval(cmd);
        disp(cmd);
    end
end
disp(toc);

disp(r(t))

%% Postscript

end