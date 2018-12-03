function [ J_opt, u_opt_ind ] = ValueIteration( P, G )
%VALUEITERATION Value iteration
%   Solve a stochastic shortest path problem by Value Iteration.
%
%   [J_opt, u_opt_ind] = ValueIteration(P, G) computes the optimal cost and
%   the optimal control input for each state of the state space.
%
%   Input arguments:
%
%       P:
%           A (K x K x L)-matrix containing the transition probabilities
%           between all states in the state space for all control inputs.
%           The entry P(i, j, l) represents the transition probability
%           from state i to state j if control input l is applied.
%
%       G:
%           A (K x L)-matrix containing the stage costs of all states in
%           the state space for all control inputs. The entry G(i, l)
%           represents the cost if we are in state i and apply control
%           input l.
%
%   Output arguments:
%
%       J_opt:
%       	A (K x 1)-matrix containing the optimal cost-to-go for each
%       	element of the state space.
%
%       u_opt_ind:
%       	A (K x 1)-matrix containing the index of the optimal control
%       	input for each element of the state space.

% put your code here
K = size(P,1);
L = size(P,3);

% Initialize costs and optimal policy
J = ones(K,1);
FVal = ones(K,1);
costToGo = zeros(K,1);
cost_to_minimize = zeros(1,L);
% Iterate until cost has converged
err = 1e-5;
iter = 0;
while (1)
    iter = iter + 1;
    for i = 1:K
        for j = 1:L
            cost_to_minimize(1,j) = G(i,j) + P(i,:,j)*J(:) ;
        end
        [costToGo(i),FVal(i)] = min(cost_to_minimize);

    end

    % Check if cost has converged
    if (max(abs(J-costToGo))/max(abs(costToGo)) < err)
        % update cost and break
        J_opt = costToGo(:);
        u_opt_ind = FVal(:);
        break;
    else
        % update cost
        J = costToGo;
    end
end


end