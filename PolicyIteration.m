function [ J_opt, u_opt_ind ] = PolicyIteration( P, G )
%POLICYITERATION Value iteration
%   Solve a stochastic shortest path problem by Policy Iteration.
%
%   [J_opt, u_opt_ind] = PolicyIteration(P, G) computes the optimal cost and
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

Jmu = zeros(K,1);
FVal = 5*ones(K,1);
costToGo = zeros(K,1);
cost_to_minimize = zeros(K,L);
Pmu = zeros(K,K);
Gmu = zeros(K,1);

err = 1e-5;
iter = 0;

while(1)
    Jold = Jmu;
    for i=1:K
        Pmu(i,:) = P(i,:,FVal(i));
        Gmu(i,1) = G(i,FVal(i));
    end
    Jmu = (eye(K)-Pmu)\Gmu;
    
    iter = iter + 1;
        for j = 1:L
            cost_to_minimize(:,j) = G(:,j) + P(:,:,j)*Jmu(:) ;
        end
        [costToGo, FVal] = min(cost_to_minimize(:,:),[],2);
        if ( max(abs(Jmu - Jold))<err)
            J_opt = costToGo;
            u_opt_ind = FVal;
            break
        end

end
    
end
