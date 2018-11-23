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
[K, L] = size(G);
V1 = zeros(L,1);
V = zeros(K,1);
u_opt_ind = zeros(K,1);
error = 100;
eps = 1e-5;
while(error>eps)
    V_old = V;
    for k1=1:K
        for l = 1:L
            V1(l) = G(k1,l);
            for k2=1:K
                V1(l) = V1(l) + P(k1,k2,l)*V(k2);
            end
        end
        [V(k1), u_opt_ind(k1)] = min(V1);
    end
    [error, I] = max(V-V_old);      
end
J_opt = V;


end