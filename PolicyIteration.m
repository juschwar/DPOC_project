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
[K, L] = size(G);

mu_hat = 5*ones(K,1);
mu_hat_old = mu_hat*0;

G_hat = zeros(K,1);
P2 = zeros(K);
eps = 1e-5;
error=100;
J_mu=zeros(K,1);
while isequal(mu_hat,mu_hat_old)~=1
    J_mu_old = J_mu;
    for k1=1:K
        G_hat(k1) = G(k1,mu_hat(k1));
        P2(k1,:) = P(k1,:,mu_hat(k1));
    end
    J_mu = (eye(K)-P2)\G_hat;

    mu_hat_old = mu_hat;
    for k1=1:K
        J_hat = G(k1,:);
        for u=1:L
            J_hat(u) = J_hat(u)+P(k1,:,u)*J_mu(:);
        end
        [~,mu_hat(k1)] = min(J_hat);
    end
    error = max(abs(J_mu-J_mu_old));
end
u_opt_ind = mu_hat;
J_opt = J_mu;

end

