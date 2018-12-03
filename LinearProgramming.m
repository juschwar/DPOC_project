function [ J_opt, u_opt_ind ] = LinearProgramming( P, G )
%LINEARPROGRAMMING Value iteration
%   Solve a stochastic shortest path problem by Linear Programming.
%
%   [J_opt, u_opt_ind] = LinearProgramming(P, G) computes the optimal cost
%   and the optimal control input for each state of the state space.
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
A=zeros(K*L,K);
b=zeros(K*L,1);
f = -1*ones(K,1);
for i=1:L
    A(1+(i-1)*K:K+(i-1)*K,:) = eye(K)-P(:,:,i);
    b(1+(i-1)*K:K+(i-1)*K,1) = G(:,i);
end
J_opt = linprog(f,A,b);
cost_to_minimize=zeros(K,L);
for i=1:L
    cost_to_minimize(:,i) = G(:,i)+P(:,:,i)*J_opt;
end
[~,u_opt_ind] = min(cost_to_minimize,[],2);
end

