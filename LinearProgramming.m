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
[K, L] = size(G);
f = -1*ones(K,1);
A=zeros(K*L,K);
b=zeros(K*L,1);
for u=1:L
    A(1+(u-1)*K:K+(u-1)*K,:) = eye(K)-P(:,:,u);
    b(1+(u-1)*K:K+(u-1)*K,1) = G(:,u);
end
lb = zeros(K,1);
Aeq = [];
beq = [];
J_opt = linprog(f,A,b,Aeq,beq,lb);
c=zeros(K,L);
for u=1:L
    c(:,u) = G(:,u)+P(:,:,u)*J_opt;
end
[~,u_opt_ind] = min(c.');
end