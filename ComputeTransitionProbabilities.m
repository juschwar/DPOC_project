function P = ComputeTransitionProbabilities( stateSpace, controlSpace, map, gate, mansion, cameras )
%COMPUTETRANSITIONPROBABILITIES Compute transition probabilities.
% 	Compute the transition probabilities between all states in the state
%   space for all control inputs.
%
%   P = ComputeTransitionProbabilities(stateSpace, controlSpace,
%   map, gate, mansion, cameras) computes the transition probabilities
%   between all states in the state space for all control inputs.
%
%   Input arguments:
%
%       stateSpace:
%           A (K x 2)-matrix, where the i-th row represents the i-th
%           element of the state space.
%
%       controlSpace:
%           A (L x 1)-matrix, where the l-th row represents the l-th
%           element of the control space.
%
%       map:
%           A (M x N)-matrix describing the terrain of the estate map.
%           Positive values indicate cells that are inaccessible (e.g.
%           trees, bushes or the mansion) and negative values indicate
%           ponds or pools.
%
%   	gate:
%          	A (1 x 2)-matrix describing the position of the gate.
%
%    	mansion:
%          	A (F x 2)-matrix indicating the position of the cells of the
%           mansion.
%
%    	cameras:
%          	A (H x 3)-matrix indicating the positions and quality of the 
%           cameras.
%
%   Output arguments:
%
%       P:
%           A (K x K x L)-matrix containing the transition probabilities
%           between all states in the state space for all control inputs.
%           The entry P(i, j, l) represents the transition probability
%           from state i to state j if control input l is applied.

% put your code here
K = length(stateSpace(:,1));
L = length(controlSpace(:,1));
[M, N] = size(map);
F = length(mansion(:,1));
H = length(cameras(:,1));
P = zeros(K,K,L);
Pr_gc = zeros(M,N);
for h=1:H
    n=cameras(h,1);
    m=cameras(h,2);
    gamma=cameras(h,3);
    ind = true;
    m1=m;
    while(ind&&m1>1)
        m1=m1-1;
        m1
        n
        if map(m1,n)>0
            ind = false;
        else
            Pr_gc(m1,n)=Pr_gc(m1,n)+gamma/(abs(m1-m));
        end
    end
    ind=true;
    m2=m;
    while(ind&&m2<M)
        m2=m2+1;
        if map(m2,n)>0
            ind = false;
        else
            Pr_gc(m2,n)=Pr_gc(m2,n)+gamma/(abs(m2-m));
        end
    end
    ind = true;
    n1=n;
    while(ind&&n1>1)
        n1=n1-1;
        if map(m,n1)>0
            ind = false;
        else
            Pr_gc(m,n1)=Pr_gc(m,n1)+gamma/(abs(n-n1));
        end
    end
    ind=true;
    n2=n;
    while(ind&&n2<N)
        n2=n2+1;
        if map(m,n2)>0
            ind = false;
        else
            Pr_gc(m,n2)=Pr_gc(m,n2)+gamma/(abs(n2-n));
        end
    end
    
end
figure
bar3(Pr_gc)  

Pr_tp = zeros(M,N);
for f=1:F
    n=mansion(f,1);
    m=mansion(f,2);
    gamma=gamma_p;
    ind = true;
    m1=m;
    while(ind&&m1>1)
        m1=m1-1;
        m1
        n
        if map(m1,n)>0
            ind = false;
        else
            Pr_tp(m1,n)=gamma/(abs(m1-m));
        end
    end
    ind=true;
    m2=m;
    while(ind&&m2<M)
        m2=m2+1;
        if map(m2,n)>0
            ind = false;
        else
            Pr_tp(m2,n)=gamma/(abs(m2-m));
        end
    end
    ind = true;
    n1=n;
    while(ind&&n1>1)
        n1=n1-1;
        if map(m,n1)>0
            ind = false;
        else
            Pr_tp(m,n1)=gamma/(abs(n-n1));
        end
    end
    ind=true;
    n2=n;
    while(ind&&n2<N)
        n2=n2+1;
        if map(m,n2)>0
            ind = false;
        else
            Pr_tp(m,n2)=gamma/(abs(n2-n));
        end
    end
    
end       
figure 
bar3(Pr_tp)




end
