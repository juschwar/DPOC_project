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
global pool_num_time_steps p_c gamma_p

% Probabilities of beeing detected
Pc = zeros(size(map));
for i = 1:size(cameras,1)
    n = cameras(i,1);
    m = cameras(i,2);
    quality = cameras(i,3);
    
    %west
    x = map(m,n-1:-1:1)<=0;
    x = [x(:);0];
    ii = 1;
    while(x(ii))
        if(map(m,n-ii)==0)
            Pc(m,n-ii) = 1 - ((1 - Pc(m,n-ii)) * (1 - quality/ii));
        else
            Pc(m,n-ii) = 1 - ((1 - Pc(m,n-ii)) * (1 - quality/ii)^pool_num_time_steps);
        end
        ii = ii + 1;
    end
    %south
    x = map(m-1:-1:1,n)<=0;
    x = [x(:);0];
    ii = 1;
    while(x(ii))
        if(map(m-ii,n)==0)
            Pc(m-ii,n) = 1 - ((1 - Pc(m-ii,n)) * (1 - quality/ii));
        else
            Pc(m-ii,n) = 1 - ((1 - Pc(m-ii,n)) * (1 - quality/ii)^4);
        end
        ii = ii + 1;
    end
    %east
    x = map(m,n+1:end)<=0;
    x = [x(:);0];
    ii = 1;
    while(x(ii))
        if(map(m,n+ii)==0)
            Pc(m,n+ii) = 1 - ((1 - Pc(m,n+ii)) * (1 - quality/ii));
        else
            Pc(m,n+ii) = 1 - ((1 - Pc(m,n+ii)) * (1 - quality/ii)^4);
        end
        ii = ii + 1;
    end
    %north
    x = map(m+1:end,n)<=0;
    x = [x(:);0];
    ii = 1;
    while(x(ii))
        if(map(m+ii,n)==0)
            Pc(m+ii,n) = 1 - ((1 - Pc(m+ii,n)) * (1 - quality/ii));
        else
            Pc(m+ii,n) = 1 - ((1 - Pc(m+ii,n)) * (1 - quality/ii)^4);
        end
        ii = ii + 1;
    end
    
end

% Probabilities of taking a photo
Pp = zeros(size(map));
for i = 1:size(mansion,1)
    n = mansion(i,1);
    m = mansion(i,2);
    quality = gamma_p;
    
    %west
    x = map(m,n-1:-1:1)<=0;
    x = [x(:);0];
    ii = 1;
    while(x(ii))
        Pp(m,n-ii) = max([p_c quality/ii]);
        ii = ii + 1;
    end
    %south
    x = map(m-1:-1:1,n)<=0;
    x = [x(:);0];
    ii = 1;
   while(x(ii))
        Pp(m-ii,n) = max([p_c quality/ii]);
        ii = ii + 1;
    end
    %east
    x = map(m,n+1:end)<=0;
    x = [x(:);0];
    ii = 1;
    while(x(ii))
        Pp(m,n+ii) = max([p_c quality/ii]);
        ii = ii + 1;
    end
    %north
    x = map(m+1:end,n)<=0;
    x = [x(:);0];
    ii = 1;
   while(x(ii))
        Pp(m+ii,n) = max([p_c quality/ii]);
        ii = ii + 1;
    end
    
end

% build the transition probability Matrix
K = size(stateSpace,1);
L = size(controlSpace,1);
P = zeros(K,K,L);

for i = 1:K
     from = stateSpace(i,:);
     [~, idx_gate] = max(stateSpace(:,1)==gate(1) & stateSpace(:,2)==gate(2));
     
     % u=n
     to = from + [0 1];
     [v_to, idx_to] = max(stateSpace(:,1)==to(1) & stateSpace(:,2)==to(2));
     if (v_to~=0)
         if(idx_to ~= idx_gate)
             P(i,idx_gate,1) = Pc(from(2),from(1));
             P(i,idx_to,1) = 1-P(i,idx_gate,1);
         else
             P(i,idx_to,1) = 1;
         end
     else
         if (i~=idx_gate)
             P(i,idx_gate,1) = Pc(from(2),from(1));
             P(i,i,1) = 1-P(i,idx_gate,1);
         else
             P(i,idx_gate,1) = 1;
         end
     end
     
     % u=w
     to = from + [-1 0];
     [v_to, idx_to] = max(stateSpace(:,1)==to(1) & stateSpace(:,2)==to(2));
     if (v_to~=0)
         if(idx_to ~= idx_gate)
             P(i,idx_gate,2) = Pc(from(2),from(1));
             P(i,idx_to,2) = 1-P(i,idx_gate,2);
         else
             P(i,idx_to,2) = 1;
         end
     else
         if (i~=idx_gate)
             P(i,idx_gate,2) = Pc(from(2),from(1));
             P(i,i,2) = 1-P(i,idx_gate,2);
         else
             P(i,idx_gate,2) = 1;
         end
     end
     
     % u=s
     to = from + [0 -1];
     [v_to, idx_to] = max(stateSpace(:,1)==to(1) & stateSpace(:,2)==to(2));
     if (v_to~=0)
         if(idx_to ~= idx_gate)
             P(i,idx_gate,3) = Pc(from(2),from(1));
             P(i,idx_to,3) = 1-P(i,idx_gate,3);
         else
             P(i,idx_to,3) = 1;
         end
     else
         if (i~=idx_gate)
             P(i,idx_gate,3) = Pc(from(2),from(1));
             P(i,i,3) = 1-P(i,idx_gate,3);
         else
             P(i,idx_gate,3) = 1;
         end
     end
     % u=e
     to = from + [1 0];
     [v_to, idx_to] = max(stateSpace(:,1)==to(1) & stateSpace(:,2)==to(2));
     if (v_to~=0)
         if(idx_to ~= idx_gate)
             P(i,idx_gate,4) = Pc(from(2),from(1));
             P(i,idx_to,4) = 1-P(i,idx_gate,4);
         else
             P(i,idx_to,4) = 1;
         end
     else
         if (i~=idx_gate)
             P(i,idx_gate,4) = Pc(from(2),from(1));
             P(i,i,4) = 1-P(i,idx_gate,4);
         else
             P(i,idx_gate,4) = 1;
         end
     end
     %u=p
     if (i~=idx_gate)
         P(i,idx_gate,5) = Pc(from(2),from(1))*(1-Pp(from(2),from(1)));
         P(i,i,5) = 1.0-Pp(from(2),from(1))-P(i,idx_gate,5);  
     else
         P(i,i,5) = 1-Pp(from(2),from(1));
     end
     
end
% tmp = round(sum(P,2),3)
% [v, idx] = min(tmp(:,:,5));
% disp(idx)
% disp(nnz(tmp(:,:,1)~=1))
% disp(nnz(tmp(:,:,2)~=1))
% disp(nnz(tmp(:,:,3)~=1))
% disp(nnz(tmp(:,:,4)~=1))
% disp(nnz(tmp(:,:,5)~=1))

