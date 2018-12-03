function G = ComputeStageCosts( stateSpace, controlSpace, map, gate, mansion, cameras )
%COMPUTESTAGECOSTS Compute stage costs.
% 	Compute the stage costs for all states in the state space for all
%   control inputs.
%
%   G = ComputeStageCosts(stateSpace, controlSpace, map, gate, mansion,
%   cameras) computes the stage costs for all states in the state space
%   for all control inputs.
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
%       G:
%           A (K x L)-matrix containing the stage costs of all states in
%           the state space for all control inputs. The entry G(i, l)
%           represents the expected stage cost if we are in state i and 
%           apply control input l.

% put your code here
global p_c gamma_p pool_num_time_steps detected_additional_time_steps
K = size(stateSpace,1);
L = size(controlSpace,1);
P = ComputeTransitionProbabilities( stateSpace, controlSpace, ...
        map, gate, mansion, cameras );
[~, idx_gate] = max(stateSpace(:,1)==gate(1) & stateSpace(:,2)==gate(2));
     
G = zeros(K,L);
for i = 1:K
    n = stateSpace(i,1);
    m = stateSpace(i,2);
    % n
    [v_to, idx_to] = max(stateSpace(:,1)==n & stateSpace(:,2)==m+1);
    if (v_to~=0)
        if (map(m+1,n)==0)
            G(i,1) = P(i,idx_to,1)*1 + ...
                P(i,idx_gate,1) * (detected_additional_time_steps+1);
        else
            G(i,1) = P(i,idx_to,1)*pool_num_time_steps + ...
                P(i,idx_gate,1) * (pool_num_time_steps + detected_additional_time_steps);
        end
    else
            G(i,1) = P(i,i,1)*1 + ...
                P(i,idx_gate,1) * (1 + detected_additional_time_steps);
    end
    
    % w
    [v_to, idx_to] = max(stateSpace(:,1)==n-1 & stateSpace(:,2)==m);
    if (v_to~=0)
        if (map(m,n-1)==0)
            G(i,2) = P(i,idx_to,2)*1 + ...
                P(i,idx_gate,2) * (detected_additional_time_steps+1);
        else
            G(i,2) = P(i,idx_to,2)*pool_num_time_steps + ...
                P(i,idx_gate,2) * (pool_num_time_steps + detected_additional_time_steps);
        end
    else
            G(i,2) = P(i,i,2)*1 + ...
                P(i,idx_gate,2) * (1 + detected_additional_time_steps);
    end
    
    % s
    [v_to, idx_to] = max(stateSpace(:,1)==n & stateSpace(:,2)==m-1);
    if (v_to~=0)
        if (map(m-1,n)==0)
            G(i,3) = P(i,idx_to,3)*1 + ...
                P(i,idx_gate,3) * (detected_additional_time_steps+1);
        else
            G(i,3) = P(i,idx_to,3)*pool_num_time_steps + ...
                P(i,idx_gate,3) * (pool_num_time_steps + detected_additional_time_steps);
        end
    else
            G(i,3) = P(i,i,3)*1 + ...
                P(i,idx_gate,3) * (1 + detected_additional_time_steps);
    end
    
    % e
    [v_to, idx_to] = max(stateSpace(:,1)==n+1 & stateSpace(:,2)==m);
    if (v_to~=0)
        if (map(m,n+1)==0)
            G(i,4) = P(i,idx_to,4)*1 + ...
                P(i,idx_gate,4) * (detected_additional_time_steps+1);
        else
            G(i,4) = P(i,idx_to,4)*pool_num_time_steps + ...
                P(i,idx_gate,4) * (pool_num_time_steps + detected_additional_time_steps);
        end
    else
            G(i,4) = P(i,i,4)*1 + ...
                P(i,idx_gate,4) * (1 + detected_additional_time_steps);
    end
    
    %u=p
    G(i,5) = P(i,i,5)*1 + P(i,idx_gate,5)*(1+detected_additional_time_steps) + ...
        (1-P(i,i,5)-P(i,idx_gate,5))*1;
end

end
