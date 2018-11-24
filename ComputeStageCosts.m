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
global pool_num_time_steps detected_additional_time_steps
P = ComputeTransitionProbabilities( stateSpace, controlSpace, ...
        map, gate, mansion, cameras );
K = length(stateSpace(:,1));
L = length(controlSpace(:,1));
G = zeros(K,L);
ofsx = [0 -1 0 1];
ofsy = [1 0 -1 0];

[ind,gate_idx] = max(stateSpace(:,1)==gate(1)&stateSpace(:,2)==gate(2));

for k=1:K
    x1 = stateSpace(k,1); %x-coordinate of the k-th state
    y1 = stateSpace(k,2); %y-coordinate of the k-th state
    for l=1:4
        %Depending on ofsx and ofsy find the state to the north, west,
        %south, east. If this state does not exist the ind yields zero, else
        %one. We get the index idx of the state to the n,w,s,e.
        [ind,idx] = max(stateSpace(:,1)==x1+ofsx(l)&stateSpace(:,2)==y1+ofsy(l));
        %We only transist to the state if it exists
        if(ind)
            x2 = x1+ofsx(l);
            y2 = y1+ofsy(l);
            ts = 1;
            %If the paparazzi lands in the pool this takes us more time steps
            if(map(y2,x2)<0)
                ts = pool_num_time_steps;
            end
            %If we get caught this takes us additional timesteps. As we
            %transist first to another cell before we get caught we also
            %have the timesteps to transist to this next cell
            G(k,l) = G(k,l)+P(k,gate_idx,l)*(detected_additional_time_steps+ts);
            %In case the paparazzi does not get caught in the new state,
            %this just takes us the timesteps to transist
            G(k,l) = G(k,l)+P(k,idx,l)*ts;
        else
            %In case the paparazzi wants to move n,w,s,e but the way is
            %blocked he stays at the current state. If he gets caught this
            %costs him again some timesteps
            G(k,l) = G(k,l)+P(k,gate_idx,l)*(detected_additional_time_steps+1);
            %If the Paparazzi does not get caught, he just stays at the
            %current state and this costs him 1 timestep
            G(k,l) = G(k,l)+P(k,k,l)*1.0;
        end
    end
    %In case the paparazzi takes a picture but the celebrity is not on the
    %picture and he gets caught, this takes him the time to take the photo 
    %plus the time to get to the gate
    G(k,5) = G(k,5)+P(k,gate_idx,5)*(6+1);
    %If the paparazzi takes a picture but the celebrity is not on the
    %picture but he does not get caught either
    G(k,5) = G(k,5)+P(k,k,5)*1.0;
    %If the paparazzi takes successfully a picture
    G(k,5) = G(k,5)+(1-P(k,k,5)-P(k,gate_idx,5))*1;
end


end
