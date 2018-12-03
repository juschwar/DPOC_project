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
global p_c gamma_p pool_num_time_steps
K = length(stateSpace(:,1));
L = length(controlSpace(:,1));
[M, N] = size(map);
F = length(mansion(:,1));
H = length(cameras(:,1));
P = zeros(K,K,L);
Pr_gc = zeros(M,N);
ofsx = [0 -1 0 1];
ofsy = [1 0 -1 0];

%Get probability of getting caught by iterating over all cameras
for h=1:H
    n=cameras(h,1); %get x-coordinate
    m=cameras(h,2); %get y-coordinate
    gamma=cameras(h,3); %get quality of the camera
    for l=1:4
        ind = true; %As soon as there is an obstacle, this indicator is set to false
        m1=m+ofsy(l); %We go in y-direction from the camera until an obstacle arives
        n1=n+ofsx(l); %We go in x-direction from the camera until an obstacle arives
        while(ind&&m1>0&&n1>0&&m1<M+1&&n1<N+1) %We also stop as soon as we reach the end of the map
            if map(m1,n1)>0
                ind = false; %The indicator is set to false when an obstacle arives
            else
                %We have to calculate the probability of not getting caught
                %by all cameras that see that state and then calculate from
                %that the probablity of getting caught
                Pr_gc(m1,n1)=1.0-(1.0-Pr_gc(m1,n1))*(1.0-gamma/(abs(m1-m+n1-n)));
            end
            %We go one step further in the direction the camera sees
            m1=m1+ofsy(l);
            n1=n1+ofsx(l);
        end
    end
end

%Plot of the propabilities of getting caught by the cameras
% figure
% h=bar3(Pr_gc,1);  
% shading interp
% for i = 1:length(h)
%      zdata = get(h(i),'Zdata');
%      set(h(i),'Cdata',zdata)
%      set(h,'EdgeColor','k')
% end
% set(gca, 'Xdir', 'reverse')
% title('Probability of getting caught')
% xlabel x
% ylabel y


%Calculation of the probability of taking a picture of celebrity in the
%mansion. This does not yet include the rare case of taking picture outside
%of the mansion. This is the same problem that the paparazzi is seen by the
%mansion and this is the same problem as beeing caught by the camera.
Pr_tp = zeros(M,N);
%We iterate over every part of the mansion
for f=1:F
    n=mansion(f,1); %x-coordinate of the part of the mansion
    m=mansion(f,2); %y-coordinate of the part of the mansion
    gamma=gamma_p;
    for l=1:4
        ind = true;
        n1=n;
        m1=m;
        m1=m1+ofsy(l);
        n1=n1+ofsx(l);
        while(ind&&m1>0&&n1>0&&m1<M+1&&n1<N+1)
            if map(m1,n1)>0
                ind = false;
            else
                %It is not possible to be in line of sight with several
                %parts of the mansion
                Pr_tp(m1,n1)=gamma/(abs(m1-m+n1-n));
            end
            m1=m1+ofsy(l);
            n1=n1+ofsx(l);
        end
    end
end

%Plot the probability of taking picture of the celebrity inside the mansion
% figure 
% h=bar3(Pr_tp,1);
% shading interp
% for i = 1:length(h)
%      zdata = get(h(i),'Zdata');
%      set(h(i),'Cdata',zdata)
%      set(h,'EdgeColor','k')
% end
% set(gca, 'Xdir', 'reverse')
% title('Probability of taking picture')
% xlabel x
% ylabel y

%Get Index of the gate in stateSpace
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
            %The probability of getting caught is 1 - the probability of
            %not getting caught over all time steps
            P(k,gate_idx,l) = 1.0-(1-Pr_gc(y2,x2))^ts;
            %In case the paparazzi does not get caught in the new state, he
            %transists to the new state
            P(k,idx,l) = P(k,idx,l)+1.0-P(k,gate_idx,l);
        else
            %In case the paparazzi wants to move n,w,s,e but the way is
            %blocked he stays at the current state. And the probability of
            %getting caught is that of the current state.
            P(k,gate_idx,l) = Pr_gc(y1,x1);
            %If the Paparazzi does not get caught, he just stays at the
            %current state
            P(k,k,l) = P(k,k,l)+1-P(k,gate_idx,l);
        end
    end
    %In case the paparazzi takes a picture, the probability of getting
    %caught is that of the current state.
    P(k,gate_idx,5) = Pr_gc(y1,x1)*(1-max(p_c,Pr_tp(y1,x1)));
    %The Probability of staying at the current state is 1 - the probability
    %of taking a picture and finishing the task - The probability of
    %getting caught.
    P(k,k,5) = P(k,k,5)+1.0-max(p_c,Pr_tp(y1,x1))-P(k,gate_idx,5);
end


% Only leave this assert functions as long as we develop the code

%The transition probabilities have to sum to one if there is no possibility
%to finish the task. As we can transist to the terminal state when we take
%a picture this does not hold for the control input 5.
if(sum(sum(P(:,:,1),2).*sum(P(:,:,2),2).*sum(P(:,:,3),2).*sum(P(:,:,4),2)<1)>0)
    save('errorMap.mat','map');
    save('errorP.mat','P');
    save('errorSs.mat','stateSpace');
end
    
assert(sum(sum(P(:,:,1),2).*sum(P(:,:,2),2).*sum(P(:,:,3),2).*sum(P(:,:,4),2)<1)==0)

%All Probabilities have to be greater than 0
assert(sum((P(:,:,1)>=0).*(P(:,:,2)>=0).*(P(:,:,3)>=0).*(P(:,:,4)>=0)<1,'all')==0)





end
