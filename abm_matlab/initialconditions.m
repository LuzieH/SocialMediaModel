function [x, media, In1, In2, In3, In4, influencer, followers, state, I1, I2, Net] = initialconditions(n)
    x = 4*rand(2,n)-2*ones(2,n); %random initial conditions

    %initialize media and influencer data and follower/influencer relationship
    media =[-1 1; -1 1];
    In1 = find((x(1,:)>0) & (x(2,:)>0));
    In2 = find((x(1,:)<=0) & (x(2,:)>0));
    In3 = find((x(1,:)>0) & (x(2,:)<=0));
    In4 = find((x(1,:)<=0) & (x(2,:)<=0));
    %initial opinions of influencers
    influencer(:,1) = sum(x(:,In1)')'/length(In1);
    influencer(:,2) = sum(x(:,In2)')'/length(In2);
    influencer(:,3) = sum(x(:,In3)')'/length(In3);
    influencer(:,4) = sum(x(:,In4)')'/length(In4);
    %follower network
    followers=zeros(4,n);
    followers(1,In1)=1;
    followers(2,In2)=1;
    followers(3,In3)=1;
    followers(4,In4)=1;

    % initialization of political attitude state of all agents
    state = rand(1,n)-0.5;  
    state = sign(state);
    I1 = find(state==-1);
    I2 = find(state==1);

    % initialization of interaction network between agents 
    Net = ones(n,n);  % everyone connected to everyone including self-connections