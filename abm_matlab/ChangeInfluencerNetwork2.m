function [follow,In1,In2,In3,In4] = ChangeInfluencerNetwork2(state,x,n,followers,influencer,dt);

eta = 50;
theta =0.1
%%%Version based on rates
follow = followers;

%compute happiness = fraction of followers with same state
for i=1:4
    fraction(i)= sum(followers(i,:).*state)/sum(followers(i,:));
end

% fraction
% length(find(state==1))
% length(find(state==-1))
%compute distance of followers to influencers
for i=1:4
    for j=1:n
        d = norm(x(:,j)-influencer(:,i));
        distance(j,i)= exp(-d);
    end
end

%compute attractiveness of influencer for followers
for j=1:n
    for i=1:4
        g = state(j)*fraction(i);
        if g<0 
            g=theta;
        else
            g+=theta
        end
        attractive(j,i)= eta * distance(j,i)*g;
    end
    r=rand;
    alpha=-log(1-r); %random number distributed due to exp(1)
    lambda = sum(attractive(j,:)); %total jump rate
    if lambda*dt>alpha %CHECK should be lamda*dt > r ??
        p = attractive(j,:)/lambda;
        r2=rand;
        k=1;
        while sum(p(1:k))<r2
            k=k+1;
        end
        follow(:,j)=[0 0 0 0]';        
        follow(k,j)=1;
    end
end

    
In1 = find(follow(1,:)==1);
In2 = find(follow(2,:)==1);
In3 = find(follow(3,:)==1);
In4 = find(follow(4,:)==1);
followers = follow;