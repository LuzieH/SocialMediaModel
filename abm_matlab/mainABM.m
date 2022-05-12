clear all
% Prepare grid for distribution over many simulations
hGrid = 0.2;
Nsim=10;
maxNr = fix(4/hGrid);
Distri1 = zeros(maxNr,maxNr);
Distri2 = zeros(maxNr,maxNr);
DistriI = zeros(maxNr,maxNr);
x1=-2+hGrid/2:hGrid:2-hGrid/2; %TODO: better with full interval (-2,2)
x2=-2+hGrid/2:hGrid:2-hGrid/2;

% repetition of ABM simulations
for count=1:Nsim
    
% initialize number and opinions of agents
n =128;
x=4*rand(2,n)-2*ones(2,n); %random initial conditions

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
state = rand(1,n)-0.5; %check if still correct!
state = sign(state);
I1 = find(state==-1);
I2 = find(state==1);

% initialization of interaction network between agents 
%Net = InitialNet(n,x);  % if based on closeness of opinions
Net = ones(n,n);  % everyone connected to everyone, CHECK selfconnections??

% setting model and simulation parameters
sigma=0.5; % noise on individual agents %CHECK
sigmahat=0; % noise on influencers
sigmatilde=0; %noise on media
dt=0.01;  % simulation stepsize
NT=100;  % simulation time T=NT*dt
a=1; % interaction strength between agents CHECK
b=2; % interaction strength between agents and influencers
c=2; % interaction strength between agents and media

xx=zeros(2,n,NT);
xx(:,:,1)=x;  % xx will store the entire trajectory of opinions

% Display of initial configuration
figure(1)
clf
hold on
plot(x(1,I1),x(2,I1),'bo');
plot(x(1,I2),x(2,I2),'ro');
plot(media(1,:),media(2,:),'k+');
plot(influencer(1,:),influencer(2,:),'go');
axis([-2 2 -2 2]);
if count<9
    print('-dpng',['img/initial' int2str(count) '.png']);
end

% performing the simulation loop
for k=2:NT
    % opinions change due to opinions of friends, influencers and media
    force = a * attraction(xx(:,:,k-1),Net,n) + influence(xx(:,:,k-1),media,influencer,followers,n,state,b,c); %first term deleted
    xx(:,1:n,k) = xx(:,1:n,k-1) + dt*force(:,1:n) + sqrt(dt*sigma)*randn(2,n); % new individual opinions
    
    % influencer opinions adapt slowly to opinions of followers with friction
    % depending on number of followers
    masscenter(:,1) = sum(xx(:,In1,k)')'/length(In1);
    masscenter(:,2) = sum(xx(:,In2,k)')'/length(In2);
    masscenter(:,3) = sum(xx(:,In3,k)')'/length(In3);
    masscenter(:,4) = sum(xx(:,In4,k)')'/length(In4);
    for i=1:4
        if sum(followers(i,:))>0 
            frictionI = 25; %%%%% CHECK
            influencer(:,i) =  influencer(:,i)  + dt/frictionI * (-influencer(:,i)+masscenter(:,i)) + 1/frictionI*sqrt(dt*sigmahat)*randn(2,1);
        end
    end
    
    % media opinions change very slowly based on opinions of followers with friction
    % depending on number of followers
    massmedia(:,1) = sum(xx(:,I1,k)')'/length(I1); %CHECK k-1??
    massmedia(:,2) = sum(xx(:,I2,k)')'/length(I2);
    friction = 100; %%%%%
    for i=1:2
        media(:,i) = media(:,i)  + dt/friction * (-media(:,i)+massmedia(:,i)) + 1/friction * sqrt(dt*sigmatilde)*randn(2,1);
    end
    % exp(-dt/friction) * (media(:,1)-massmedia(:,1)) + massmedia(:,1); %CHECK Friction Term
    

    %% individual may jump from one influencer to another
    % jumps according to rate model
    [followers,In1,In2,In3,In4] = ChangeInfluencerNetwork2(state,xx(:,:,k),n,followers,influencer,dt);

    %graphical display of configuration after each timestep
    figure(1)
    hold off
    clf
    plot(xx(1,I1,k),xx(2,I1,k),'bo'); hold on
    plot(xx(1,I2,k),xx(2,I2,k),'ro');
    %plot(xx(1,In1,k),xx(2,In1,k),'k+');
    plot(influencer(1,:),influencer(2,:),'go');
    plot(media(1,:),media(2,:),'k+');
    axis([-2 2 -2 2]);
    pause(0.1)
    
  %% displaying influencer statistics: individual agents are colored wrt the influencer they are following
figure(37)
clf
In1b=find(state(In1)==-1);
subplot(4,2,1); plot(xx(1,In1(In1b),k),xx(2,In1(In1b),k),'bo',media(1,1),media(2,1),'k+',influencer(1,1),influencer(2,1),'go'); title('state=-1,Influencer 1');
axis([-2 2 -2 2]);
In1r=find(state(In1)==1);
subplot(4,2,2); plot(xx(1,In1(In1r),k),xx(2,In1(In1r),k),'ro',media(1,2),media(2,2),'k+',influencer(1,1),influencer(2,1),'go'); title('state=+1,Influencer 1');
axis([-2 2 -2 2]);

In2b=find(state(In2)==-1);
subplot(4,2,3); plot(xx(1,In2(In2b),k),xx(2,In2(In2b),k),'bo',media(1,1),media(2,1),'k+',influencer(1,2),influencer(2,2),'go'); title('state=-1,Influencer 2');
axis([-2 2 -2 2]);
In2r=find(state(In2)==1);
subplot(4,2,4); plot(xx(1,In2(In2r),k),xx(2,In2(In2r),k),'ro',media(1,2),media(2,2),'k+',influencer(1,2),influencer(2,2),'go'); title('state=+1,Influencer 2');
axis([-2 2 -2 2]);

In3b=find(state(In3)==-1);
subplot(4,2,5); plot(xx(1,In3(In3b),k),xx(2,In3(In3b),k),'bo',media(1,1),media(2,1),'k+',influencer(1,3),influencer(2,3),'go'); title('state=-1,Influencer 3');
axis([-2 2 -2 2]);
In3r=find(state(In3)==1);
subplot(4,2,6); plot(xx(1,In3(In3r),k),xx(2,In3(In3r),k),'ro',media(1,2),media(2,2),'k+',influencer(1,3),influencer(2,3),'go'); title('state=+1,Influencer 3');
axis([-2 2 -2 2]);

In4b=find(state(In4)==-1);
subplot(4,2,7); plot(xx(1,In4(In4b),k),xx(2,In4(In4b),k),'bo',media(1,1),media(2,1),'k+',influencer(1,4),influencer(2,4),'go'); title('state=-1,Influencer 4');
axis([-2 2 -2 2]);
In4r=find(state(In4)==1);
subplot(4,2,8); plot(xx(1,In4(In4r),k),xx(2,In4(In4r),k),'ro',media(1,2),media(2,2),'k+',influencer(1,4),influencer(2,4),'go'); title('state=+1,Influencer 4');
axis([-2 2 -2 2]);
pause(0.1)
end %(end of simulation loop)


% the following lines do only produce statistics across the final configurations of all simulations

%update histogram of opinions
[Distri1,Distri2,DistriI] = updateDistribution(xx(:,:,k),influencer,n,I1,I2,hGrid,Distri1,Distri2,DistriI);

if count<9
    figure(1)
    print('-dpng',['img/final' int2str(count) '.png']); 
    xxx = xx(:,:,NT);
    save(['data/finalstate' int2str(count) '.mat'],'xxx','influencer','media','followers','state','Net','In1','In2','In3','In4');
    figure(37)
    print('-dpng',['img/finalInfluencerDistri' int2str(count) '.png']); 
end

figure(10)
clf
surf(x1,x2,Distri1');
title(['c=' int2str(count)]);
view(2);
colorbar
print('-dpng','img/histblue.png'); 
figure(11)
clf
surf(x1,x2,Distri2');
view(2);
colorbar
title(['c=' int2str(count)]);
print('-dpng','img/histred.png'); 
figure(12)
clf
surf(x1,x2,DistriI');
view(2);
colorbar
title(['c=' int2str(count)]);
print('-dpng','img/histinf.png'); 

end