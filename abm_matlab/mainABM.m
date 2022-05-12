clear all
% Prepare grid for distribution over many simulations
hGrid = 0.2;
Nsim=100;
maxNr = fix(4/hGrid);
gridx=-2:hGrid:2;
gridy=-2:hGrid:2;
gridcentersx = -2+hGrid/2:hGrid:2-hGrid/2;
gridcentersy = -2+hGrid/2:hGrid:2-hGrid/2;
histInd = zeros(2,4,maxNr, maxNr);
histInf = zeros(maxNr, maxNr);
histMed = zeros(maxNr, maxNr);

% repetition of ABM simulations
for count=1:Nsim
    
    % initialize number and opinions of agents
    n = 128;
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
    %Net = InitialNet(n,x);  % if based on closeness of opinions
    Net = ones(n,n);  % everyone connected to everyone, WHAT ABOUT SELFCONNECTIONS

    % setting model and simulation parameters
    sigma=0.5; % noise on individual agents 
    sigmahat=0; % noise on influencers
    sigmatilde=0; %noise on media
    dt=0.01;  % simulation stepsize
    NT=100;  % simulation time T=NT*dt
    a=1; % interaction strength between agents 
    b=2; % interaction strength between agents and influencers
    c=2; % interaction strength between agents and media

    xx=zeros(2,n,NT);
    xx(:,:,1)=x;  % xx will store the entire trajectory of opinions

    % Display of initial configuration
    if count<9
        %fig = plotStates(xx(:,:,1), media, influencer, I1, I2);
        fig = plotStatesMatrix(xx(:,:,1), media, influencer, I1, I2 ,In1, In2, In3, In4);
        print(fig,[ 'img/initial' int2str(count)],'-dpng');
    end

    % performing the simulation loop
    for k=2:NT
        % opinions change due to opinions of friends, influencers and media
        force = a * attraction(xx(:,:,k-1),Net,n) + influence(xx(:,:,k-1),media,influencer,followers,n,state,b,c); %first term deleted
        xx(:,1:n,k) = xx(:,1:n,k-1) + dt*force(:,1:n) + sqrt(dt*sigma)*randn(2,n); % new individual opinions
        %%%BOUNDARY CONDITIONS FOR AGENTS?
        % influencer opinions adapt slowly to opinions of followers with friction
        % depending on number of followers
        masscenter(:,1) = sum(xx(:,In1,k-1)')'/length(In1);
        masscenter(:,2) = sum(xx(:,In2,k-1)')'/length(In2);
        masscenter(:,3) = sum(xx(:,In3,k-1)')'/length(In3);
        masscenter(:,4) = sum(xx(:,In4,k-1)')'/length(In4);
        for i=1:4
            if sum(followers(i,:))>0 
                frictionI = 25;  
                influencer(:,i) =  influencer(:,i)  + dt/frictionI * (-influencer(:,i)+masscenter(:,i)) + 1/frictionI*sqrt(dt*sigmahat)*randn(2,1);
            end
        end
        
        % media opinions change very slowly based on opinions of followers with friction
        % depending on number of followers
        massmedia(:,1) = sum(xx(:,I1,k-1)')'/length(I1); 
        massmedia(:,2) = sum(xx(:,I2,k-1)')'/length(I2);
        friction = 100; %%%%%
        for i=1:2
            media(:,i) = media(:,i)  + dt/friction * (-media(:,i)+massmedia(:,i)) + 1/friction * sqrt(dt*sigmatilde)*randn(2,1);
        end
        

        %% individual may jump from one influencer to another
        % jumps according to rate model
        [followers,In1,In2,In3,In4] = ChangeInfluencerNetwork2(state,xx(:,:,k-1),n,followers,influencer,dt);

        %graphical display of configuration after each timestep
        fig = plotStates(xx(:,:,k), media, influencer, I1, I2);
        %pause(0.1)
        
        %displaying influencer statistics: individual agents are colored wrt the influencer they are following
        fig = plotStatesMatrix(xx(:,:,k), media, influencer, I1, I2 ,In1, In2, In3, In4);
        %pause(0.1)
    
    end %(end of simulation loop)

    %final distribution
    if count<9
        xxx = xx(:,:,NT);
        %fig = plotStates(xxx, media, influencer, I1, I2);
        fig = plotStatesMatrix(xxx, media, influencer, I1, I2 ,In1, In2, In3, In4);
        print(fig, '-dpng',['img/final' int2str(count) '.png']); 
        save(['data/finalstate' int2str(count) '.mat'],'xxx','influencer','media','followers','state','Net','In1','In2','In3','In4');
    end

    % update ensemble histogram
    [histInd, histInf, histMed] = updateHist(xx(:,:,k), media, influencer, I1, I2 ,In1, In2, In3, In4, gridx, gridy, histInd, histInf, histMed);

end

%plot ensemble histogram
fig = plotHistMatrix(gridcentersx, gridcentersy, histInd, histInf, histMed);
print(fig, '-dpng',['img/ensemble.png']); 