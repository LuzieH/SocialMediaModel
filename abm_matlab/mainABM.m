clear all

%%% PARAMETERS

% Prepare grid for distribution over many simulations
hGrid = 0.2;
maxNr = fix(4/hGrid);
gridx=-2:hGrid:2;
gridy=-2:hGrid:2;
gridcentersx = -2+hGrid/2:hGrid:2-hGrid/2;
gridcentersy = -2+hGrid/2:hGrid:2-hGrid/2;
%histograms for Individuals, Influencer and Media
histInd = zeros(2,4,maxNr, maxNr); 
histInf = zeros(maxNr, maxNr);
histMed = zeros(maxNr, maxNr);

% setting model and simulation parameters
Nsim=100;
dt=0.01;  % simulation stepsize
NT=100;  % simulation time frame T=NT*dt
n = 128;
sigma=0.5; % noise on individual agents 
sigmahat=0; % noise on influencers
sigmatilde=0; % noise on media
a=1; % interaction strength between agents 
b=2; % interaction strength between agents and influencers
c=2; % interaction strength between agents and media
frictionI = 25;  % friction for influencers
friction = 100; % friction for media
eta = 50; % rate constant for changing influencer 

% repetition of ABM simulations
for count=1:Nsim
    % initial conditions
    [x, media, In1, In2, In3, In4, influencer, followers, state, I1, I2, Net] = initialconditions(n);

    xx=zeros(2,n,NT);
    xx(:,:,1)=x;  % xx will store the entire trajectory of opinions

    % Display of initial configuration
    if count<9
        %fig = plotStates(xx(:,:,1), media, influencer, I1, I2); %plot all states in one figure
        fig = plotStatesMatrix(xx(:,:,1), media, influencer, I1, I2 ,In1, In2, In3, In4);
        print(fig,[ 'img/initial' int2str(count)],'-dpng');
    end

    % performing the simulation loop
    for k=2:NT
        % opinions change due to opinions of friends, influencers and media
        force = a * attraction(xx(:,:,k-1),Net,n) + influence(xx(:,:,k-1),media,influencer,followers,n,state,b,c);
        xx(:,1:n,k) = xx(:,1:n,k-1) + dt*force(:,1:n) + sqrt(dt*sigma)*randn(2,n); % new individual opinions
        % note that there are no boundary conditions, agents could escape [-2,2]x[-2,2]


        % influencer opinions adapt slowly to opinions of followers with friction
        % depending on number of followers
        In = {In1, In2, In3, In4};
        for i=1:4
            if sum(followers(i,:))>0 
                masscenter(:,i) = sum(xx(:,In{i},k-1)')'/length(In{i});
                influencer(:,i) =  influencer(:,i)  + dt/frictionI * (-influencer(:,i)+masscenter(:,i)) + 1/frictionI*sqrt(dt*sigmahat)*randn(2,1);
            end
        end
        
        % media opinions change very slowly based on opinions of followers with friction
        % depending on number of followers
        I = {I1, I2};
        for i=1:2
            massmedia(:,i) = sum(xx(:,I{i},k-1)')'/length(I{i});
            media(:,i) = media(:,i)  + dt/friction * (-media(:,i)+massmedia(:,i)) + 1/friction * sqrt(dt*sigmatilde)*randn(2,1);
        end
        

        %% individual may jump from one influencer to another
        % jumps according to rate model
        [followers,In1,In2,In3,In4] = ChangeInfluencerNetwork2(state,xx(:,:,k-1),n,followers,influencer,dt, eta);

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
        fig = plotStatesMatrix(xxx, media, influencer, I1, I2 ,In1, In2, In3, In4);
        print(fig, '-dpng',['img/final' int2str(count) '.png']); 
    end

    % update ensemble histogram
    [histInd, histInf, histMed] = updateHist(xx(:,:,k), media, influencer, I1, I2 ,In1, In2, In3, In4, gridx, gridy, histInd, histInf, histMed);

    disp("Simulation " + int2str(count) + " finished")
end

%plot ensemble histogram
fig = plotHistMatrix(gridcentersx, gridcentersy, histInd, histInf, histMed);
print(fig, '-dpng',['img/ensemble.png']); 

%save ensemble histograms
save(['data/ensemble_finalhist.mat'], 'histInd', 'histInf', 'histMed','n', 'Nsim');