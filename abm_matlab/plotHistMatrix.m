function fig = plotHistMatrix(gridcentersx, gridcentersy, histInd, histInf, histMed)
    fig=figure(1);
    clf
    c = 0;
    states =[-1,1];
    for i = 1:4
        for j = 1:2
            c = c+1;
            subplot(5,2,c); 
            surf(gridcentersx, gridcentersy,squeeze(histInd(j,i,:,:)));
            title(['state=' int2str(states(j)) ',Influencer ' int2str(i)])
            set(gca,'XTickLabel',[]);
            view(2); colorbar
        end
    end

    subplot(5,2,c+1)
    surf(gridcentersx, gridcentersy,histInf);
    title(['Influencer']);
    view(2); colorbar;


    subplot(5,2,c+2)
    surf(gridcentersx, gridcentersy,histMed);
    title(['Media']);
    view(2); colorbar;