function fig = plotHistMatrix(gridcentersx, gridcentersy, histInd, histInf, histMed)
    fig=figure(1);
    clf
    count = 0;
    states =[-1,1];
    for i = 1:4
        for j = 1:2
            count = count+1;
            subplot(5,2,count); 
            surf(gridcentersx, gridcentersy,squeeze(histInd(j,i,:,:))');
            title(['state=' int2str(states(j)) ',Influencer ' int2str(i)])
            set(gca,'XTickLabel',[]);
            view(2); colorbar
        end
    end

    subplot(5,2,count+1)
    surf(gridcentersx, gridcentersy,histInf');
    title(['Influencer']);
    view(2); colorbar;


    subplot(5,2,count+2)
    surf(gridcentersx, gridcentersy,histMed');
    title(['Media']);
    view(2); colorbar;