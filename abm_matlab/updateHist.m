function [histInd, histInf, histMed] = updateHist(xx_current, media, influencer, I1, I2 ,In1, In2, In3, In4, gridx, gridy, histInd, histInf, histMed)
    In = {In1, In2, In3, In4};
    I = {I1, I2};
    %individuals
    for i = 1:4
        for j = 1:2
            entries = intersect(I{j}, In{i});
            counts = histcounts2(xx_current(1,entries),xx_current(2,entries), gridx, gridy);
            histInd(j,i,:,:)= squeeze(histInd(j,i,:,:))+counts;
        end
    end

    %influencer
    for i=1:4
        counts = histcounts2(influencer(1,i),influencer(2,i), gridx, gridy);
        histInf = histInf + counts;
    end

    %media
    for j=1:2
        counts = histcounts2(media(1,j),media(2,j), gridx, gridy);
        histMed = histMed + counts;
    end        
  