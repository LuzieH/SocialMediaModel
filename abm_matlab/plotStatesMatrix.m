function fig = plotStatesMatrix(xx_current, media, influencer, I1, I2 ,In1, In2, In3, In4)
    fig=figure(1);
    clf
    c = 1;
    In = {In1, In2, In3, In4};
    I = {I1, I2};
    color = {'bo', 'ro'};
    states =[-1,1];
    for i = 1:4
        for j = 1:2
            entries = intersect(I{j}, In{i});
            subplot(4,2,c); plot(xx_current(1,entries),xx_current(2,entries),color{j},media(1,j),media(2,j),'k+',influencer(1,i),influencer(2,i),'go'); title(['state=' int2str(states(j)) ',Influencer ' int2str(i)]);
            axis([-2 2 -2 2]);
            c = c+1;
        end
    end