function fig = plotStates(xx_current, media, influencer, I1, I2);
%graphical display of configuration after each timestep
fig = figure(1);
hold off
clf
plot(xx_current(1,I1),xx_current(2,I1),'bo'); 
hold on
plot(xx_current(1,I2),xx_current(2,I2),'ro');
plot(influencer(1,:),influencer(2,:),'go');
plot(media(1,:),media(2,:),'k+');
axis([-2 2 -2 2]);