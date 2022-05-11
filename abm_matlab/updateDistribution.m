function [d1,d2,dI] = updateDistribution(x,influencer,n,I1,I2,hGrid,Distri1,Distri2,DistriI)

d1=Distri1;
d2=Distri2;
maxNr = fix(4/hGrid);
boxnr = boxNummer(x,n,hGrid);

for i=1:length(I1)
    j=I1(i);
    if (boxnr(1,j)>=1)&&(boxnr(1,j)<=maxNr)&& (boxnr(2,j)>=1)&&(boxnr(2,j)<=maxNr)
        d1(boxnr(1,j),boxnr(2,j)) = d1(boxnr(1,j),boxnr(2,j)) + 1;
    end
end
for i=1:length(I2)
    j=I2(i);
    if (boxnr(1,j)>=1)&&(boxnr(1,j)<=maxNr)&& (boxnr(2,j)>=1)&&(boxnr(2,j)<=maxNr)
        d2(boxnr(1,j),boxnr(2,j)) = d2(boxnr(1,j),boxnr(2,j)) + 1;
    end
end

dI=DistriI;
boxnrI = boxNummerI(influencer,hGrid);
for i=1:4
    if (boxnrI(1,i)>=1)&&(boxnrI(1,i)<=maxNr)&& (boxnrI(2,i)>=1)&&(boxnrI(2,i)<=maxNr)
        dI(boxnrI(1,i),boxnrI(2,i)) = dI(boxnrI(1,i),boxnrI(2,i)) + 1;
    end
end



