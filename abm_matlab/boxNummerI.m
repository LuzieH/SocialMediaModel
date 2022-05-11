function boxnr = boxNummerI(influencer,h)

for i=1:4
    boxnr(1,i) = fix((influencer(1,i)+2)/h)+1;
    boxnr(2,i) = fix((influencer(2,i)+2)/h)+1;
end