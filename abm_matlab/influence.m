function force = influence(x,media,influencer,B,n,state,b,c)
force1 =zeros(size(x));
force2 =zeros(size(x));
for j=1:n
    if state(j)==1
        force1(:,j) = -x(:,j)+media(:,2);
    else %when state==-1
        force1(:,j) = -x(:,j)+media(:,1);
    end
    for kk=1:4
        if B(kk,j)==1
        force2(:,j)=-x(:,j)+influencer(:,kk);
        end
    end
end
force = c*force1 + b*force2;