function A=InitialNet(n,x)
A=zeros(n,n);
for j=1:n
    for i=1:n
        if i~=j
            dij = norm(x(:,j)-x(:,i));
            if dij<0.5
                A(i,j)=1;
            end
        end
    end
end
