function A=InitialNet_bots(n,m,x,Net)
A=zeros(n+m,n+m);
A(1:n,1:n)=Net;
for j=n+1:n+m
    for i=1:n+m
        if i~=j
            dij = norm(x(:,j)-x(:,i));
            if dij<1
                A(i,j)=1;
            end
        end
    end
end
for j=1:n+m
    for i=n+1:n+m
        if i~=j
            dij = norm(x(:,j)-x(:,i));
            if dij<1
                A(i,j)=1;
            end
        end
    end
end
% A(1,2)=1; A(1,3)=1;
% A(n,n-1)=1; A(n,n-2)=1;
% for i=2:n-1
%     A(i,i-1)=1;
%     A(i,i+1)=1;
% end