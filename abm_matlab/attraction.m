function f = attraction(x,Net,n)
for j=1:n
    J=find(Net(j,:)==1);
    if isempty(J)
        f(:,j)=[0 0]';
    else
        ff = [0 0]'; ww=0;
        for m=1:length(J)
            d = norm(x(:,j)-x(:,J(m)));
            w = exp(-d/2);
            ff = ff + w*(-x(:,J(m))+x(:,j));
            ww = ww+w;
        end
        f(:,j) = ff/ww;
    end

end