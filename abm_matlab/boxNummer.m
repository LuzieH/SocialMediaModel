function boxnr = boxNummer(x,n,h)

for i=1:n
    boxnr(1,i) = fix((x(1,i)+2)/h)+1;
    boxnr(2,i) = fix((x(2,i)+2)/h)+1;
end