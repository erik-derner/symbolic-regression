function mu = tmpgrade(x,c,BFtype)

if BFtype == 1 || BFtype == 6
    n1 = ones(length(x),1); 
    x = x(:,ones(1,length(c)-1));
    c1 = c(n1,1:end-1); c2 = c(n1,2:end); d = c2 - c1;
    mu = max(0,min([n1,(x-c1)./d],[(c2-x)./d,n1]));
elseif BFtype == 2
    B = -(0.5*mean(diff(c)))^2/log(0.5);
    for i = 1 : length(c)
        mu(:,i) = exp(-(x-c(i)).^2/B);
    end;
    mu = mu./(sum(mu,2)*ones(1,length(c)));
elseif BFtype == 3
    c2 = (c(2) - c(1))/2;
    for i = 1 : length(c)
        mu(:,i) = 0.5*(x > c(i)-2*c2).*(x <= c(i)-c2).*((x-c(i)+2*c2)/c2).^2 + ...
                  (x > c(i)-c2).*(x <= c(i)+c2).*(1-0.5*((x-c(i))/c2).^2) + ...
                  0.5*(x > c(i)+c2).*(x <= c(i)+2*c2).*((x-c(i)-2*c2)/c2).^2;
    end;
    mu(x<c(1),1) = 1;
    mu(x>c(end),end) = 1;
else
    error('Unknown type of basis function.');
end;