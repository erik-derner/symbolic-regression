function dof1 = tmprodvect(x,mfs,BFtype)
% Calculate Cartesian product of membership degrees
% for a set of points x (x is a matrix)
global Vfilename

if BFtype == 4,                                     % polynomial
    dof1 = polyregr(x);                             % basis function values for all current states (polynomial)
elseif BFtype == 7,                                 % symbolic
    dof1 = feval(Vfilename,x);
else
    [n1,n2] = size(x);                              % aux var
    n3 = 1;
    for j = n2 : -1 : 1,
        Mu{n2-j+1} = tmpgrade(x(:,j),mfs{j},BFtype); % pre-computer membership degrees
        n3 = n3*length(mfs{j});
    end;
    
%     dof1 = sparse(n3,n3);
    for i = n1 : -1 : 1,
        dof = Mu{end}(i,:)';
        for j = n2-1 : -1 : 1,
            dof = dof*Mu{j}(i,:);                   % calculate product
            dof = dof(:);                           % flatten in an array
        end;
        dof1(i,:) = dof';
    end;
    dof1 = sparse(dof1);
end;