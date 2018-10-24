function r = rho(x,u)

global MFcores rho_option xr ur Qdiag Pdiag wrapflag rho_offset

aux = ones(size(x,1),1);
if wrapflag, x(:,1) = abs(x(:,1)); end;       % the top position can be approached from both sides
xr1 = aux*xr(:)';
ur1 = aux*ur(:)';

if rho_option == 1
    r = rho_offset - (xr1 - x).^2*Qdiag' - (ur1 - u).^2*Pdiag';
elseif rho_option == 2
    r = rho_offset - abs(xr1 - x)*Qdiag' - abs(ur1 - u)*Pdiag';
elseif rho_option == 3
    r = rho_offset - sqrt(abs(xr1 - x))*Qdiag' - sqrt(abs(ur1 - u))*Pdiag';
else
    error('Unknown reward option.');
end;
% r = max(0, r/100+1);