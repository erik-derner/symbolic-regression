function [xnext,urep] = f(x,u,flag)
% for a single state, x is a row vector, 
%   otherwise it is a matrix with the individual state vectors in rows
% for a single input u is a row vector (or scalar for SI systems), 
%   otherwise it is a matrix with individual input vectors in rows

global Ts eomflag wrapflag velflag

x(:,1) = x(:,1);

if nargin < 3, flag = 0; end;

N = size(x,1);
if flag == 0,                   % calculate next state for all states and all inputs
    for i = size(u,1) : -1 : 1
        for j = N : -1 : 1
%             y = ode4_ti('eompend',0:Ts/2:Ts,[x(j,:)';u(i,:)']);
            y = [x(j,:)';u(i,:)'];
            k1 = eompend(y);
            k2 = eompend(y+0.25*Ts*k1);
            k3 = eompend(y+0.25*Ts*k2);
            k4 = eompend(y+0.5*Ts*k3);
            y = y + (Ts/12)*(k1+2*k2+2*k3+k4);
            k1 = eompend(y);
            k2 = eompend(y+0.25*Ts*k1);
            k3 = eompend(y+0.25*Ts*k2);
            k4 = eompend(y+0.5*Ts*k3);
            y = y + (Ts/12)*(k1+2*k2+2*k3+k4);
            if velflag == 2
                y(2) = (y(1)-x(j,1))/Ts;    % approximate velocity as backward difference
            end;
            if wrapflag,
%                 if y(1) > pi,  y(1) = y(1) - 2*pi; end; % wrapping
%                 if y(1) < -pi, y(1) = y(1) + 2*pi; end; % wrapping
                if y(1) > 2*pi,  y(1) = y(1) - 2*pi; end; % wrapping
                if y(1) < 0, y(1) = y(1) + 2*pi; end; % wrapping
            end;
            xnext((i-1)*N+j,:) = y(1:2)';
        end;
    end;
    urep = repelem(u,N,1);
    %xnext(:,1) = xnext(:,1) - 2.0*pi*(xnext(:,1) > pi) + 2.0*pi*(xnext(:,1) < -pi);
else
    for j = N : -1 : 1
%         y = ode4_ti('eompend',0:Ts/2:Ts,[x(j,:)';u(j,:)']);
        y = [x(j,:)';u(j,:)'];
        k1 = eompend(y);
        k2 = eompend(y+0.25*Ts*k1);
        k3 = eompend(y+0.25*Ts*k2);
        k4 = eompend(y+0.5*Ts*k3);
        y = y + (Ts/12)*(k1+2*k2+2*k3+k4);
        k1 = eompend(y);
        k2 = eompend(y+0.25*Ts*k1);
        k3 = eompend(y+0.25*Ts*k2);
        k4 = eompend(y+0.5*Ts*k3);
        y = y + (Ts/12)*(k1+2*k2+2*k3+k4);
        if velflag == 2
            y(2) = (y(1)-x(j,1))/Ts;    % approximate velocity as backward difference
        end;
        if wrapflag,
%             if y(1) > pi,  y(1) = y(1) - 2*pi; end; % wrapping
%             if y(1) < -pi, y(1) = y(1) + 2*pi; end; % wrapping
              if y(1) > 2*pi,  y(1) = y(1) - 2*pi; end; % wrapping
              if y(1) < 0, y(1) = y(1) + 2*pi; end; % wrapping
        end;
        xnext(j,:) = y(1:2)';
    end;
    urep = u;    
end;