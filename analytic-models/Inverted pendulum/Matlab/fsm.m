function [xnext,urep] = fsm(x,u,flag)
%FSM Symbolic process model for 1DOF pendulum
% Input:
%  x - [N x 2] current state
%  u - [N x 1] control action (input)
%  flag
% Output:
%  xnext - [N x 2] next state(s) in row(s)
%  urep

global Ts wrapflag velflag

if nargin < 3, flag = 0; end;

N = size(x,1);
if flag == 0,                   % calculate next state for all states and all inputs
    for i = size(u,1) : -1 : 1
        for j = N : -1 : 1
            y = [x(j,:), u(i,:)];
            x1_next = fsm_model_x1(y);
            x2_next = fsm_model_x2(y);
            y = [x1_next; x2_next; y(3)];
            if velflag == 2
                y(2) = (y(1)-x(j,1))/Ts;    % approximate velocity as backward difference
            end;
            if wrapflag,
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
        y = [x(j,:), u(j,:)];
        x1_next = fsm_model_x1(y);
        x2_next = fsm_model_x2(y);
        y = [x1_next; x2_next; y(3)];
        if velflag == 2
            y(2) = (y(1)-x(j,1))/Ts;    % approximate velocity as backward difference
        end;
        if wrapflag,
            if y(1) > 2*pi,  y(1) = y(1) - 2*pi; end; % wrapping
            if y(1) < 0, y(1) = y(1) + 2*pi; end; % wrapping
        end;
        xnext(j,:) = y(1:2)';
    end;
    urep = u;    
end;