% save VI animation as a movie 
set(0,'DefaultAxesFontName','times');
set(0,'DefaultTextFontName','times');
set(0,'DefaultAxesFontSize',14);
set(0,'DefaultTextFontSize',14);
set(0,'DefaultAxesLineWidth',0.2);
set(0,'DefaultFigurePaperOrientation','portrait');

% combo parameters
rhoCoeff = 0.845321134919033; maxVsym = -16.0176219294322;      % combo
% rhoCoeff = 0.83863174903311; maxVsym = -16.2195096968864;       % model2
% rhoCoeff = 0.79645329487492; maxVsym = -23.3609972995661;       % model3

options = optimset('Display','iter','TolFun',1e-100,'TolX',1e-20);

Tsim = 3;                       % simulation time
hillflag = 1;                   % set to 0 - policy, 1 - HC on V, 2 - proxy, 3 - combo, 4 - DNN, 5 - spline, 6 - SR approx of V
% U = [-2:0.01:2]';             % refined input domain

% filename = 'Ref=180_Umax=2_Qvel=0,1_Grid=30_gamma=0,98_quad_lin.mat';
%load(filename);

movie_flag = 0;

if movie_flag 
    v = VideoWriter(['Response-V_' filename '.avi']); 
    v.FrameRate = 5;
    open(v); 
end;

thetaind = size(theta,2);
% thetaind = 1;

V(:) = tmprodvect(XX,MFcores,BFtypeV)*theta(:,end);
% V(:) = feval(Vfilename, XX)*theta(:,thetaind);
% V(:) = VFunSNGPRandom(XX);
P(:) = tmprodvect(XX,MFcores,BFtypeP)*p(:,end);
R = 0;
%x(1,:) = [-2.3471 0]';
x(1,:) = [coordflag*pi 0]';
clear u
figure(1); clf; 
subplot(221); contour(xx1,xx2,V,150); hold on; plot(X(:,1),X(:,2),'r.')
title('Value function and state trajectory');
xlabel('Angle [rad]'); ylabel('Angular velocity [rad/s]'); 
subplot(223); contour(xx1,xx2,P,150); hold on; plot(X(:,1),X(:,2),'r.');
title('Policy and state trajectory')
xlabel('Angle [rad]'); ylabel('Angular velocity [rad/s]'); 
for k = 1 : floor(Tsim/Ts)
    if eomflag == 0
        xnext = f(x(k,:),U);
    else
        xnext = fsm(x(k,:),U);
    end;        
    if xnextflag == 1
        reward = rho(xnext,U);
    else
        reward = rho(x(k,:),U);
    end;
% HC
    if hillflag == 1        
        [maxV(k),maxUind] = max(reward + gamma*tmprodvect(xnext,MFcores,BFtypeV)*theta(:,end));
        u(k) = U(maxUind);
% Proxy
    elseif hillflag == 2
        [maxV(k),maxUind] = max(feval(Vfilename,xnext));
        u(k) = U(maxUind);
% Combo
    elseif hillflag == 3
        [maxV(k),maxUind] = max(rhoCoeff*reward + gamma*(model_combo(xnext)-maxVsym));
        u(k) = U(maxUind);
% DNN
    elseif hillflag == 4
%         [maxV(k),maxUind] = max(rho(x(k,:),U) + gamma*graph.compute(xnext));
%         u(k) = U(maxUind);
% Spline
    elseif hillflag == 5
        [maxV,maxUind] = max(reward + gamma*cbsplinen_eval(PHI, TRI-1, BCoefs, d, gridres, (xnext - ones(size(xnext))*diag(minX)) ./ (ones(size(xnext))* diag(rangeX))));
        u(k) = U(maxUind);
% SR approximation of V
    elseif hillflag == 6
        [maxV(k),maxUind] = max(reward + gamma*feval(Vfilename, xnext)*theta(:,thetaind));
        u(k) = U(maxUind);
    else
        u(k) = tmprodvect(x(k,:),MFcores,BFtypeP)*p(:,end);
    end;

        u(k) = min(max(u(k),-Umax),Umax);
        x(k+1,:) = f(x(k,:),u(k));
        xp(k+1,:) = x(k+1,:);                   % just for plotting
        R(k,:) = rho(x(k+1,:), u(k));

        % unwrap data for plot
        if k > 1
            dxp = xp(k,1) - xp(k-1,1);
            if abs(dxp)>pi & coordflag == 0
                xp(k,1) = xp(k,1) + (dxp < 0)*2*pi - (dxp > 0)*2*pi;
            end;
        end;

        subplot(221);
        if k > 1
            plot(x(k-1:k,1),x(k-1:k,2),'k.', 'markersize', 20);
        else
            plot(x(k,1),x(k,2),'k.', 'linewidth', 2, 'markersize', 20);
        end;
        subplot(223);
        if k > 1
            plot(x(k-1:k,1),x(k-1:k,2),'k.', 'markersize', 20);
        else
            plot(x(k,1),x(k,2),'k.', 'linewidth', 2, 'markersize', 20);
        end;

        subplot(222)
        stairs(Ts*(0:k-1), xp(1:k,1), '-', 'linewidth', 2, 'markersize', 10); hold on;
        stairs(Ts*(0:k-1), xr(1)*ones(k,1), 'r'); hold off;
        grid
        if xr(1) == pi
            hold on; stairs(Ts*(0:k-1), -xr(1)*ones(k,1), 'r'); hold off;
        end;
        axis([0 Tsim -4 4]);
        xlabel('Time [s]'); ylabel('Angle [rad]'); title('Time-domain response')
        subplot(224)
        stairs(Ts*(0:k-1), u(1:k), '-', 'linewidth', 2, 'markersize', 10); grid
        axis([0 Tsim -Umax-.2 Umax+.2]);
        xlabel('Time [s]'); ylabel('Control input [V]');
        drawnow
        if movie_flag
            frame = getframe(gcf);
            writeVideo(v,frame);
        end;
        
%         if abs(abs(x(k,1)) - xr(1)) < 0.0001, break; end;
    end;
    x = x(1:length(u),:);

if movie_flag
    frame = getframe(gcf);
    for j = 1 : 1
        writeVideo(v,frame); 
    end;
    close(v);
end;

sum(R)/k