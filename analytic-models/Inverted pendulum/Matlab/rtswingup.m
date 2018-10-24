%% flags
hillflag = 1;                   % set to 0 - policy, 1 - HC on BE RHS
Tf = 5;                         % final time
h = Ts;                         % sample time
n = 2;                          % system order
Ae = 0.0;                       % exploration amplitude
N = floor(Tf/Ts);               % number of samples

%% History variables
y = zeros(N,1);                 % initialize output history
x = zeros(N,n);                 % initialize state history
u = zeros(N,1);                 % initialize input history
ys = y;                         % sensor data history
 
%% initialize DCMOTOR - USB init
% sensors to read [status, time, angle, velocity, current, voltage, dangle/dt]
fugiboard('CloseAll');                          % close port
H = fugiboard('Open', 'mops1');                 % open port
H.WatchdogTimeout = 5;                          % watchdog timeout
fugiboard('SetParams', H);                      % set the parameter
fugiboard('Write', H, 1, 1, 0, 0);              % dummy write to sync interface board
fugiboard('Write', H, 0, 1, 0, 0);              % end reset
for i = 1 : 3,                                  % read a couple of times
        D = fugiboard('Read',H);                    % dummy read sensor data
end;

%% Control loop
tic;                                            % reset Matlab's tic-toc timer
for k = 1 : N,                                  % run 
    D = fugiboard('Read',H);                    % read sensor data
    ys(k) = D(3);                               % measured angle
    y(k,:) = mod(ys(k)+2*pi,2*pi);              % wrapping
    x(k,:) = [y(k) D(4)];                       % measured velocity
    x(k,1) = x(k,1) + coordflag*pi;
    if eomflag == 0
        xnext = f(x(k,:),U);                    % use physical model (EOM)
    else
        xnext = fsm(x(k,:),U);                  % use symbolic model
    end;        
    if hillflag == 1
        [maxV(k),maxUind] = max(rho(xnext,U) + gamma*tmprodvect(xnext,MFcores,BFtypeV)*theta(:,end));
        u(k) = U(maxUind);
    else
        u(k) = tmprodvect(x(k,:),MFcores,BFtypeP)*p(:,end);
    end;
    u(k) = u(k) + Ae*randn;                    % add exploration noise
    fugiboard('Write',H,0,1,u(k),0);           % send control input to process    
    tt = toc;
    while tt < Ts, tt = toc; end;              % synchronize with real time
    tic;
end;
% reset system, close the port
fugiboard('Write', H, 0, 1, .5, 0.0);           % send the pendulum down
pause(1);                                       % wait a while
fugiboard('Write', H, 0, 1, 0, 0.0);            % reset input to zero
pause(2);                                       % wait for the swings to damp
fugiboard('CloseAll');                          % close the port

plotres  