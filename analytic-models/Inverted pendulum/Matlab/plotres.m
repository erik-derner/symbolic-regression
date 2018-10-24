% clear all;

% filename = 'near_perfect2_RT';
% filename = 'near_perfect2_Coordflag_1_RT';
% filename = 'rt-sr-1';

% D = dir('*sim.mat');
% for jj = 1 : length(D),
% filename = D(jj).name(1:end-4);

printflag = 0;
%load(filename, 'x', 'u', 'Ts', 'xr', 'ur', 'rho_option', 'rho_offset', 'Qdiag', 'Pdiag', 'wrapflag', 'gamma'); printflag = 1;
u = u(:);
filename(filename == '_') = '-';            % replace underscore by dash for nice plotting

% unclear, what the setting was in the different experiments, so this may have to be revised
% Pdiag = 0;

set(0,'DefaultAxesFontName','times');
set(0,'DefaultTextFontName','times');
set(0,'DefaultAxesFontSize',14);
set(0,'DefaultTextFontSize',14);
%set(0,'DefaultAxesLineWidth',0.1);
offs = 0.0;
lw = 1;

% calculate reward and return
u(isnan(u)) = [];
x(isnan(x(:,1)),:) = [];
r = rho(x(2:end,:), u(1:end-1));
R = sum(r.*(gamma.^[0:length(r)-1]'));
rmse = rms(xr(1)-x(:,1));

if coordflag == 1
    x(:,1) = x(:,1) - pi;
end;
% unwrap x for nice plotting
for k = 2 : size(x,1)
    dx = x(k,1) - x(k-1,1);
    if abs(dx) > pi
        x(k,1) = x(k,1) + (dx < 0)*2*pi - (dx > 0)*2*pi;
    end;
end;

% just for plotting convenience
if ~isnan(u(end))
    x = [x;NaN*x(end,:)];
    u = [u;NaN*u(end,:)];
    xtick = [0:0.5:5]';
end;

t = Ts*(0:size(x,1)-1)';
ref = xr(1)*ones(size(t));

figure(1); clf;
subplot(211);
stairs(t,ref,'r','LineWidth',1); hold on;
stairs(t,-ref,'r','LineWidth',1);
stairs(t,x(:,1),'b','LineWidth',lw); hold off;
% axis([0 max(t) 1.1*[-xr(1) xr(1)]]);
grid on
pos = get(gca,'position');
set(gca,'position',[pos(1)+offs pos(2) pos(3)-2*offs  pos(4)]);
xlabel('Time [s]');
ylabel('Angle [rad]');
set(gca,'XTick',xtick)
title([filename ', Discounted return: ' num2str(R) ', RMSE: ' num2str(rmse)])

subplot(212);
stairs(t,u,'b','LineWidth',lw);
axis([0 max(t) 1.1*[min(u) max(u)]]);
grid on
pos = get(gca,'position');
set(gca,'position',[pos(1)+offs pos(2) pos(3)-2*offs  pos(4)]);
xlabel('Time [s]');
ylabel('Control input [V]');
% set(gca,'XTickLabel','')

if printflag
    print('-djpeg', filename);
end;

% end;