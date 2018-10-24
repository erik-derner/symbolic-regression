function dx = eompend(x)
% Equations of motion for DCSC inverted pendulum setup, including Coulomb
% friction

global coordflag   

J = 1.7937e-4;  % Pendulum inertia
M = 5.5e-2;     % Pendulum mass
l = 4.2e-2;     % Pendulum length
b = 1.94e-5;    % Viscous damping
c = 8.50e-4;    % Coulomb friction
K = 5.36e-2;    % Torque constant
R = 9.5;        % Rotor resistance
g = 9.81;       % acceleration due to gravity
 
dx = [ x(2);
    (K/R*x(3) + M*g*l*sin(x(1))*(2*coordflag-1) - (b + K^2/R)*x(2) - c*sign(x(2)))/J;
    0];

