%% Curtise2020, Problem 4.14
%
% Method-1: Solve the problem by 
%   1. converting ECI RV to COE;
%   2. using time t and true anomaly conversion to add additional duration;
%   3. converting back to a new set of COE (only true anomaly is differetn);
%   4. converting COE to ECI RV.
%
% Method-2: Numerical propagation using ode45

clear all; clc;

muEarth = 398600; % [km^3/s^2]

r1 = [-5000 -8000 -2100].'; % [km/s]
v1 = [-4 3.5 -3].'; % [km/s]
dt = 50 * 60; % [s]

timer1 = tic;

%% Method-1: Step-1
coe1 = ConvertRvToCoe(r1, v1);

%% Method-1: Step-2
% prepare useful variables
h = coe1(1);
e = coe1(2);
a = h^2 / muEarth / (1-e^2); % semi-major axis [km]
T = 2*pi/sqrt(muEarth) * a^(3/2); % period [s] 
n = 2*pi / T; % mean motion [rad/s]

trueAnomaly1 = coe1(3);

eccentricAnomaly1 = 2 * atan( sqrt((1-e)/(1+e)) * tan(trueAnomaly1 / 2) ); % [rad]
meanAnomaly1 = eccentricAnomaly1 - e * sin(eccentricAnomaly1); % Kepler's equation [rad]
t1 = meanAnomaly1 / n; % time since perigee [s]

t2 = t1 + dt;
meanAnomaly2 = n * t2;
eccentricAnomaly2 = MeanToEccentricAnomaly(e, meanAnomaly2);
trueAnomlay2 = 2 * atan( sqrt((1+e)/(1-e)) * tan(eccentricAnomaly2 / 2) ); % [rad]

%% Method-1: Step-3
coe2 = coe1;
coe2(3) = trueAnomlay2;

%% Method-1: Step-4
[r2, v2] = ConvertCoeToRv(coe2);

%% Method-1: display
disp(['# Hao: r1 = [' num2str(r1.', '%+15.4e') '] km'])
disp(['# Hao: v1 = [' num2str(v1.', '%+15.4e') '] km/s'])
disp(['# Hao: dt = ' num2str(dt, '%.2f') ' s']);
disp(['# Hao: r2 = [' num2str(r2.', '%+15.4e') '] km'])
disp(['# Hao: v2 = [' num2str(v2.', '%+15.4e') '] km/s'])
fprintf(['# Hao: --------- ']);
toc(timer1);

%% Method-2
timer2 = tic;
[tt, xx] = ode45(@(t,x)OdeTwoBody(muEarth,x), [0, dt], [r1; v1], odeset('RelTol',1e-8,'AbsTol',1e-8));
r2Numerical = xx(end, 1:3).';
v2Numerical = xx(end, 4:6).';
% disp(['# Hao: r1 = [' num2str(r1.', '%+15.4e') '] km'])
% disp(['# Hao: v1 = [' num2str(v1.', '%+15.4e') '] km/s'])
% disp(['# Hao: dt = ' num2str(dt, '%.2f') ' s']);
disp(['# Hao: t1r1v1 --ode45--> t2r2v2']);
disp(['# Hao: r2 = [' num2str(r2Numerical.', '%+15.4e') '] km'])
disp(['# Hao: v2 = [' num2str(v2Numerical.', '%+15.4e') '] km/s'])
fprintf(['# Hao: --------- ']);
toc(timer2);
