%% Using `ode45` to propage the two-body problem to find out the true anomaly at the end
%
% OFT-Example-02 (based on Curtis 2021 Example 3.1)
%   A geocentric elliptical orbit has a perigee radius of 9600 km and an 
%   apogee radius of 21,000 km. 
%   Calculate the true anomaly 3 hour after perigee. 


%% constants
muEarth = 398600; % [km^3/s^2]

%% parameters
rp = 9600; % [km]
ra = 21000; % [km]

%% using formulas to calculate essential parameters for trajectory equation.
% semi-major
a = (rp + ra) / 2; % [km]
% eccentricity
ecc = (ra - rp) / (ra + rp);
% orbital energy
epsilon = - muEarth / 2 / a; % [km^2/s^2]
% perigee and apogee velocities
vp = sqrt( 2 * (epsilon + muEarth / rp ) ); % [km/s]
va = sqrt( 2 * (epsilon + muEarth / ra ) ); % [km/s]
% period
T = 2*pi/sqrt(muEarth) * a^(3/2); % [s]
T/3600


%% 01. analytical solutions
t = 3*3600; % [s]
Me = 2*pi/T*t
E = MeanToEccentricAnomaly(ecc, Me)
theta = 2 * atand(sqrt(ecc+1)/sqrt(1-ecc) * tan(E/2)) + 360

%% 02. numerical integration for a give time interval
rpVec = [rp, 0, 0].'; % [km]
vpVec = [0, vp, 0].'; % [km/s]

% the eccentricity vector
eVec = cross(vpVec, cross(rpVec, vpVec)) / muEarth - rpVec / norm(rpVec);

% propagation duration
tspan = [0, 3] * 3600; % [s] % this is OFT-Example-02 on slides 
% tspan = [0, 1.132] * 3600; % [s] % this leads to ~120 deg. OFT-Example-0 on slides
% tspan = [0, 4077.04]; % [s] % this leads to ~120 deg. OFT-Example-0 on slides

% initial value
X0 = [rpVec; vpVec]; % semi-colon connects vertically

% integration options, control the accuracy and outputs
options = odeset('RelTol',1e-6, 'AbsTol',1e-6, 'Stats','on'); % try lower the accuracy and see what happens?

% call ode45
[TT, XX] = ode45(@(t,X)OdeTwoBody(muEarth, X), tspan, X0, options);

% calculate true anomaly at the end directly
rEndVec = XX(end, 1:3).';
vEndVec = XX(end, 4:6).';
thetaEnd = CalculateTrueAnomaly(rEndVec, vEndVec, eVec);
% display results
num2str(rad2deg(thetaEnd), '%.2f')
fprintf('# HP: the final true anomaly is %.2f deg.\n', rad2deg(thetaEnd));

% %% visualize the trajectory
figure(53);
clf;
plot(XX(:,1), XX(:,2), 'r.-')
hold on;
% entire elliptic orbit
if ecc < 1
    [tmpTT, tmpXX] = ode45(@(t,X)OdeTwoBody(muEarth, X), [0,2*pi/sqrt(muEarth)*a^(3/2)], X0, options);
    plot(tmpXX(:,1), tmpXX(:,2), ':k')
end
% focus
plot(0, 0, 'rx', 'MarkerSize',15)
% rEnd
plot([0, rEndVec(1)], [0, rEndVec(2)], 'b-');
% rBegin
plot([0, rpVec(1)], [0, rpVec(2)], 'b-');
% angle
tmpN = 10;
tmpAngles = linspace(0, thetaEnd, 2*tmpN+1);
tmpCosA = 0.2*rp*cos(tmpAngles);
tmpSinA = 0.2*rp*sin(tmpAngles);
plot(tmpCosA, tmpSinA, 'b-');
tmpString1 = ['\theta = ' num2str(rad2deg(thetaEnd), '%.2f') ' deg'];
tmpString2 = ['t = ' num2str(TT(end), '%.2f') ' s'];
text(tmpCosA(tmpN), tmpSinA(tmpN), {tmpString2, tmpString1}, 'HorizontalAlignment','left', 'VerticalAlignment','bottom')
%
axis equal;
grid on;
title(['$e=$ ' num2str(ecc, '%.8f') '; only \texttt{tspan} used for \texttt{ode45}.'], 'Interpreter','latex')





%% 03. numerical integration with an event to stop at precisely a given true anomaly

% desired true anomaly
trueAnomalyDesired = deg2rad(120); % [rad]

% large enough time interval
period = 2*pi/sqrt(muEarth)*a^(3/2); % [s]
tspan = [0, 1] * period; % [s]

% initial value
X0 = [rpVec; vpVec]; % semi-colon connects vertically

% integration options, control the accuracy and outputs
%   an `event` handle is added now
optionsWithEvent = odeset('RelTol',1e-6, 'AbsTol',1e-6, 'Stats','on', 'Events',@(a,b)OdeEventExactTrueAnomaly(a,b,eVec,trueAnomalyDesired)); % try lower the accuracy and see what happens?

% call ode45
[TTWithEvent, XXWithEvent] = ode45(@(t,X)OdeTwoBody(muEarth, X), tspan, X0, optionsWithEvent);

% calculate true anomaly at the end directly
rWithEventEndVec = XXWithEvent(end, 1:3).';
vWithEventEndVec = XXWithEvent(end, 4:6).';
thetaWithEventEnd = CalculateTrueAnomaly(rWithEventEndVec, vWithEventEndVec, eVec);
% display results
num2str(rad2deg(thetaWithEventEnd), '%.2f')
fprintf('# HP: the final true anomaly is %.2f deg.\n', rad2deg(thetaWithEventEnd));


% %% visualize the trajectory
figure(54);
clf;
plot(XXWithEvent(:,1), XXWithEvent(:,2), 'r.-')
hold on;
% entire orbit
[tmpTT, tmpXX] = ode45(@(t,X)OdeTwoBody(muEarth, X), [0,2*pi/sqrt(muEarth)*a^(3/2)], X0, options);
plot(tmpXX(:,1), tmpXX(:,2), ':k')
% focus
plot(0, 0, 'rx', 'MarkerSize',15)
% rEnd
plot([0, XXWithEvent(end, 1)], [0, XXWithEvent(end, 2)], 'b-');
% rBegin
plot([0, rpVec(1)], [0, rpVec(2)], 'b-');
% angle
tmpN = 10;
tmpAngles = linspace(0, trueAnomalyDesired, 2*tmpN+1);
tmpCosA = 0.2*rp*cos(tmpAngles);
tmpSinA = 0.2*rp*sin(tmpAngles);
plot(tmpCosA, tmpSinA, 'b-');
tmpString1 = ['\theta_{\rm desired} = ' num2str(rad2deg(trueAnomalyDesired), '%.2f') ' deg'];
tmpString2 = ['\theta_{\rm integrated} = ' num2str(rad2deg(thetaWithEventEnd), '%.2f') ' deg'];
tmpString3 = ['t = ' num2str(TTWithEvent(end), '%.2f') ' s'];
text(tmpCosA(tmpN), tmpSinA(tmpN), {tmpString3, tmpString2, tmpString1}, 'HorizontalAlignment','left', 'VerticalAlignment','bottom')
%
axis equal;
grid on;
title(['$e=$ ' num2str(ecc, '%.8f') '; \texttt{tspan} and \texttt{event} used for \texttt{ode45}.'], 'Interpreter','latex')





%% 04. numerical integration with an event to stop at the next periapsis or apoapsis

% pick a random initial condition from the previously integrated trajectory
rng(2024) % control random number generator seed
idRandom = randi(size(XXWithEvent, 1));
t0 = TTWithEvent(idRandom);
X0 = XXWithEvent(idRandom, :);

% desired true anomaly
styApsisType = 'periapsis';
styApsisType = 'apoapsis';

% large enough time interval
period = 2*pi/sqrt(muEarth)*a^(3/2); % [s]
tspan = [t0, 2*period]; % [s]

% integration options, control the accuracy and outputs
%   an `event` handle is added now
optionsWithEvent = odeset('RelTol',1e-6, 'AbsTol',1e-6, 'Stats','on', 'Events',@(a,b)OdeEventApsis(a,b,styApsisType));

% call ode45
[TTWithEventApsis, XXWithEventApsis] = ode45(@(t,X)OdeTwoBody(muEarth, X), tspan, X0, optionsWithEvent);

% calculate true anomaly at the end directly
rWithEventEndVec = XXWithEventApsis(end, 1:3).';
vWithEventEndVec = XXWithEventApsis(end, 4:6).';
thetaWithEventApsisEnd = CalculateTrueAnomaly(rWithEventEndVec, vWithEventEndVec, eVec);
% display results
num2str(rad2deg(thetaWithEventApsisEnd), '%.2f')
fprintf('# HP: the final true anomaly is %.2f deg.\n', rad2deg(thetaWithEventApsisEnd));


% %% visualize the trajectory
figure(55);
clf;
plot(XXWithEventApsis(:,1), XXWithEventApsis(:,2), 'r.-')
hold on;
% entire orbit
[tmpTT, tmpXX] = ode45(@(t,X)OdeTwoBody(muEarth, X), [0,2*pi/sqrt(muEarth)*a^(3/2)], X0, options);
plot(tmpXX(:,1), tmpXX(:,2), ':k')
% focus
plot(0, 0, 'rx', 'MarkerSize',15)
% rEnd
plot([0, XXWithEventApsis(end, 1)], [0, XXWithEventApsis(end, 2)], 'b-');
% rBegin
plot([0, X0(1)], [0, X0(2)], 'b-');
% angle
tmpN = 10;
thetaBegin = CalculateTrueAnomaly(X0(1:3), X0(4:6), eVec);
tmpAngles = linspace(thetaBegin, thetaWithEventApsisEnd, 2*tmpN+1);
tmpCosA = 0.2*rp*cos(tmpAngles);
tmpSinA = 0.2*rp*sin(tmpAngles);
plot(tmpCosA, tmpSinA, 'b-');
tmpString2 = ['\Delta \theta_{\rm integrated} = ' num2str(rad2deg(thetaWithEventApsisEnd - thetaBegin), '%.2f') ' deg'];
tmpString3 = ['\Delta t = ' num2str(TTWithEventApsis(end) - t0, '%.2f') ' s'];
text(tmpCosA(tmpN), tmpSinA(tmpN), {tmpString3, tmpString2}, 'HorizontalAlignment','left', 'VerticalAlignment','bottom')
%
axis equal;
grid on;
title(['$e=$ ' num2str(ecc, '%.4f') '; \texttt{tspan} and ' styApsisType ' \texttt{event} used for \texttt{ode45}.'], 'Interpreter','latex')


