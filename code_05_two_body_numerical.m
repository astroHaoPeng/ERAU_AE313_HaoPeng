%% Using `ode45` to propage the two-body problem to find out the true anomaly at the end

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

%% numerical integration for a give time interval
rpVec = [rp, 0, 0].'; % [km]
raVec = [ra, 0, 0].'; % [km]
vpVec = [0, vp, 0].'; % [km/s]
vaVec = [0, -va, 0].'; % [km/s]

% propagation duration
% tspan = [0, 3] * 3600; % [s]
% tspan = [0, 1.132] * 3600; % [s] % this leads to ~120 deg.
% tspan = [0, TTWithEvent(end)]; % [s] % this leads to more accurate 120 deg. (Must be called after getting `TTWithEvent`.)

% initial value
X0 = [rpVec; vpVec]; % semi-colon connects vertically

% integration options, control the accuracy and outputs
options = odeset('RelTol',1e-6, 'AbsTol',1e-6, 'Stats','on'); % try lower the accuracy and see what happens?

% call ode45
[TT, XX] = ode45(@(t,X)OdeTwoBody(muEarth, X), tspan, X0, options);

%% calculate true anomaly at the end directly
% restore to r, v vectors
rEndVec = XX(end, 1:3).';
vEndVec = XX(end, 4:6).';
% unit vector of apex line, also the eccentricity vector direction
iHatVec = [1, 0, 0];
% using dot-product to get the true anomlay
thetaEnd = acos(dot(rEndVec, iHatVec) / norm(rEndVec));
% test whether in upper or lower half of the ellipse
if dot(rEndVec, vEndVec) < 0
    thetaEnd = 2*pi - thetaEnd; % adjust if in the lower half
end
% display results
num2str(rad2deg(thetaEnd), '%.2f')
fprintf('# HP: the final true anomaly is %.2f deg.\n', thetaEnd);

%% visualize the trajectory
figure(33);
clf;
plot(XX(:,1), XX(:,2), 'r.-')
hold on;
% entire orbit
[tmpTT, tmpXX] = ode45(@(t,X)OdeTwoBody(muEarth, X), [0,2*pi/sqrt(muEarth)*a^(3/2)], X0, options);
plot(tmpXX(:,1), tmpXX(:,2), ':k')
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
title(['$e=$ ' num2str(ecc, '%.8f')], 'Interpreter','latex')





%% numerical integration with an event to stop at precisely a given true anomaly

% desired true anomaly
trueAnomalyDesired = deg2rad(120); % [rad]

% large enough time interval
period = 2*pi/sqrt(muEarth)*a^(3/2); % [s]
tspan = [0, 1] * period; % [s]

% integration options, control the accuracy and outputs
%   an `event` handle is added now
optionsWithEvent = odeset('RelTol',1e-6, 'AbsTol',1e-6, 'Stats','on', 'Events',@(a,b)eventExactTrueAnomaly(a,b,trueAnomalyDesired)); % try lower the accuracy and see what happens?

% call ode45
[TTWithEvent, XXWithEvent] = ode45(@(t,X)OdeTwoBody(muEarth, X), tspan, X0, optionsWithEvent);

%% visualize the trajectory
figure(34);
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
tmpString1 = ['\theta = ' num2str(rad2deg(trueAnomalyDesired), '%.2f') ' deg'];
tmpString2 = ['t = ' num2str(TTWithEvent(end), '%.2f') ' s'];
text(tmpCosA(tmpN), tmpSinA(tmpN), {tmpString2, tmpString1}, 'HorizontalAlignment','left', 'VerticalAlignment','bottom')
%
axis equal;
grid on;
title(['$e=$ ' num2str(ecc, '%.8f')], 'Interpreter','latex')





%% private functions
function dXdt = OdeTwoBody(muEarth, X)

r = X(1:3);
v = X(4:6);

drdt = v;

rCube = (sum(r.^2)).^(3/2);
dvdt = - muEarth / rCube * r;

dXdt = [drdt(:); dvdt(:)]; % trick: A(:) will change any vector to a column vector, equivalent to `reshape(A, [], 1)`.
end


function [value,isterminal,direction] = eventExactTrueAnomaly(t, X, trueAnomalyDesired)
% restore to r, v vectors
rVec = X(1:3).';
vVec = X(4:6).';
% unit vector of apex line, also the eccentricity vector direction
iHatVec = [1, 0, 0];
% using dot-product to get the true anomlay
thetaEnd = acos(dot(rVec, iHatVec) / norm(rVec));
% test whether in upper or lower half of the ellipse
if dot(rVec, vVec) < 0
    thetaEnd = 2*pi - thetaEnd; % adjust if in the lower half
end
% outputs
value = thetaEnd - trueAnomalyDesired;
isterminal = 1;
direction = 0;
end
