%% Using `ode45` to propage the two-body problem to find out the true anomaly at the end
%
%   OFT-Example-03 (simplified Curtis 2021 Example 3.5)
%

clear; clc;

%% constants
muEarth = 398600; % [km^3/s^2]
rE = 6378; % [km]

%% parameters
rp = 300 + rE; % [km]
vp = 15; % [km/s]

% other parameters deduced
epsilon = 0.5*vp^2 - muEarth/rp;
h = rp * vp;
ecc = sqrt(1 + 2 * epsilon * h^2 / muEarth^2);

%% 01. numerical integration for a give time interval
rpVec = [rp, 0, 0].'; % [km]
vpVec = [0, vp, 0].'; % [km/s]

% the eccentricity vector
eVec = cross(vpVec, cross(rpVec, vpVec)) / muEarth - rpVec / norm(rpVec);

% propagation duration
tspan = [0, 3] * 3600; % [s]

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
fprintf('# HP: the final true anomaly is %.2f deg.\n', thetaEnd);

% %% visualize the trajectory
figure(73);
clf;
plot(XX(:,1), XX(:,2), 'r.-')
hold on;
% focus
plot(0, 0, 'rx', 'MarkerSize',15)
% rEnd
plot([0, rEndVec(1)], [0, rEndVec(2)], 'b-');
text(rEndVec(1)*0.5, rEndVec(2)*0.5, ['r = ' num2str(norm(rEndVec), '%.2f') ' km'], 'HorizontalAlignment','right', 'VerticalAlignment','bottom')
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





%% 02. numerical integration with an event to stop at precisely a given true anomaly

% desired true anomaly
trueAnomalyDesired = deg2rad(100); % [rad]

% large enough time interval
tspan = [0, 24] * 3600; % [s]

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
fprintf('# HP: the final true anomaly is %.2f deg.\n', thetaWithEventEnd);


% %% visualize the trajectory
figure(74);
clf;
plot(XXWithEvent(:,1), XXWithEvent(:,2), 'r.-')
hold on;
% entire orbit
if ecc < 1
    [tmpTT, tmpXX] = ode45(@(t,X)OdeTwoBody(muEarth, X), tspan, X0, options);
    plot(tmpXX(:,1), tmpXX(:,2), ':k')
end
% focus
plot(0, 0, 'rx', 'MarkerSize',15)
% rEnd
plot([0, XXWithEvent(end, 1)], [0, XXWithEvent(end, 2)], 'b-');
text(rWithEventEndVec(1)*0.5, rWithEventEndVec(2)*0.5, ['r = ' num2str(norm(rWithEventEndVec), '%.2f') ' km'], 'HorizontalAlignment','right', 'VerticalAlignment','bottom')
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
title(['$e=$ ' num2str(ecc, '%.8f')], 'Interpreter','latex')




