function [value,isterminal,direction] = OdeEventExactTrueAnomaly(t, X, eVec, trueAnomalyDesired)
%% The event that the integration reaches the desired true anomaly.
%
% Inputs:
%   t: time (not used in this function indeed)
%   X: [x y z vx vy vz]^T, row vector of position and velocity [km km/s]
%   eVec: used to calculate the true anomlay
%   trueAnomalyDesired: [rad] desired true anomlay, must be in [0, 2*pi)
%   
% Outputs:
%   dXdt: derivative
%
% Author: Hao Peng (hao.peng@erau.edu)

%% calculate the current true anomaly
% % restore to r, v vectors
% rVec = X(1:3).';
% vVec = X(4:6).';
% % using dot-product to get the true anomlay
% thetaEnd = acos(dot(rVec, eVec) / norm(rVec) / norm(eVec));
% % test whether in upper or lower half of the ellipse
% if dot(rVec, vVec) < 0
%     thetaEnd = 2*pi - thetaEnd; % adjust if in the lower half
% end

%% calculate the current true anomaly using a function call
% restore to r, v vectors
rVec = X(1:3).';
vVec = X(4:6).';
thetaEnd = CalculateTrueAnomaly(rVec, vVec, eVec);

%% outputs
% the value that ode45 monitors
value = mod(thetaEnd, 2*pi) - mod(trueAnomalyDesired, 2*pi);
% stop the integration at the detection of the event, otherwise the integration continues but just logs the event
isterminal = 1;
% 1: even only when the value increases; 
% -1: even only when value decreases
% 0: even for both cases
direction = 0;

end