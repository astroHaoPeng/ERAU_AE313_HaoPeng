function [value,isterminal,direction] = OdeEventApsis(t, X, strType)
%% Stop the integration exactly at the desired true anomaly.
%
% Inputs:
%   t: time (not used in this function indeed)
%   X: [x y z vx vy vz]^T, row vector of position and velocity [km km/s]
%   strType: 'periapsis' or 'apoapsis'
%   
% Outputs:
%   value: true anomlay [rad], within [0, 2*pi)
%   isterminal: 0 or 1
%   direction: -1, 1, or 0.
%
% Author: Hao Peng (hao.peng@erau.edu)

%% At apsis, r and v are perpendicular
% restore to r, v vectors
rVec = X(1:3).';
vVec = X(4:6).';
%
tmp = dot(rVec, vVec) / norm(rVec) / norm(vVec);

%% outputs
switch strType
    case 'periapsis'
        % the value that ode45 monitors
        value = tmp;
        % stop the integration at the detection of the event, otherwise the integration continues but just logs the event
        isterminal = 1;
        % 1: even only when the value increases;
        % -1: even only when value decreases
        % 0: even for both cases
        direction = 1;
    case 'apoapsis'
        value = tmp;
        isterminal = 1;
        direction = -1;
end


end