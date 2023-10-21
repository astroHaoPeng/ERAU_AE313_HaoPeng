function coeECI = RvToCoe(rVecECI, vVecECI)
%% convert ECI coordinates to classical orbital elements.
%
% Inputs:
%   rVecECI: [km]
%   vVecECI: [km/s]
%
% Outputs:
%   coeECI: classical orbital elements consisting of [h, e, theta, RAAN, i, AP]
%
% Author: Hao Peng (hao.peng@erau.edu)

%% default input 
if nargin == 0
    % testing data during development
    rVecECI = [7706.3257;645.9415;-1468.8057]; % position at t1 [km]
    vVecECI = [-1.1701;6.1123;-3.4510];     % velocity at t1 [km/s]
end

%% constants
muEarth = 398600;         % gravitational parameter [km^3/s^2]

% ECI coordinate of ECI bases
xECI = [1 0 0].';
yECI = [0 1 0].';
zECI = [0 0 1].';

%% calculate
rNorm = norm(rVecECI);    % magnitude of r_vec

% % calculate semi-major axis
% vNorm = norm(vVecECI);    % magnitude of v_vec
% epsilon = vNorm^2/2 - muEarth/rNorm;   % specific mechanical energy [km^2/s^2]
% a = - muEarth / (2*epsilon);      % semi-major axis [km]

hVec = cross(rVecECI, vVecECI);  % specific angular momentum vector [km^2/s]
hNorm = norm(hVec);  % specific angular momentum [km^2/s]

eccVec = cross(vVecECI, hVec) / muEarth - rVecECI/rNorm; % eccentricity vector
eccNorm = norm(eccVec);    % eccentricity

theta = acos(dot(rVecECI, eccVec) / rNorm / eccNorm);
if dot(rVecECI, vVecECI) < 0
    % should be in 3rd/4th quadrant
    theta = 2*pi - theta;
end

NVec = cross(zECI, hVec);   % node line vector

OmegaRAAN = acos(NVec(1) / norm(NVec)); % RAAN Omega [rad], range: [0, 2*pi)
if NVec(2) < 0
    % should be in 3rd/4th quadrant
    OmegaRAAN = 2*pi - OmegaRAAN;
end

inc = acos(dot(hVec, zECI) / norm(hNorm)); % inclination [rad], range: [0, pi)

omegaAP = acos(dot(NVec, eccVec) / norm(NVec) / eccNorm); % Argument of Perigee omega [rad], range: [0, 2*pi)
if eccVec(3) < 0
    % should be in 3rd/4th quadrant
    omegaAP = 2*pi - omegaAP;
end

%% outputs
% a 1x6 row vector
coeECI = [hNorm, eccNorm, theta, OmegaRAAN, inc, omegaAP];

end