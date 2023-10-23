function [rVecECI, vVecECI] = ConvertCoeToRv(coe)
%% convert ECI coordinates to classical orbital elements.
%
% Inputs:
%   coeECI: classical orbital elements consisting of [h, e, theta, RAAN, i, AP]
%
% Outputs:
%   rVecECI: [km]
%   vVecECI: [km/s]
%
% Author: Hao Peng (hao.peng@erau.edu)

%% default input 
% if nargin == 0
%     % testing data during development
%     rVecECI = [7706.3257;645.9415;-1468.8057]; % position at t1 [km]
%     vVecECI = [-1.1701;6.1123;-3.4510];     % velocity at t1 [km/s]
% end

%% constants
muEarth = 398600;         % gravitational parameter [km^3/s^2]

%$ restore variables
hNorm = coe(1);
eccNorm = coe(2);
theta = coe(3);
OmegaRAAN = coe(4);
inc = coe(5);
omegaAP = coe(6);

%% perifocal frame coordinates
rNorm = hNorm^2 / muEarth / (1 + eccNorm * cos(theta));
rVecPQW = rNorm * [cos(theta), sin(theta), 0].';
vVecPQW = muEarth / hNorm * [-sin(theta), eccNorm + cos(theta), 0].';

%% change of basis
rotationMatrix = pvtRotationMatrix3D(3, OmegaRAAN) * pvtRotationMatrix3D(1, inc) * pvtRotationMatrix3D(3, omegaAP);

rVecECI = rotationMatrix * rVecPQW;
vVecECI = rotationMatrix * vVecPQW;

%% output
% already done

end

function rotationMatrix3D = pvtRotationMatrix3D(idAxis, angle)
% auxillary matrix to get the rotation matrix along axis-`idAxis` for an `angle`
% Inputs:
%   idAxis: 1, 2, 3
%   angle: [rad]
% Outputs:
%   rotationMatrix3D
%

% this definition of rotation matrix rotates a vector
rotationMatrix2D = [
    cos(angle), -sin(angle)
    sin(angle),  cos(angle)];

rotationMatrix3D = eye(3);
switch idAxis
    case 1
        tmpIds = [2, 3];
    case 2
        tmpIds = [3, 1];
    case 3
        tmpIds = [1, 2];
end
rotationMatrix3D(tmpIds, tmpIds) = rotationMatrix2D;
end











