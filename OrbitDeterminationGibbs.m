function [oe2, v2Vec, oe1, v1Vec, oe3, v3Vec] = OrbitDeterminationGibbs(r1Vec, r2Vec, r3Vec)
%% Gibbs method for perliminary orbit determination using three ECI position vectors.
%   Inputs are not verified for the coplannar assumptionl.
%   Numerical stabilities are not considered.
%
% Inputs:
%   r1Vec, r2Vec, r3Vec: three position vectors [km]
%
% Outputs:
%   oe2: at the 2nd location, classical orbital elements consisting of 
%        [h, e, theta, RAAN, i, AP] in [km^2/s, 1, rad, rad, rad, rad]
%   v2Vec: [km/s] at the 2nd location, Gibbs approximation of the velocity at the center point
%   oe1, v1Vec: similar but at the 1st location
%   oe3, v3Vec: similar but at the 3rd location
%
% Author: Hao Peng (hao.peng@erau.edu)
% Updates: Initially used as a class demo of AE 313 in Fall 2023.
%          20240317, updated with more output for comparison among using different velocities


%% constants
muEarth = 398600;         % gravitational parameter [km^3/s^2]

r1Norm = sum(r1Vec.^2)^(1/2);
r2Norm = sum(r2Vec.^2)^(1/2);
r3Norm = sum(r3Vec.^2)^(1/2);

D = cross(r1Vec, r2Vec) + cross(r2Vec, r3Vec) + cross(r3Vec, r1Vec);
N = r1Norm * cross(r2Vec, r3Vec) + r2Norm * cross(r3Vec, r1Vec) + r3Norm * cross(r1Vec, r2Vec);
S = (r2Norm - r3Norm) * r1Vec + (r3Norm - r1Norm) * r2Vec + (r1Norm - r2Norm) * r3Vec;

% test if the three given vectors are coplannar enough
tmp = norm(cross(D/norm(D), N/norm(N)));
if tmp > 1e-2
    warning('# Hao: D/|D| X N/|N|) = %.5e > 1e-2, not quite coplannar, so Gibbs method may not function properly!', tmp)
end

DNorm = sum(D.^2)^(1/2);
NNorm = sum(N.^2)^(1/2);

v2Vec = sqrt(muEarth/NNorm/DNorm) * (cross(D, r2Vec)/r2Norm + S);

oe2 = ConvertRvToCoe(r2Vec, v2Vec);

% only generate additional outputs when required for
if nargout > 2
    v1Vec = sqrt(muEarth/NNorm/DNorm) * (cross(D, r1Vec)/r1Norm + S);
    oe1 = ConvertRvToCoe(r1Vec, v1Vec);
    v3Vec = sqrt(muEarth/NNorm/DNorm) * (cross(D, r3Vec)/r3Norm + S);
    oe3 = ConvertRvToCoe(r3Vec, v3Vec);
    % some codes for debugging
    % num2str([oe1; oe2; oe3]) % show values
    % num2str(diff([oe1; oe2; oe3], 1, 1)) % difference of oe at three places, should be similar except for true anomaly
end

end