function [oe, v2Vec] = OrbitDeterminationGibbs(r1Vec, r2Vec, r3Vec)
%% Gibbs method for perliminary orbit determination using three ECI position vectors.
%   Inputs are not verified for the coplannar assumptionl.
%   Numerical stabilities are not considered.
%
%   Initially used as a class demo of AE 313 in Fall 2023.
%
% Inputs:
%   r1Vec, r2Vec, r3Vec: three position vectors [km]
%
% Outputs:
%   oe: classical orbital elements consisting of [h, e, theta, RAAN, i, AP]
%   v2Vec: [km/s] Gibbs approximation of the velocity at the center point
%
% Author: Hao Peng (hao.peng@erau.edu)

%% constants
muEarth = 398600;         % gravitational parameter [km^3/s^2]

r1Norm = sum(r1Vec.^2)^(1/2);
r2Norm = sum(r2Vec.^2)^(1/2);
r3Norm = sum(r3Vec.^2)^(1/2);

D = cross(r1Vec, r2Vec) + cross(r2Vec, r3Vec) + cross(r3Vec, r1Vec);
N = r1Norm * cross(r2Vec, r3Vec) + r2Norm * cross(r3Vec, r1Vec) + r3Norm * cross(r1Vec, r2Vec);
S = (r2Norm - r3Norm) * r1Vec + (r3Norm - r1Norm) * r2Vec + (r1Norm - r2Norm) * r3Vec;
if norm(cross(D, N)) > 1e-5
    warning('# Hao: dot(D, N) = %.5e > 1e-5, Gibbs method may not function properly!', norm(cross(D, N)))
end

DNorm = sum(D.^2)^(1/2);
NNorm = sum(N.^2)^(1/2);

v2Vec = sqrt(muEarth/NNorm/DNorm) * (cross(D, r2Vec)/r2Norm + S);

oe = ConvertRvToCoe(r2Vec, v2Vec);
end