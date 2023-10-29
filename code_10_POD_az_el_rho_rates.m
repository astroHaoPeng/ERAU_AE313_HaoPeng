%% Az, El, rho, and their rates
clear all; clc;

tic;

% constants
muEarth = 398600; % [km^3/s^2]
rEarth = 6378; % [km]
omegaEarth = deg2rad(360) / 86164.0905; % Earth rotation rate [rad/s] (note: must using sidereal time of a day)

% radar station site parameters
thetaLst = deg2rad(240.7); % local sidereal time [rad]
latitude = deg2rad(42.9); % latitude [rad]
height = 379.476 / 1e3; % altitude [km] % use `height` to avoid potentail typo or ambiguity with `latitude`

% measurements from a radar station
az = deg2rad(135.4); % azimuth [rad]
el = deg2rad(62.5); % elevation [rad]
rho = 668.3; % range [km]
azDot = deg2rad(-0.38); % azimuth rate [rad/s]
elDot = deg2rad(-0.65); % elevation rate [rad/s]
rhoDot = 2.39; % range rate [km/s]

% site position vector, ECI coordinates.
rSite = height + rEarth;
rSiteVecECI = rSite * [cos(latitude)*cos(thetaLst), cos(latitude)*sin(thetaLst), sin(latitude)].';

% Earth rotation vector, ECI coordinates.
omegaEarthVecECI = [0, 0, omegaEarth].';

% rotation matrix that rotates {IJK} to {ijk} by two Euler angle rotations
% ECI --rotated-to--> topocentric ENZ
Phi_3 = [
    cos(thetaLst+pi/2) -sin(thetaLst+pi/2) 0;
    sin(thetaLst+pi/2)  cos(thetaLst+pi/2) 0;
    0                   0                  1];
Phi_1 = [
    1 0              0;
    0 cos(pi/2-latitude) -sin(pi/2-latitude);
    0 sin(pi/2-latitude)  cos(pi/2-latitude)];
Phi = Phi_3 * Phi_1;

% 
unitRhoVecENZ = [cos(el)*sin(az), cos(el)*cos(az), sin(el)].';
unitRhoVecECI = Phi * unitRhoVecENZ;
unitRhoVecDotECI = Phi * [- sin(el)*elDot*sin(az) + cos(el)*cos(az)*azDot,...
                          - sin(el)*elDot*cos(az) - cos(el)*sin(az)*azDot,...
                            cos(el)*elDot].' ...
                   + cross(omegaEarthVecECI, unitRhoVecECI);
rECI = rSiteVecECI + rho * unitRhoVecECI;
vECI = cross(omegaEarthVecECI, rSiteVecECI) ...
        + rhoDot * unitRhoVecECI ...
        + rho * unitRhoVecDotECI;

% display results
disp(['# Hao: ECI position = [' num2str(rECI.', '%+15.6e') '] km'])
disp(['# Hao: ECI velocity = [' num2str(vECI.', '%+15.6e') '] km/s'])
toc;

% ECI RV -> OE
oe = ConvertRvToCoe(rECI, vECI)



