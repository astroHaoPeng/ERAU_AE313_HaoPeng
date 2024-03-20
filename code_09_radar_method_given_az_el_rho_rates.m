%% Radar method
%   Az, El, rho, and their rates
%   These measurements are given in the local horizon frame, aka, topocentric frame,
%   or E(ast)-N(orth)-Z(enith) frame.

clear all; clc;

tic;

%% constants
muEarth = 398600; % [km^3/s^2]
rEarth = 6378; % [km]
omegaEarth = deg2rad(360) / 86164.0905; % Earth rotation rate [rad/s] (note: must using sidereal time of a day)

%% radar station site parameters and measurements from the radar station
% one example
thetaLst = deg2rad(240.7); % local sidereal time [rad]
latitude = deg2rad(42.9); % latitude [rad]
height = 379.476 / 1e3; % altitude [km] % use `height` to avoid potentail typo or ambiguity with `latitude`
az = deg2rad(135.4); % azimuth [rad]
el = deg2rad(62.5); % elevation [rad]
rho = 668.3; % range [km]
azDot = deg2rad(-0.38); % azimuth rate [rad/s]
elDot = deg2rad(-0.65); % elevation rate [rad/s]
rhoDot = 2.39; % range rate [km/s]

% another example in textbook: EXAMPLE 5.10
% thetaLst = deg2rad(300); % local sidereal time [rad]
% latitude = deg2rad(60); % latitude [rad]
% height = 0 / 1e3; % altitude [km] % use `height` to avoid potentail typo or ambiguity with `latitude`
% az = deg2rad(90); % azimuth [rad]
% el = deg2rad(30); % elevation [rad]
% rho =  2551; % range [km]
% azDot = deg2rad(0.1130); % azimuth rate [rad/s]
% elDot = deg2rad(0.05651); % elevation rate [rad/s]
% rhoDot = 0; % range rate [km/s]


%% site position vector, ECI coordinates.
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

Psi = [
    1 0              0;
    0 cos(-(pi/2-latitude)) -sin(-(pi/2-latitude));
    0 sin(-(pi/2-latitude))  cos(-(pi/2-latitude))] * ...
    [
    cos(-(thetaLst+pi/2)) -sin(-(thetaLst+pi/2)) 0;
    sin(-(thetaLst+pi/2))  cos(-(thetaLst+pi/2)) 0;
    0                   0                  1];

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
disp( '# Peng: Results from given radar measurements')
disp(['#       ECI position = [' num2str(rECI.', '%+15.6e') '] km'])
disp(['#       ECI velocity = [' num2str(vECI.', '%+15.6e') '] km/s'])

% ECI RV -> OE
oe = ConvertRvToCoe(rECI, vECI);

PrintResults(thetaLst, oe)

% timing
toc;

%%
function PrintResults(lst_, oe_)
fprintf('#----------------------------------------------------------\n')
fprintf('# At local sidereal time (LST) = %.3f [deg] \n', rad2deg(lst_))
fprintf('#            h = %.10f [km^2/s]\n', oe_(1))
fprintf('#            e = %.10f \n', oe_(2))
fprintf('#        theta = %.10f [rad] = %.10f [deg]\n', oe_(3), rad2deg(oe_(3)))
fprintf('#         RAAN = %.10f [rad] = %.10f [deg]\n', oe_(4), rad2deg(oe_(4)))
fprintf('#  inclination = %.10f [rad] = %.10f [deg]\n', oe_(5), rad2deg(oe_(5)))
fprintf('#           AP = %.10f [rad] = %.10f [deg]\n', oe_(6), rad2deg(oe_(6)))
end


