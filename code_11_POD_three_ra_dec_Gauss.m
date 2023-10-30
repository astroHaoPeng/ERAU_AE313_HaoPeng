%% Gauss method
%
% no detailed comments are added to the script to encourage reviewing the lecture slide
%

clear all; clc;

tic;

% constants
muEarth = 398600; % [km^3/s^2]
rEarth = 6378; % [km]
omegaEarth = deg2rad(360) / 86164.0905; % Earth rotation rate [rad/s] (note: must using sidereal time of a day)

% radar station site parameters
latitude = deg2rad(40.0); % latitude [rad]
height = 1.0; % altitude [km] % use `height` to avoid potentail typo or ambiguity with `latitude`

% measurements: three sets of (t, ra, dec)
time = [0 118.10 237.58]; % [s]
thetaLst = deg2rad([44.506 45.000 45.499]); % local sidereal time [rad]
rightAscension = deg2rad([43.537 54.420 64.318]); % right ascension [rad]
declination = deg2rad([-8.7833 -12.074 -15.105]); % declination [rad]

% 
tau = time(3) - time(1);
tau1 = time(1) - time(2);
tau3 = time(3) - time(2);

unitRhoVecTopo1 = [cos(declination(1))*cos(rightAscension(1)), cos(declination(1))*sin(rightAscension(1)), sin(declination(1))].';
unitRhoVecTopo2 = [cos(declination(2))*cos(rightAscension(2)), cos(declination(2))*sin(rightAscension(2)), sin(declination(2))].';
unitRhoVecTopo3 = [cos(declination(3))*cos(rightAscension(3)), cos(declination(3))*sin(rightAscension(3)), sin(declination(3))].';

matrixP = [unitRhoVecTopo1, unitRhoVecTopo2, unitRhoVecTopo3];

sitePositionVecECI1 = (rEarth + height) * [cos(latitude)*cos(thetaLst(1)), cos(latitude)*sin(thetaLst(1)), sin(latitude)].';
sitePositionVecECI2 = (rEarth + height) * [cos(latitude)*cos(thetaLst(2)), cos(latitude)*sin(thetaLst(2)), sin(latitude)].';
sitePositionVecECI3 = (rEarth + height) * [cos(latitude)*cos(thetaLst(3)), cos(latitude)*sin(thetaLst(3)), sin(latitude)].';

% debugging data, copied from Curtis2020 example 5.11
% sitePositionVecECI1 = [3489.8, 3430.2, 4078.5].';
% sitePositionVecECI2 = [3460.1, 3460.1, 4078.5].';
% sitePositionVecECI3 = [3429.9, 3490.1, 4078.5].';


matrixR = [sitePositionVecECI1, sitePositionVecECI2, sitePositionVecECI3];

M = matrixP \ matrixR; %inv(matrixP) * matrixR;

A = M(2, 1) * tau3/tau - M(2, 2) - M(2, 3) * tau1/tau;
B = M(2, 1) * tau3/6/tau * (tau^2 - tau3^2) - M(2, 3) * tau1/6/tau * (tau^2 - tau1^2);

a = - ( A^2 + 2*A*dot(sitePositionVecECI2, unitRhoVecTopo2) + sum(sitePositionVecECI2.^2) );
b = - 2 * muEarth * B * ( A + dot(sitePositionVecECI2, unitRhoVecTopo2) );
c = - muEarth^2 * B^2;

potential_r2s = roots([1, 0, a, 0, 0, b, 0, 0, c]);
real_r2s = potential_r2s(imag(potential_r2s) == 0 & real(potential_r2s) > 0);

if length(real_r2s) > 1
    fprintf('# Hao: all the real roots from the 8-degree polynominal equation of `r2` is:\n')
    for ii = 1 : length(real_r2s)
        fprintf('# Hao:   %1d : %.2f\n', ii, real_r2s(ii))
    end
    id_r2_input = input('# Hao: please enter the id of desired `r2`, then press `Enter` to continue: ');
    r2 = real_r2s(id_r2_input);
else
    r2 = real_r2s;
end

c1 = tau3/tau * ( 1 + 1/6 * muEarth / r2^3 * (tau^2 - tau3^2) );
c3 = - tau1/tau * ( 1 + 1/6 * muEarth / r2^3 * (tau^2 - tau1^2) );

tmpVec = M * [-c1, 1, -c3].';
rho1 = tmpVec(1) / c1;
rho2 = tmpVec(2) / (-1);
rho3 = tmpVec(3) / c3;

r1VecECI = sitePositionVecECI1 + rho1 * unitRhoVecTopo1;
r2VecECI = sitePositionVecECI2 + rho2 * unitRhoVecTopo2;
r3VecECI = sitePositionVecECI3 + rho3 * unitRhoVecTopo3;

f1 = 1 - muEarth /2/r2^3 * tau1^2;
g1 = tau1 - muEarth /6/r2^3 * tau1^3;
f3 = 1 - muEarth /2/r2^3 * tau3^2;
g3 = tau3 - muEarth /6/r2^3 * tau3^3;

v2VecECI_from_r1Vec = (r1VecECI - f1 * r2VecECI) / g1;
v2VecECI_from_r3Vec = (r3VecECI - f3 * r2VecECI) / g3;
v2VecECI_from_both = (f3*r1VecECI - f1*r3VecECI) / (f3*g1 - f1*g3);


% display preliminary results
disp(['# Hao: ECI position            = [' num2str(r2VecECI.', '%+15.6e') '] km'])
disp(['# Hao: ECI velocity from r1Vec = [' num2str(v2VecECI_from_r1Vec.', '%+15.6e') '] km/s'])
disp(['# Hao: ECI velocity from r3Vec = [' num2str(v2VecECI_from_r3Vec.', '%+15.6e') '] km/s'])
disp(['# Hao: ECI velocity from both  = [' num2str(v2VecECI_from_both.', '%+15.6e') '] km/s'])
toc;



%% an example of iteration
typeIteration = 'default';
% typeIteration = 'Gibbs';
if strcmpi(typeIteration, 'Gibbs')
    oe2 = OrbitDeterminationGibbs(r1VecECI, r2VecECI, r3VecECI);
    oe2Iter = oe2;
else
    oe2 = ConvertRvToCoe(r2VecECI, v2VecECI_from_both);
    oe2Iter = oe2;
end
c1Iter = c1;
c3Iter = c3;
r1VecECIIter = r1VecECI;
r2VecECIIter = r2VecECI;
r3VecECIIter = r3VecECI;
v2VecECIIter = v2VecECI_from_both;
kk = 0;
while true
    % calculate true anomalies at t1 and t2
    hIter = oe2Iter(1);
    eIter = oe2Iter(2);
    if eIter < 1 % elliptical orbits
        aIter = hIter^2 / muEarth / (1-eIter^2);
        meanMotion = sqrt(muEarth) / aIter^(3/2);
        trueAnomaly2 = oe2Iter(3);
        eccentricAnomaly2 = 2 * atan( sqrt((1-eIter)/(1+eIter)) * tan(trueAnomaly2/2) );
        meanAnomaly2 = eccentricAnomaly2 - eIter * sin(eccentricAnomaly2);
        timeSincePerigee2 = meanAnomaly2 / meanMotion;
    else % hyperbolic orbits
        trueAnomaly2 = oe2Iter(3);
        eccentricAnomaly2 = 2 * atanh( sqrt((eIter-1)/(eIter+1)) * tan(trueAnomaly2/2) );
        meanAnomaly2 = eIter * sinh(eccentricAnomaly2) - eccentricAnomaly2;
        timeSincePerigee2 = meanAnomaly2 / (muEarth^2/hIter^3*(eIter^2-1)^(3/2));
    end
    %
    timeSincePerigee1 = timeSincePerigee2 + tau1;
    if eIter < 1
        meanAnomaly1 = meanMotion * timeSincePerigee1;
        eccentricAnomaly1 = MeanToEccentricAnomaly(eIter, meanAnomaly1);
        trueAnomaly1 = 2 * atan( sqrt((1+eIter)/(1-eIter)) * tan(eccentricAnomaly1/2) );
    else
        meanAnomaly1 = (muEarth^2/hIter^3*(eIter^2-1)^(3/2)) * timeSincePerigee1;
        eccentricAnomaly1 = MeanToEccentricAnomaly(eIter, meanAnomaly1);
        trueAnomaly1 = 2 * atan( sqrt((eIter+1)/(eIter-1)) * tanh(eccentricAnomaly1/2) );
    end
    oe1Iter = [oe2Iter(1:2), trueAnomaly1, oe2Iter(4:6)];
    %
    timeSincePerigee3 = timeSincePerigee2 + tau3;
    if eIter < 1
        meanAnomaly3 = meanMotion * timeSincePerigee3;
        eccentricAnomaly3 = MeanToEccentricAnomaly(eIter, meanAnomaly3);
        trueAnomaly3 = 2 * atan( sqrt((1+eIter)/(1-eIter)) * tan(eccentricAnomaly3/2) );
    else
        meanAnomaly3 = (muEarth^2/hIter^3*(eIter^2-1)^(3/2)) * timeSincePerigee3;
        eccentricAnomaly3 = MeanToEccentricAnomaly(eIter, meanAnomaly3);
        trueAnomaly3 = 2 * atan( sqrt((eIter+1)/(eIter-1)) * tanh(eccentricAnomaly3/2) );
    end
    oe3Iter = [oe2Iter(1:2), trueAnomaly3, oe2Iter(4:6)];

    % calculate f1, g1, f3, g3 using expressions of deltaTheta
    deltaTheta1 = trueAnomaly1 - trueAnomaly2;
    [tmp_r1VecECI, tmp_v1VecECI] = ConvertCoeToRv(oe1Iter);
    [tmp_r2VecECI, tmp_v2VecECI] = ConvertCoeToRv(oe2Iter);
    [tmp_r3VecECI, tmp_v3VecECI] = ConvertCoeToRv(oe3Iter);
    f1Iter = 1 - muEarth * norm(tmp_r1VecECI) / hIter^2 * (1 - cos(deltaTheta1));
    g1Iter = norm(tmp_r1VecECI) * norm(tmp_r2VecECI) / hIter * sin(deltaTheta1);
    deltaTheta3 = trueAnomaly3 - trueAnomaly2;
    f3Iter = 1 - muEarth * norm(tmp_r3VecECI) / hIter^2 * (1 - cos(deltaTheta3));
    g3Iter = norm(tmp_r3VecECI) * norm(tmp_r2VecECI) / hIter * sin(deltaTheta3);

    % calculate c1, c3
    tmpDenominator = f1Iter*g3Iter - f3Iter*g1Iter;
    c1IterOld = c1Iter;
    c3IterOld = c3Iter;
    c1Iter =   g3Iter / tmpDenominator;
    c3Iter = - g1Iter / tmpDenominator;
    c1Iter = (c1IterOld + c1Iter) / 2;
    c3Iter = (c3IterOld + c3Iter) / 2;

    % calculate rho1, rho2, rho3 and others
    tmpVec = M * [-c1Iter, 1, -c3Iter].';
    rho1Iter = tmpVec(1) / c1Iter;
    rho2Iter = tmpVec(2) / (-1);
    rho3Iter = tmpVec(3) / c3Iter;

    r1VecECIIterOld = r1VecECIIter;
    r2VecECIIterOld = r2VecECIIter;
    r3VecECIIterOld = r3VecECIIter;
    r1VecECIIter = sitePositionVecECI1 + rho1Iter * unitRhoVecTopo1;
    r2VecECIIter = sitePositionVecECI2 + rho2Iter * unitRhoVecTopo2;
    r3VecECIIter = sitePositionVecECI3 + rho3Iter * unitRhoVecTopo3;

    if strcmpi(typeIteration, 'Gibbs')
        [oe2IterNext, v2VecECIIter] = OrbitDeterminationGibbs(r1VecECIIter, r2VecECIIter, r3VecECIIter);
    else
        v2VecECIIter = (f3Iter*r1VecECIIter - f1Iter*r3VecECIIter) / (f3Iter*g1Iter - f1Iter*g3Iter);
        oe2IterNext = ConvertRvToCoe(r2VecECIIter, v2VecECIIter);
    end

    % disp('# Hao: true anomalies:')
    % disp([trueAnomaly1, trueAnomaly2, trueAnomaly3])
    fprintf('# Hao: ~~~~~~~~~~ Iter-%d ~~~~~~~~~~\n', kk)
    fprintf('# Hao: rho comparison (km):\n')
    fprintf('# Hao:   before   %15.6e, %15.6e, %15.6e\n', rho1, rho2, rho3);
    fprintf('# Hao:   now      %15.6e, %15.6e, %15.6e\n', rho1Iter, rho2Iter, rho3Iter);
    fprintf('# Hao: v2Vec comparison  (km/s):\n')
    fprintf('# Hao:   before   %15.6e, %15.6e, %15.6e\n', v2VecECI_from_both(1), v2VecECI_from_both(2), v2VecECI_from_both(3));
    fprintf('# Hao:   now      %15.6e, %15.6e, %15.6e\n', v2VecECIIter(1), v2VecECIIter(2), v2VecECIIter(3));
    fprintf('# Hao: OE comparison  (km^2/s^2, 1, deg, deg, deg, deg):\n')
    tmp = [oe2(1:2), rad2deg(oe2(3:6))];
    fprintf('# Hao:   before   %13.6e, %8.5f, %7.2f, %7.2f, %7.2f, %7.2f\n',...
        tmp(1), tmp(2), tmp(3), tmp(4), tmp(5), tmp(6));
    tmp = [oe2Iter(1:2), rad2deg(oe2Iter(3:6))];
    fprintf('# Hao:   now      %13.6e, %8.5f, %7.2f, %7.2f, %7.2f, %7.2f\n',...
        tmp(1), tmp(2), tmp(3), tmp(4), tmp(5), tmp(6));
    
    if max(abs(oe2Iter - oe2IterNext)) < 1
        disp('# Hao: Iteration correction done!')
        % display refined results
        disp(['# ------------ After ' num2str(kk) ' iterations -------------------------'])
        disp(['# Hao: ECI position         = [' num2str(r2VecECI.', '%+15.6e') '] km'])
        disp(['# Hao: ECI position refined = [' num2str(r2VecECIIter.', '%+15.6e') '] km/s'])
        disp(['# Hao: ECI velocity         = [' num2str(v2VecECI_from_both.', '%+15.6e') '] km/s'])
        disp(['# Hao: ECI velocity refined = [' num2str(v2VecECIIter.', '%+15.6e') '] km/s'])
        toc;
        break
    else
        oe2Iter = oe2IterNext;
        kk = kk + 1;
    end

end