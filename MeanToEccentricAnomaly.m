function [eccentricAnomaly, info] = MeanToEccentricAnomaly(ecc, meanAnomaly, typeVerbose, tol)
%% Convert mean anomlay to eccentric anomaly, by solving Kepler's equation using Newton's method.
%
% Inputs:
%   ecc: eccentricity in (0, 1)
%   meanAnomaly: [rad] mean anomaly
%   verbose: control displayment of iterations, 'none' (default), 'step', or 'final'.
%   tol: convergence tolerance (default: 1e-8 rad)
%
% Outputs:
%   eccentricAnomaly: eccentric anomaly [rad]
%   info: other information such as steps used, diff of the last step, etc.
%
% Author: Hao Peng (hao.peng@erau.edu)
% Versions:
%   - 202309281217, Initial development for AE 313 Fall 2023.

%% default values
if nargin < 4
    % if `tol` is not an input, we use the default value.
    tol = 1e-10;
    if nargin < 3
        % by default, we don't generate output
        typeVerbose = 'none';
        % typeVerbose = 'step';
        % typeVerbose = 'final';
    end
end
% max steps of iteration
iiStepMax = 2000000;
% only support eccentric orbits
if ecc < 0
    error('# HP: e < 0, wrong definition.')
elseif ecc == 0
    error('# HP: e = 0, no such conversion for circular orbits.')
elseif ecc == 1
    error('# HP: e = 1, no such conversion for parabolic orbits.')
end

%% we define two functions to simplify the code (handle both elliptic and hyperbolic orbits)

if ecc < 1 % elliptic
    % f(E) = E - e * sin(E) - Me
    funFx = @(E) E - ecc * sin(E) - meanAnomaly;
    % f'(E) = 1 - e * cos(E)
    funDFDx = @(E) 1 - ecc * cos(E);
else % hyperbolic
    % f(F) = e * sinh(F) - F - Mh
    funFx = @(F) ecc * sinh(F) - F - meanAnomaly;
    % f'(E) = e * cosh(F) - 1
    funDFDx = @(F) ecc * cosh(F) - 1;
end

%% the choice of the initial value matters a lot (handle both elliptic and hyperbolic orbits)
if ecc < 1 % elliptic
    
    % default in Curtis2020 Algorithm 3.1
    if meanAnomaly < pi
        eccentricAnomaly_0 = meanAnomaly + ecc / 2;
    else
        eccentricAnomaly_0 = meanAnomaly - ecc / 2;
    end

    % % a naive choice using just the input mean anomaly
    % eccentricAnomaly_0 = meanAnomaly; % this leads to more iterations
    % warning('# HP: naive initial values are used.')

    % % a even more naive choice of using just 0
    % eccentricAnomaly_0 = 0; % this leads to even more iterations, and even diverge.
    % warning('# HP: naive initial values are used.')

else % hyperbolic

    % default in Curtis2020 Algorithm 3.2
    eccentricAnomaly_0 = meanAnomaly;

end


%% the iteration loop
% prepare
iiStep = 0;
flagDone = false;

% remaining steps until converged
while ~flagDone
    % apply Newton's step
    eccentricAnomaly_1 = eccentricAnomaly_0 - funFx(eccentricAnomaly_0) / funDFDx(eccentricAnomaly_0);
    % test if converged already
    diff = abs(eccentricAnomaly_1 - eccentricAnomaly_0);
    if diff <= tol
        flagDone = true;
    end
    % verbose: print each iteration
    if strcmpi(typeVerbose, 'step')
        fprintf('# HP: step-%05d | E0 = %13.9f deg | E1 = %13.9f deg | diff = %13.9e deg \n', ...
            iiStep, rad2deg(eccentricAnomaly_0), rad2deg(eccentricAnomaly_1), rad2deg(eccentricAnomaly_1 - eccentricAnomaly_0));
    end
    % update values for the next step
    eccentricAnomaly_0 = eccentricAnomaly_1;
    % count iterations
    iiStep = iiStep + 1;
    % avoid infinite loop
    if iiStep > iiStepMax
        error('# HP: exceeds the maximum step of %d but does not converged with tol = %.2e.', iiStep, tol);
    end
end

% verbose: print the final step
if strcmpi(typeVerbose, 'final') || strcmpi(typeVerbose, 'step')
    fprintf("# ------\n");
    fprintf("# HP: Kepler's equation solved after [%d] steps\n", iiStep);
    fprintf("# HP:      eccentricity = %6.2f\n", ecc);
    fprintf("# HP:         tolerance = %10.2e\n", tol);
    fprintf("# HP:      mean anomaly = %13.9f rad\n", meanAnomaly);
    fprintf("# HP:                   = %13.9f deg\n", rad2deg(meanAnomaly));
    fprintf("# HP: eccentric anomaly = %13.9f rad\n", eccentricAnomaly_0);
    fprintf("# HP:                   = %13.9f deg\n", rad2deg(eccentricAnomaly_0));
end

%% the output
eccentricAnomaly = eccentricAnomaly_0;
if nargout == 2
    info.M = meanAnomaly;
    info.E = eccentricAnomaly;
    info.tol = tol;
    info.steps = iiStep;
    info.diffLastStep = diff; % [rad]
    info.errorKeplersEquation = funFx(eccentricAnomaly);
end

end