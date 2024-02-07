function trueAnomaly = TimeToTrueAnomalyHyperbolic(ecc, time, h, muEarth, typeVerbose, tol)

%% default values
if nargin < 6
    % if `tol` is not an input, we use the default value.
    tol = 1e-10;
    if nargin < 5
        % by default, we don't generate output
        typeVerbose = 'none';
        % typeVerbose = 'step';
        % typeVerbose = 'final';
    end
end

%% t --> Mh
meanAnomaly = muEarth^2 / h^3 * (ecc^2 - 1)^(3/2) * time;

%% Mh --> F
eccentricAnomaly = zeros(size(time));
for ii = 1:length(time)
    eccentricAnomaly(ii) = MeanToEccentricAnomaly(ecc, meanAnomaly(ii), typeVerbose, tol);
end

%% E --> theta
trueAnomaly = 2 * atan( sqrt((1+ecc)/(-1+ecc)) .* tanh(eccentricAnomaly/2) );
trueAnomaly = mod(trueAnomaly, 2*pi); % bring to [0, 2*pi) range

end