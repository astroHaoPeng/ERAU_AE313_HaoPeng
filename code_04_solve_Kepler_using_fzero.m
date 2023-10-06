%% Solve Kepler's equation using MATALB built-in function fzero
clear; clc;

%% parameters
ecc = 10.95;
tol = 1e-8;

%% compare results
meanAnomlayList = deg2rad(0 : 0.1 : 360);

tmpNewton = zeros(1, length(meanAnomlayList));
tmpFzero = zeros(1, length(meanAnomlayList));

options = optimset('TolX', tol);

for ii = 1:length(meanAnomlayList)
    meanAnomaly = meanAnomlayList(ii);
    if ecc < 1
        tmpFzero(ii) = fzero( @(E) E - ecc * sin(E) - meanAnomaly, meanAnomaly, options); % elliptic
    elseif ecc > 1
        tmpFzero(ii) = fzero( @(F) ecc * sinh(F) - F - meanAnomaly, meanAnomaly, options); % hyperbolic
    end
    tmpNewton(ii) = MeanToEccentricAnomaly(ecc, meanAnomaly, 'none', tol);
end

% %% visualize
figure(43);
clf;

subplot (2,1,1)
plot(meanAnomlayList, rad2deg(tmpFzero), 'r', 'DisplayName','fzero');
hold on;
plot(meanAnomlayList, rad2deg(tmpNewton), 'b', 'DisplayName',"Newton's method");
l = legend('Location','best');
axis tight;
grid on; 
grid minor;

subplot (2,1,2)
plot(meanAnomlayList, rad2deg(tmpFzero - tmpNewton), 'k');
axis tight;
title("differecen between `fzero` and `Newton's method` solutions");

%% compare speed of `fzero` and our `MeanToEccentricAnomaly`
tic;
for ii = 1:length(meanAnomlayList)
    meanAnomaly = meanAnomlayList(ii);
    if ecc < 1
        tmpFzero(ii) = fzero( @(E) E - ecc * sin(E) - meanAnomaly, meanAnomaly, options); % elliptic
    elseif ecc > 1
        tmpFzero(ii) = fzero( @(F) ecc * sinh(F) - F - meanAnomaly, meanAnomaly, options); % hyperbolic
    end
end
toc;

tic;
for ii = 1:length(meanAnomlayList)
    meanAnomaly = meanAnomlayList(ii);
    tmpNewton(ii) = MeanToEccentricAnomaly(ecc, meanAnomaly, 'none', tol);
end
toc;