%% Solve Kepler's equation using MATALB built-in function fzero

%% parameters
ecc = 0.5;
tol = 1e-8;

%% compare results
meanAnomlayList = deg2rad(0 : 0.1 : 360);

tmpNewton = zeros(1, length(meanAnomlayList));
tmpFzero = zeros(1, length(meanAnomlayList));

options = optimset('TolX', tol);

for ii = 1:length(meanAnomlayList)
    meanAnomaly = meanAnomlayList(ii);
    tmpFzero(ii) = fzero( @(E) E - ecc * sin(E) - meanAnomaly, meanAnomaly, options);
    tmpNewton(ii) = MeanToEccentricAnomaly(ecc, meanAnomaly, 'none', tol);
end

% visualize
figure(33);
clf;

subplot (2,1,1)
plot(meanAnomlayList, tmpFzero, 'r', 'DisplayName','fzero');
hold on;
plot(meanAnomlayList, tmpNewton, 'b', 'DisplayName',"Newton's method");
l = legend('Location','best');
axis tight;
grid on; 
grid minor;

subplot (2,1,2)
plot(meanAnomlayList, tmpFzero - tmpNewton, 'k');
axis tight;

%% compare speed
tic;
for ii = 1:length(meanAnomlayList)
    meanAnomaly = meanAnomlayList(ii);
    tmpFzero(ii) = fzero( @(E) E - ecc * sin(E) - meanAnomaly, meanAnomaly, options);
end
toc;

tic;
for ii = 1:length(meanAnomlayList)
    meanAnomaly = meanAnomlayList(ii);
    tmpNewton(ii) = MeanToEccentricAnomaly(ecc, meanAnomaly, 'none', tol);
end
toc;