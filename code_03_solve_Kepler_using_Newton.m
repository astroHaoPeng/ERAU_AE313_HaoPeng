%% Solve Kepler's equation using Newton's method in textbook

% elliptic
ecc = 0.99;
tol = 1e-8;
meanAnomlayListDeg = 0 : 0.1 : 360;

% % hyperbolic
% ecc = 5;
% tol = 1e-8;
% meanAnomlayListDeg = -1e4 : 1 : 1e4;

%% solve
tmp = zeros(length(meanAnomlayListDeg), 3);

for ii = 1 : length(meanAnomlayListDeg)
    [~, info] = MeanToEccentricAnomaly(ecc, deg2rad(meanAnomlayListDeg(ii)), 'none', tol);
    tmp(ii, 1) = info.steps;
    tmp(ii, 2) = info.diffLastStep;
    tmp(ii, 3) = info.errorKeplersEquation;
    tmp(ii, 4) = info.E;
end

%% visualize results
lw = 1;
fs = 12;

figure(33);
set(gcf, 'Position',[488 726 560 800])

ax(1) = subplot(4, 1, 1);
plot(meanAnomlayListDeg, rad2deg(tmp(:,4)), 'b', 'LineWidth',lw);
ylabel('eccentric anomaly $E$ [deg]', 'Interpreter','latex', 'FontSize',fs)
text(60, 250, ['$e = ' num2str(ecc, '%.4f') '$, tol = ' num2str(tol, '%.2e')], 'Interpreter','latex', 'FontSize',fs*1.2)

ax(2) = subplot(4, 1, 2);
plot(meanAnomlayListDeg, tmp(:,1), '-', 'LineWidth',lw);
ylabel('Steps used', 'Interpreter','latex', 'FontSize',fs)

ax(3) = subplot(4, 1, 3);
plot(meanAnomlayListDeg, rad2deg(tmp(:,2)), 'r', 'LineWidth',lw);
ylabel('last $|E_{i+1} - E_i|$ [deg]', 'Interpreter','latex', 'FontSize',fs)

ax(4) = subplot(4, 1, 4);
plot(meanAnomlayListDeg, tmp(:,3), 'b', 'LineWidth',lw);
ylabel('final $|E - e\sin E - M_e|$', 'Interpreter','latex', 'FontSize',fs)
xlabel('input: mean anomaly $M_e$ [deg]', 'Interpreter','latex', 'FontSize',fs)
% axis touchup
for ii = 1:4
    axes(ax(ii))
    grid on;
    axis tight;
end