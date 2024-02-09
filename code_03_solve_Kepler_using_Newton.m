%% Solve Kepler's equation using Newton's method in textbook

%% single point solution
ecc = 0.95;
tol = 1e-8;
meanAnomalyDeg = 120; % [deg]
eccentricAnomaly = MeanToEccentricAnomaly(ecc, deg2rad(meanAnomalyDeg), 'step', tol); % [rad]


%% solution for the entire domain [0, 2*pi]
ecc = 0.1;
tol = 1e-8;
meanAnomalyListDeg = 0 : 0.1 : 360;

% solve for each point
eccentricAnomalyList = zeros(length(meanAnomalyListDeg), 1);
tmp = zeros(length(meanAnomalyListDeg), 3); % initialize the storage space
for ii = 1 : length(meanAnomalyListDeg)
    [eccentricAnomaly, info] = MeanToEccentricAnomaly(ecc, deg2rad(meanAnomalyListDeg(ii)), 'none', tol);
    tmp(ii, 1) = info.steps;
    tmp(ii, 2) = info.diffLastStep;
    tmp(ii, 3) = info.errorKeplersEquation;
    if ecc < 1
        eccentricAnomalyList(ii) = eccentricAnomaly;
    elseif ecc > 1
        eccentricAnomalyList(ii) = eccentricAnomaly;
    end
end

% visualize results
lw = 1;
fs = 12;

figure(33);
clf;
set(gcf, 'Position',[10 10 560 800])

ax(1) = subplot(4, 1, 1);
plot([0 360], [0 360], 'k:', 'LineWidth',lw); hold on
plot(meanAnomalyListDeg, rad2deg(eccentricAnomalyList), 'b', 'LineWidth',lw);
ylabel('eccentric anomaly $E$ [deg]', 'Interpreter','latex', 'FontSize',fs)
text(60, 250, ['$e = ' num2str(ecc, '%.4f') '$, tol = ' num2str(tol, '%.2e')], 'Interpreter','latex', 'FontSize',fs*1.2)

ax(2) = subplot(4, 1, 2);
plot(meanAnomalyListDeg, tmp(:,1), '-', 'LineWidth',lw);
ylabel('Steps used', 'Interpreter','latex', 'FontSize',fs)

ax(3) = subplot(4, 1, 3);
plot(meanAnomalyListDeg, rad2deg(tmp(:,2)), 'r', 'LineWidth',lw);
ylabel('last $|E_{i+1} - E_i|$ [deg]', 'Interpreter','latex', 'FontSize',fs)

ax(4) = subplot(4, 1, 4);
plot(meanAnomalyListDeg, tmp(:,3), 'b', 'LineWidth',lw);
ylabel('final $|E - e\sin E - M_e|$', 'Interpreter','latex', 'FontSize',fs)
xlabel('input: mean anomaly $M_e$ [deg]', 'Interpreter','latex', 'FontSize',fs)
% axis touchup
for ii = 1:4
    axes(ax(ii))
    grid on;
    axis tight;
end