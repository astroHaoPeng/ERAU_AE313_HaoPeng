%% Solve Kepler's equation using Newton's method in textbook

% hyperbolic
ecc = 3.0;
tol = 1e-8;
meanAnomlayListDeg = -1e3 : 1 : 1e3; % This can to to infinity for hyperbolas.

%% solve
tmp = zeros(length(meanAnomlayListDeg), 3);

for ii = 1 : length(meanAnomlayListDeg)
    [~, info] = MeanToEccentricAnomaly(ecc, deg2rad(meanAnomlayListDeg(ii)), 'none', tol);
    tmp(ii, 1) = info.steps;
    tmp(ii, 2) = info.diffLastStep;
    tmp(ii, 3) = info.errorKeplersEquation;
    tmp(ii, 4) = info.F;
end

%% visualize results
lw = 1;
fs = 12;

figure(33);
set(gcf, 'Position',[10 10 560 800])

ax(1) = subplot(4, 1, 1);
plot(meanAnomlayListDeg, rad2deg(tmp(:,4)), 'b', 'LineWidth',lw);
ylabel('eccentric anomaly $E$ [deg]', 'Interpreter','latex', 'FontSize',fs)
text(60, -100, {['$e = ' num2str(ecc, '%.4f') '$'], ['tol = ' num2str(tol, '%.2e')]}, 'Interpreter','latex', 'FontSize',fs*1.2)
ylim(ax(1), [-360, 360]);
yticks(ax(1), [-360:120:360])

ax(2) = subplot(4, 1, 2);
plot(meanAnomlayListDeg, tmp(:,1), '-', 'LineWidth',lw);
ylabel('Steps used', 'Interpreter','latex', 'FontSize',fs)

ax(3) = subplot(4, 1, 3);
plot(meanAnomlayListDeg, rad2deg(tmp(:,2)), 'r', 'LineWidth',lw);
ylabel('last $|F_{i+1} - F_i|$ [deg]', 'Interpreter','latex', 'FontSize',fs)

ax(4) = subplot(4, 1, 4);
plot(meanAnomlayListDeg, tmp(:,3), 'b', 'LineWidth',lw);
ylabel('final $|e\sinh F - F - M_e|$', 'Interpreter','latex', 'FontSize',fs)
xlabel('input: mean anomaly $M_h$ [deg]', 'Interpreter','latex', 'FontSize',fs)
% axis touchup
for ii = 1:4
    axes(ax(ii))
    grid on;
    if ii > 1
        axis tight;
    end
end
