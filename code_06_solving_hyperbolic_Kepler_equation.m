%% Solve Kepler's equation using Newton's method in textbook

% hyperbolic
ecc = 3.0;
tol = 1e-8;
meanAnomlayListRad = -1e2 : 1 : 1e2; % This can to to infinity for hyperbolas.

%% solve for the entire range
stepsList = zeros(length(meanAnomlayListRad), 1);
diffLastStepList = zeros(length(meanAnomlayListRad), 1);
errorKeplersEquationList = zeros(length(meanAnomlayListRad), 1);
eccentricAnomalyListRad = zeros(length(meanAnomlayListRad), 1);

for ii = 1 : length(meanAnomlayListRad)
    [~, info] = MeanToEccentricAnomaly(ecc, meanAnomlayListRad(ii), 'none', tol);
    stepsList(ii) = info.steps;
    diffLastStepList(ii) = info.diffLastStep;
    errorKeplersEquationList(ii) = info.errorKeplersEquation;
    eccentricAnomalyListRad(ii) = info.F;
end

%% visualize results
lw = 1;
fs = 12;

figure(33);
set(gcf, 'Position',[10 10 560 800])

ax(1) = subplot(4, 1, 1);
plot(meanAnomlayListRad, eccentricAnomalyListRad, 'b', 'LineWidth',lw);
ylabel('eccentric anomaly $E$ [rad]', 'Interpreter','latex', 'FontSize',fs)
text(60, -100, {['$e = ' num2str(ecc, '%.4f') '$'], ['tol = ' num2str(tol, '%.2e')]}, 'Interpreter','latex', 'FontSize',fs*1.2)
% ylim(ax(1), [-360, 360]);
% yticks(ax(1), [-360:120:360])

ax(2) = subplot(4, 1, 2);
plot(meanAnomlayListRad, stepsList, '-', 'LineWidth',lw);
ylabel('Steps used', 'Interpreter','latex', 'FontSize',fs)

ax(3) = subplot(4, 1, 3);
plot(meanAnomlayListRad, diffLastStepList, 'r', 'LineWidth',lw);
ylabel('last $|F_{i+1} - F_i|$ [rad]', 'Interpreter','latex', 'FontSize',fs)

ax(4) = subplot(4, 1, 4);
plot(meanAnomlayListRad, errorKeplersEquationList, 'b', 'LineWidth',lw);
ylabel('final $|e\sinh F - F - M_e|$', 'Interpreter','latex', 'FontSize',fs)
xlabel('input: mean anomaly $M_h$ [rad]', 'Interpreter','latex', 'FontSize',fs)
% axis touchup
for ii = 1:4
    axes(ax(ii))
    grid on;
    if ii > 1
        axis tight;
    end
end
