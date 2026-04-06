%% define anonymous functions using the symbole @
%   Check this link for more help:
%   - Anonymous Functions: https://www.mathworks.com/help/matlab/matlab_prog/anonymous-functions.html 

clear; clc;

a1 = 1;
a2 = 1;
a3 = 1;

% an arbitrarily defined left hand side of the equation
equationLHS = @(x) 0 + 1.2.*x + 0.3*x.^2 - 0.045*x.^3;
equationLHS_latex = '$1.2 x + 0.3 x^2 - 0.045 x^3$';

% an arbitrarily defined right hand side of the equation
equationRHS = @(x, y, z) a1.*sin(x) + a2.*y.^2 + a3.*exp(z);
equationRHS_latex = '$a_1 \sin x + a_2 y^2 + a_3 z^2$';

% visualize `LHS`, `RHS`, and `LHS-RHS` to identify possible solutions
figure(401); clf;

tmp = linspace(-6, 10, 100);

subplot(2,1,1);
plot(tmp, equationLHS(tmp), 'r', 'DisplayName','LHS'); 
hold on;
plot(tmp, equationRHS(tmp, 0.5, 0.8), 'g', 'DisplayName','RHS');
ylim([-5, 15]);
legend;

subplot(2,1,2);
plot(tmp, equationLHS(tmp) - equationRHS (tmp, 0.5, 0.8), 'b', 'DisplayName','LHS - RHS'); 
hold on;
plot(tmp, tmp*0, 'b--', 'DisplayName','zero line')
ylim([-5, 15]);
legend;


%% using fzero to find the root of the equation

funToSolve = @(z) equationLHS(z) - equationRHS(z, 0.5, 0.8);
root1 = fzero(funToSolve, -4);
root2 = fzero(funToSolve, 3);
root3 = fzero(funToSolve, 10);

subplot(211);
plot([root1, root1], [-5, 15], 'k:', 'LineWidth', 2);
plot([root2, root2], [-5, 15], 'k:', 'LineWidth', 2);
plot([root3, root3], [-5, 15], 'k:', 'LineWidth', 2);
title({'Curves of LHS and RHS of the equation:', [equationLHS_latex ' $=$ ' equationRHS_latex 'with fixed $y=0.5$ and $z=0.8$']}, 'Interpreter', 'latex')

subplot(212);
plot(root1, 0, 'rx', 'MarkerSize', 20, 'LineWidth', 2)
plot(root2, 0, 'gx', 'MarkerSize', 20, 'LineWidth', 2)
plot(root3, 0, 'bx', 'MarkerSize', 20, 'LineWidth', 2)
title('Roots of the above equation.')

