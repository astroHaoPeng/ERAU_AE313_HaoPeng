%% define anonymous functions using the symbole @
%   Check this link for more help:
%   - Anonymous Functions: https://www.mathworks.com/help/matlab/matlab_prog/anonymous-functions.html 

clear; clc;

a1 = 1;
a2 = 1;
a3 = 1;

% an arbitrarily defined left hand side of the equation
equationLHS = @(x) 0 + 1.2.*x + 0.3*x.^2 - 0.045*x.^3;

% an arbitrarily defined right hand side of the equation
equationRHS = @(x, y, z) a1.*sin(x) + a2.*y.^2 + a3.*exp(z);

% visualize `LHS`, `RHS`, and `LHS-RHS` to identify possible solutions
figure(401); clf;

tmp = linspace(-6, 10, 100);

subplot(2,1,1);
plot(tmp, equationLHS(tmp), 'r', 'DisplayName','LHS'); hold on;
plot(tmp, equationRHS(tmp, 0.5, 0.8), 'g', 'DisplayName','RHS');
legend;

subplot(2,1,2);
plot(tmp, equationLHS(tmp) - equationRHS (tmp, 0.5, 0.8), 'b', 'DisplayName','LHS - RHS'); hold on;
plot(tmp, tmp*0, 'b--', 'DisplayName','zero line')
legend;


%% using fzero to find the root of the equation

funToSolve = @(z) equationLHS(z) - equationRHS(z, 0.5, 0.8);
fzero(funToSolve, -4)
fzero(funToSolve, 3)
fzero(funToSolve, 10)
