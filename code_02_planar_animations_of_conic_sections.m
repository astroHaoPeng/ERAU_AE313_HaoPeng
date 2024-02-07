% animation of different sorts of orbits
%
% The code is very rough and is not a good reference for coding sytles.
%

clear all; clc; close all;

muEarth = 398600; % [km^3/s^2]

%% common figure configs
lw = 2;
fs = 10;
ms = 10;
filenameGIF = @(tmpStr)['code_02_orbit_' tmpStr '.gif'];
outerPosition = [200, 100, 400, 400];

%% circular
N = 100;
e = 0;
a = 7000; % [km]
b = a * sqrt(1-e^2);
p = a * (1 - e^2);
% generate data to plot
ttheta = linspace(0, 2*pi, N+1);
rr = p ./ (1 + e * cos(ttheta));
xx = rr .* cos(ttheta);
yy = rr .* sin(ttheta);
% start plotting
MyFigure = figure(21);
set(MyFigure,'OuterPosition',outerPosition, 'MenuBar', 'none');
for ii = 1:N
    % plot
    clf;
    plot(xx, yy, 'k', 'LineWidth',lw); % ellipse
    hold on;
    plot([-a*(1+e), a*(1-e)]+[-0.2, 0.2]*a, [0, 0], 'k--', 'LineWidth',lw/2) % x-axis
    plot([0, 0], [-1, 1]*b+[-0.2, 0.2]*b, 'k--', 'LineWidth',lw/2) % y-axis
    plot(xx(ii), yy(ii), 'ro', 'MarkerFaceColor','r', 'MarkerSize',ms); % m2 point
    plot([0 xx(ii)], [0, yy(ii)], 'r-', 'LineWidth',lw) % position vector
    xlabel('x [km]', 'FontSize',fs)
    ylabel('y [km]', 'FontSize',fs)
    set(gca, 'FontSize',fs)
    axis equal

    % save current frame
    GifData = getframe(gcf);
    if ii == 1
        [X(:,:,1,ii),Map] = rgb2ind(GifData.cdata, 256);
    else
        X(:,:,1,ii) = rgb2ind(GifData.cdata, Map);
    end
end
% write to gif
imwrite(X,Map, filenameGIF('0_circle'), 'LoopCount',Inf,'DelayTime',0);



%% elliptical
clear X;
N = 100;
e = 0.5;
a = 7000; % [km]
b = a * sqrt(1-e^2);
p = a * (1 - e^2);
% generate data to plot
tt = linspace(0, 2*pi, N+1);
ttheta = TimeToTrueAnomalyElliptic(e, tt, 1);
rr = p ./ (1 + e * cos(ttheta));
xx = rr .* cos(ttheta);
yy = rr .* sin(ttheta);
% start plotting
MyFigure = figure(22);
set(MyFigure,'OuterPosition',outerPosition, 'MenuBar', 'none');
for ii = 1:N
    % plot
    clf;
    plot(xx, yy, 'k', 'LineWidth',lw); % ellipse
    hold on;
    plot([-a*(1+e), a*(1-e)]+[-0.2, 0.2]*a, [0, 0], 'k--', 'LineWidth',lw/2) % x-axis
    plot([0, 0], [-1, 1]*b+[-0.2, 0.2]*b, 'k--', 'LineWidth',lw/2) % y-axis
    plot(xx(ii), yy(ii), 'ro', 'MarkerFaceColor','r', 'MarkerSize',ms); % m2 point
    plot([0 xx(ii)], [0, yy(ii)], 'r-', 'LineWidth',lw) % position vector
    xlabel('x [km]', 'FontSize',fs)
    ylabel('y [km]', 'FontSize',fs)
    set(gca, 'FontSize',fs)
    axis equal

    % save current frame
    GifData = getframe(gcf);
    if ii == 1
        [X(:,:,1,ii),Map] = rgb2ind(GifData.cdata, 256);
    else
        X(:,:,1,ii) = rgb2ind(GifData.cdata, Map);
    end
end
% write to gif
imwrite(X,Map, filenameGIF('1_elliptical'), 'LoopCount',Inf,'DelayTime',0);


%% parabolic
clear X;
N = 100*5;
e = 1;
p = 5000; % [km]
% a = inf; % [km]
% b = inf;
% generate data to plot
tt = linspace(-1.5*3600, 1.5*3600, N+1);
ttheta = TimeToTrueAnomalyParabolic(tt, p, muEarth);
% ttheta = linspace(-pi, pi, N+1);
rr = p ./ (1 + e * cos(ttheta));
xx = rr .* cos(ttheta);
yy = rr .* sin(ttheta);
% start plotting
MyFigure = figure(23);
set(MyFigure,'OuterPosition',outerPosition, 'MenuBar', 'none');
kk = 1; % number of frames
for ii = 1:N
    
    % plot
    clf;
    plot(xx, yy, 'k', 'LineWidth',lw); % ellipse
    hold on;
    plot([-12000, 3000], [0, 0], 'k--', 'LineWidth',lw/2) % x-axis
    plot([0, 0], [-12000, 12000], 'k--', 'LineWidth',lw/2) % y-axis
    xlabel('x [km]', 'FontSize',fs)
    ylabel('y [km]', 'FontSize',fs)
    set(gca, 'FontSize',fs)
    xlim([-12000, 3000]);
    ylim([-12000, 12000]);
    if abs(xx(ii)) > 13000
        continue    % skip too large data
    end
    plot(xx(ii), yy(ii), 'ro', 'MarkerFaceColor','r', 'MarkerSize',ms); % m2 point
    plot([0 xx(ii)], [0, yy(ii)], 'r-', 'LineWidth',lw) % position vector
    
    % save current frame
    GifData = getframe(gcf);
    if kk == 1
        [X(:,:,1,kk),Map] = rgb2ind(GifData.cdata, 256);
    else
        X(:,:,1,kk) = rgb2ind(GifData.cdata, Map);
    end
    kk = kk + 1;
end
% write to gif
imwrite(X,Map, filenameGIF('3_parabolic'), 'LoopCount',Inf,'DelayTime',0);


%% hyperbolic
clear X;
N = 150*2;
e = 2;
a = -3000; % [km]
p = a * (1 - e^2); % [km]
b = a * sqrt(e^2 - 1); % [km]
h = sqrt(muEarth * p);
% generate data to plot
tt = linspace(-1.5*3600, 1.5*3600, N+1);
ttheta = TimeToTrueAnomalyHyperbolic(e, tt, h, muEarth);
% ttheta = linspace(-pi, pi, N+1);
rr = p ./ (1 + e * cos(ttheta)); 
xx = rr .* cos(ttheta); 
yy = rr .* sin(ttheta);
% remove vacant trajectory
tmpID = xx < abs(a*(1-e));
xx = xx(tmpID);
yy = yy(tmpID);
rr = rr(tmpID);
% start plotting
MyFigure = figure(24);
set(MyFigure,'OuterPosition',outerPosition, 'MenuBar', 'none');
kk = 1; % number of frames
for ii = 1:length(xx)
    
    % plot
    clf;
    plot(xx, yy, 'k', 'LineWidth',lw); % ellipse
    hold on;
    plot([-12000, 3000], [0, 0], 'k--', 'LineWidth',lw/2) % x-axis
    plot([0, 0], [-30000, 30000], 'k--', 'LineWidth',lw/2) % y-axis
    xlabel('x [km]', 'FontSize',fs)
    ylabel('y [km]', 'FontSize',fs)
    set(gca, 'FontSize',fs)
    xlim([-12000, 3000]);
    ylim([-30000, 30000]);
    if abs(xx(ii)) > 13000
        continue    % skip too large data
    end
    plot(xx(ii), yy(ii), 'ro', 'MarkerFaceColor','r', 'MarkerSize',ms); % m2 point
    plot([0 xx(ii)], [0, yy(ii)], 'r-', 'LineWidth',lw) % position vector
    % axis equal;
    
    % save current frame
    GifData = getframe(gcf);
    if kk == 1
        [X(:,:,1,kk),Map] = rgb2ind(GifData.cdata, 256);
    else
        X(:,:,1,kk) = rgb2ind(GifData.cdata, Map);
    end
    kk = kk + 1;
end
% write to gif
imwrite(X,Map, filenameGIF('4_hyperbolic'), 'LoopCount',Inf,'DelayTime',0);
