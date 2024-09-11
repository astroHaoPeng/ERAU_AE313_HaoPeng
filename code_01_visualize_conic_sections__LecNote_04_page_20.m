% Visualize the intersection between a cone and a plane, which gives generic conic sections.
%
%   - the algorithm is given on page 30 of LecNote_04.
%   - vector doc product is used to define cone and plane
%   - the algorithm can be further improved.
%

%% parameters for the cone
coneOrientationVector = [0, 0, 1].';
coneOrientationDirection = coneOrientationVector / norm(coneOrientationVector);
coneApexAngle = deg2rad(42);
coneApexLocation = [0 0 0].';

%% parameters for the plane
planeNormalVector = [0 0 15].'; % circle
% planeNormalVector = [5 0 15].'; % ellipse
% planeNormalVector = [5 0 1].'; % hyperbola
% planeNormalVector = [0 -cos(coneApexAngle) sin(coneApexAngle)].'; % parabola, must holds an angle of `coneApexAngle` w.r.t. `coneOrientationDirection`
planeNormalDirection = planeNormalVector / norm(planeNormalVector);
planePoint = [3 0 3].';

%% other parameters
% Visualization region [xmin, xmax, ymin, ymax, zmin, zmax]
plotRegion = [-10, 10, -10, 10, -10, 10];

% Use a boolean flag to choose how to define functions, in a vectorized way (better performance) or an intuitive way (slower)
flagVectorizedDefinition = true; 

% Set partial function definitions
%   Here @ is used to define a partial function (fewer inputs) from an original function (more inputs).
funPlaneToPlot = @(x, y, z)PlaneFunction(x, y, z, planePoint, planeNormalDirection, flagVectorizedDefinition); 
funConeToPlot = @(x, y, z)ConeFunction(x, y, z, coneApexAngle, coneApexLocation, coneOrientationDirection, flagVectorizedDefinition); 


%% plotting
figure(23); % initialize a figure window numbered 33
clf; % clear the figure

% plot cone
fimplicit3(funConeToPlot, plotRegion, ...
    'MeshDensity',50, 'EdgeColor','none', 'FaceAlpha',0.3, 'FaceColor','c');
hold on; % hold the plotted surface, otherwise it will be overwritten

% plot plane
fimplicit3(funPlaneToPlot, plotRegion, ...
    'EdgeColor','none', 'FaceAlpha',0.4, 'FaceColor','r');

% plot the normal vector
tmpVec = [planePoint, planePoint + planeNormalDirection * max(plotRegion)/3];
plot3(planePoint(1), planePoint(2), planePoint(3), 'm', 'LineWidth',2, 'Marker','s', 'MarkerSize',15)
plot3(tmpVec(1, :), tmpVec(2, :), tmpVec(3, :), 'm-', 'LineWidth',2);
text(tmpVec(1, 2), tmpVec(2, 2), tmpVec(3, 2), 'Normal Direction of Plane', 'Color','m', 'FontSize',20)

% annotations and touchups
plot3(coneApexLocation(1)+[0,plotRegion(2)/2], coneApexLocation(2)+[0,0], coneApexLocation(3)+[0,0], 'r', 'LineWidth',2) % x-axis in R
plot3(coneApexLocation(1)+[0,0], coneApexLocation(2)+[0,plotRegion(2)/2], coneApexLocation(3)+[0,0], 'g', 'LineWidth',2) % y-axis in G
plot3(coneApexLocation(1)+[0,0], coneApexLocation(2)+[0,0], coneApexLocation(3)+[0,plotRegion(2)/2], 'b', 'LineWidth',2) % z-axis in B
text(coneApexLocation(1)+plotRegion(2)/2, coneApexLocation(2)+0, coneApexLocation(3)+0, 'x', 'Color','r', 'LineWidth',2, 'FontSize',20)
text(coneApexLocation(1)+0, coneApexLocation(2)+plotRegion(2)/2, coneApexLocation(3)+0, 'y', 'Color','g', 'LineWidth',2, 'FontSize',20)
text(coneApexLocation(1)+0, coneApexLocation(2)+0, coneApexLocation(3)+plotRegion(2)/2, 'z', 'Color','b', 'LineWidth',2, 'FontSize',20)
axis equal; % same aspect ratio
xlabel('x');
ylabel('y');
zlabel('z');




%% function definitions

function outputValue = PlaneFunction(x, y, z, pointVector, normalDirection, flagVectorizedDefinition)
% The equation a point on a plance satisfying. 
if flagVectorizedDefinition
    % Vectorized definition.
    dx = x - pointVector(1);
    dy = y - pointVector(2);
    dz = z - pointVector(3);
    outputValue = dx.*normalDirection(1) + dy.*normalDirection(2) + dz.*normalDirection(3);
else
    % dot product definition
    A = [x, y, z].';
    outputValue = dot(A - pointVector, normalDirection);
end
end


function outputValue = ConeFunction(x, y, z, alpha, pointVector, orientationDirection, flagVectorizedDefinition)
% The equation a point on a cone satisfying. 
if flagVectorizedDefinition
    % Vectorized definition.
    dx = x - pointVector(1);
    dy = y - pointVector(2);
    dz = z - pointVector(3);
    tmpNorm = (dx.^2 + dy.^2 + dz.^2).^(1/2);
    outputValue = abs((dx.*orientationDirection(1) + dy.*orientationDirection(2) + dz.*orientationDirection(3)) ./ tmpNorm) - cos(alpha);
else
    % dot product definition
    A = [x, y, z].';
    tmpPA = A - pointVector;
    normTmpPA = sqrt(tmpPA(1)^2 + tmpPA(2)^2 + tmpPA(3)^2);
    outputValue = abs(dot(tmpPA, orientationDirection)) / normTmpPA  - cos(alpha);
end
end
