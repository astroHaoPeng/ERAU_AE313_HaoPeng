% Preview of using scalar and vector products to define cone and plane, 
% then visualize them to find out the conic sections.
%
% Details will be covered in the next lecture, lecture_07.

%% parameters

% parameters for the plane
planeNormalVector = [0 0 15].'; % circle
% planeNormalVector = [5 0 15].'; % ellipse
% planeNormalVector = [5 0 1].'; % hyperbola
% planeNormalVector = [0 0 0].'; % parabola
planeNormalDirection = planeNormalVector / norm(planeNormalVector);
planePoint = [3 0 3].';

% parameters for the cone
coneNormalVector = [0, 0, 1].';
coneNormalDirection = coneNormalVector / norm(coneNormalVector);
coneApexAngle = deg2rad(42);
coneApexLocation = [0 0 0].';

% visualization region [xmin, xmax, ymin, ymax, zmin, zmax]
plotRegion = [-10, 10, -10, 10, -10, 10];

% others
flagVectorizedDefinition = false; % choose how we define functions

%%%%%%%%%
% use space and brackets properly --> as possible as you can
%
% pick a naming rule:
%  - lowerCamelCased:       Variable names
%  - CamelCased:            User-defined functions and methods
%  - SCREAMING_SNAKE_CASED: CONSTANT variable
%
% be consistent and flexible
%%%%%%%%%

%% set proper function definitions
if flagVectorizedDefinition
    % vectorized definition
    funPlaneToPlot = @(x,y,z)PlaneVectorized(x,y,z, planePoint, planeNormalDirection); 
    funConeToPlot = @(x,y,z)ConeVectorized(x,y,z, coneApexAngle, coneApexLocation, coneNormalDirection); 
else
    % simple definition
    funPlaneToPlot = @(x,y,z)PlaneSimple(x,y,z, planePoint, planeNormalDirection); 
    funConeToPlot = @(x,y,z)ConeSimple(x,y,z, coneApexAngle, coneApexLocation, coneNormalDirection); 
end


%% do the actual plotting
figure(33); % initialize a figure window numbered 33
clf; % clear the figure

% plot cone
fimplicit3(funConeToPlot, plotRegion, ...
    'MeshDensity',50, 'EdgeColor','none', 'FaceAlpha',0.3, 'FaceColor','c');
hold on; % hold the plotted surface, otherwise it will be overwritten

% plot plane
fimplicit3(funPlaneToPlot, plotRegion, ...
    'EdgeColor','none', 'FaceAlpha',0.4, 'FaceColor','r');
hold on; % hold the plotted surface, otherwise it will be overwritten

% plot the normal vector
tmpVec = [planePoint, planePoint + planeNormalDirection * max(plotRegion)/3];
plot3(planePoint(1), planePoint(2), planePoint(3), 'Marker','s')
plot3(tmpVec(1, :), tmpVec(2, :), tmpVec(3, :), 'b-');

% annotations and touchups
plot3(coneApexLocation(1)+[0,4], coneApexLocation(2)+[0,0], coneApexLocation(3)+[0,0], 'r', 'LineWidth',2) % x-axis in R
plot3(coneApexLocation(1)+[0,0], coneApexLocation(2)+[0,4], coneApexLocation(3)+[0,0], 'g', 'LineWidth',2) % y-axis in G
plot3(coneApexLocation(1)+[0,0], coneApexLocation(2)+[0,0], coneApexLocation(3)+[0,4], 'b', 'LineWidth',2) % z-axis in B
text(coneApexLocation(1)+5, coneApexLocation(2)+0, coneApexLocation(3)+0, 'x', 'Color','r', 'LineWidth',2, 'FontSize',20)
text(coneApexLocation(1)+0, coneApexLocation(2)+5, coneApexLocation(3)+0, 'y', 'Color','g', 'LineWidth',2, 'FontSize',20)
text(coneApexLocation(1)+0, coneApexLocation(2)+0, coneApexLocation(3)+5, 'z', 'Color','b', 'LineWidth',2, 'FontSize',20)
axis equal; % same aspect ratio
xlabel('x');
ylabel('y');
zlabel('z');




%% function definitions

function fVal = PlaneSimple(x, y, z, point, nDirection)
% The equation a point on a plance satisfying. 
A = [x, y, z].';
fVal = dot(A - point, nDirection);
end


function fVal = PlaneVectorized(x, y, z, point, nDirection)
% The equation a point on a plance satisfying. 
% Vectorized definition.
dx = x - point(1);
dy = y - point(2);
dz = z - point(3);
fVal = dx.*nDirection(1) + dy.*nDirection(2) + dz.*nDirection(3);
end


function fVal = ConeSimple(x, y, z, alpha, point, nDirection)
% The equation a point on a cone satisfying. 
A = [x, y, z].';
tmpPA = A - point;
normTmpPA = sqrt(tmpPA(1)^2 + tmpPA(2)^2 + tmpPA(3)^2);
fVal = abs(dot(tmpPA, nDirection)) / normTmpPA  - cos(alpha);
end


function fVal = ConeVectorized(x, y, z, alpha, point, nDirection)
% The equation a point on a cone satisfying. 
% Vectorized definition.
dx = x - point(1);
dy = y - point(2);
dz = z - point(3);
tmpNorm = (dx.^2 + dy.^2 + dz.^2).^(1/2);
fVal = abs((dx.*nDirection(1) + dy.*nDirection(2) + dz.*nDirection(3)) ./ tmpNorm) - cos(alpha);
end
