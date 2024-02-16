function dXdt = OdeTwoBody(mu, X)
%% Two-body motion
%
% Inputs:
%   mu: gravitational constant
%   X: [x y z vx vy vz]^T, row vector of position and velocity [km km/s]
%   
% Outputs:
%   dXdt: derivative
%
% Author: Hao Peng (hao.peng@erau.edu)

r = X(1:3);
v = X(4:6);

drdt = v;

rCube = (sum(r.^2)).^(3/2);
dvdt = - mu / rCube * r;

dXdt = [drdt(:); dvdt(:)]; % trick: A(:) will change any vector to a column vector, equivalent to `reshape(A, [], 1)`.

end