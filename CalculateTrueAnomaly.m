function trueAnomaly = CalculateTrueAnomaly(rVec, vVec, eVec)
%% Calculate the true anomlay from position, velocity, and eccentricity vectors.

% using dot-product to get the true anomlay
trueAnomaly = acos(dot(rVec, eVec) / norm(rVec) / norm(eVec));

% test whether in upper or lower half of the ellipse
if dot(rVec, vVec) < 0
    trueAnomaly = 2*pi - trueAnomaly; % adjust if in the lower half
end

end