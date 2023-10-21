%% Randomly generate a list of OE, and test covertion between RV and OE

numTesting = 1000;

r = [2500, 16000, 4000]; % [km]
v = [-3, -1, 5]; % [km/s]
h = norm(cross(r,v));

hRandoms = h + h*0.1 * (2*rand(numTesting, 1)-1);
eRandoms = 1 * rand(numTesting, 1);
thetaRandoms = 2*pi * rand(numTesting, 1);
OmegaRAANRandoms = 2*pi * rand(numTesting, 1);
incRandoms = pi * rand(numTesting, 1);
omegaAPRandoms = 2*pi * rand(numTesting, 1);

coeList = [hRandoms, eRandoms, thetaRandoms, OmegaRAANRandoms, incRandoms, omegaAPRandoms];

tic;
for ii = 1:numTesting
    coeTest = coeList(ii, :);
    [tmpR, tmpV] = ConvertCoeToRv(coeTest);
    coeConverted = ConvertRvToCoe(tmpR, tmpV);
    disp(num2str(coeTest - coeConverted, '%+15.4e'));
    if max(abs(coeTest - coeConverted)) > 1e-6
        error('# Hao: Conversion test failed!!!') % one simple failure case is when e>10
    end
end
fprintf('# Hao: conversion test finshed. ')
toc