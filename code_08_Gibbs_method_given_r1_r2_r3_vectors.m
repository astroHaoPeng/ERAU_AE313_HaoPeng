% Gibbs example
%   Curtis2020, Example 5.1

% original values in Example 5.1
% r1Vec = [-294.32, 4265.1 + 0, 5986.7].';
% r2Vec = [-1365.5, 3637.6 + 0, 6346.8].';
% r3Vec = [-2940.3, 2473.7 + 10, 6555.8].';

% Example 5.1 manipulated by adding some errors
% r1Vec = [-294.32, 4265.1 + 0, 5986.7].';
% r2Vec = [-1365.5, 3637.6 + 0, 6346.8].';
% r3Vec = [-2940.3, 2473.7 + 10, 6555.8].';

% Problem 5.1
r1Vec = [5887 -3520 -1204].';
r2Vec = [5572 -3457 -2376].';
r3Vec = [5088 -3289+10 -3480].';

[oe2, v2Vec] = OrbitDeterminationGibbs(r1Vec, r2Vec, r3Vec);
fprintf('#------------------------------\n')
fprintf('# At the 2nd point: \n')
fprintf('#            h = %.10f [km^2/s]\n', oe2(1))
fprintf('#            e = %.10f \n', oe2(2))
fprintf('#        theta = %.10f [rad] = %.10f [deg]\n', oe2(3), rad2deg(oe2(3)))
fprintf('#         RAAN = %.10f [rad] = %.10f [deg]\n', oe2(4), rad2deg(oe2(4)))
fprintf('#  inclination = %.10f [rad] = %.10f [deg]\n', oe2(5), rad2deg(oe2(5)))
fprintf('#           AP = %.10f [rad] = %.10f [deg]\n', oe2(6), rad2deg(oe2(6)))
fprintf('#------------------------------\n')

%% compare results using different points
[oe2, v2Vec, oe1, v1Vec, oe3, v3Vec] = OrbitDeterminationGibbs(r1Vec, r2Vec, r3Vec);
fprintf('#------------------------------\n')
fprintf('# At the 1st point: \n')
PrintResults(oe1)
fprintf('# At the 2nd point: \n')
PrintResults(oe2)
fprintf('# At the 3rd point: \n')
PrintResults(oe3)
fprintf('#------------------------------\n')
function PrintResults(oe_)
fprintf('#            h = %.10f [km^2/s]\n', oe_(1))
fprintf('#            e = %.10f \n', oe_(2))
fprintf('#        theta = %.10f [rad] = %.10f [deg]\n', oe_(3), rad2deg(oe_(3)))
fprintf('#         RAAN = %.10f [rad] = %.10f [deg]\n', oe_(4), rad2deg(oe_(4)))
fprintf('#  inclination = %.10f [rad] = %.10f [deg]\n', oe_(5), rad2deg(oe_(5)))
fprintf('#           AP = %.10f [rad] = %.10f [deg]\n', oe_(6), rad2deg(oe_(6)))
end