function trueAnomaly = TimeToTrueAnomalyParabolic(time, p, muCenteral)

%% t --> Mp
Mp = sqrt(muCenteral/p^3) * time; %sqrt(muEarth/p^3) = sqrt(muEarth*p)/p^2 = h/p^2 = muEarth^2/h^3

%% Mp --> z
z = (3*Mp + (1 + (3*Mp).^2).^(1/2) ).^(1/3);

%% z --> theta
trueAnomaly = 2 * atan( z - 1./z );
trueAnomaly = mod(trueAnomaly, 2*pi); % bring to [0, 2*pi) range

end