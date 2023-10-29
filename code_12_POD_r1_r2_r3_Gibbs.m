% Gibbs example
%   Curtis2020, Example 5.1

r1Vec = [-294.32, 4265.1, 5986.7].';
r2Vec = [-1365.5, 3637.6, 6346.8].';
r3Vec = [-2940.3, 2473.7, 6555.8].';

[oe, v2Vec] = OrbitDeterminationGibbs(r1Vec, r2Vec, r3Vec)