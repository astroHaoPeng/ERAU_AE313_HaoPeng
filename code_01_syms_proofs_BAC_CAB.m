%% Use symbolic calculation to prove BAC-CAB rule for vector triple products

% declare symbolic variables as real variables.
syms Ax Ay Az Bx By Bz Cx Cy Cz real

% define three vectors
A = [Ax Ay Az].';
B = [Bx By Bz].';
C = [Cx Cy Cz].';

% built-in dot product function
dot(A, B)

% alternative ways to do dot product
A.' * B % 1x3 and 3x1 matrix multiplication
B.' * A % same, and dot product is commutative

% built-in cross product function
cross(A, B)

% for triple vector product, MATLAB does not take `cross(A, B, C)`
% but we have to express it out explicitely
cross(A, cross(B, C))

% dot product function only support inputs with same dimensions
% so `dot(B, dot(A, C))` won't work because the inner term `dot(A, C)` is a scaler
% the correct form is to use `*` for scaler multiplication, which is
B * dot(A, C) - C * dot(A, B)

% LHS simplified
simplify(expand(cross(A, cross(B, C))))

% RHS simplified
simplify(B * dot(A, C) - C * dot(A, B))

% subtract them directly, and simplify if necessary
cross(A, cross(B, C)) - ( B * dot(A, C) - C * dot(A, B) )
simplify( cross(A, cross(B, C)) - ( B * dot(A, C) + C * dot(A, B) ) )
expand( cross(A, cross(B, C)) - ( B * dot(A, C) - C * dot(A, B) ) )
