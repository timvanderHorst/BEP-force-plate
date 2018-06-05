function M2 = transform(M, R, V)
%TRANSFORM transforms N x 3 matrx M in coordinate system 1 to coordinate
%system 2
%   INPUT -
% V = Translation vector (1 x 3)
% R = Rotation matrix 
% M = N x 3 matrix.

M2 = (R*(M - V)')';
end

