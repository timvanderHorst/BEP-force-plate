function [COP, d] = getStickVector(R, v, vector)
%getStickVector Given the positions of the markers, calculate the direction
%of the force plate vector in each frame.
%   OUTPUT:
%COP - vector containing the position of the COP in the global reference frame [3x1xT frames]
%d - the stick vector in the global reference frame [3x1xT frames]
%   INPUT:
%R - Rotation matrix
%v - Translation vector
%vector - The stick vector that is to be converted
d = zeros(3, 1, size(R,3));
for i = 1 : size(R,3)
    d(:,:,i) = R(:,:,i) * vector;
end
COP = v;

end

