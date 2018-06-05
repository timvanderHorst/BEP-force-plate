function [ averageCOP, variance ] = getCOP(COP, startFrames, endFrames)
%GET_COP Calculates the average COP value from a rigid body transformation.
%   Given a 3 x 1 x T matrix containing the coordinates of the tip of an
%   instrumented pole at each frame, calculates the average X,Y,Z location
%   of this COP, and returns the variance.
%INPUT: 
%   COP - 3 x 1 x T matrix
%   frameStart - which frame to start measuring from
%   frameEnd - which frame to end the measurement
%OUTPUT:
%   averageCOP - 3 x 1 vector containing location of COP.
%   variance - 3 x 1 vector containing variance in each axis
loc = squeeze(COP);
averageCOP = zeros(3,length(startFrames));
variance = zeros(3,length(startFrames));
for i = 1:length(startFrames)
    loc_short = loc(:,startFrames(i):endFrames(i));
    averageCOP(:,i) = mean(loc_short,2);
    variance(:,i) = var(loc_short,0,2);
end

