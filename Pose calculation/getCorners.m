function [ NewCorners ] = getCorners( ForcePlateStructure, Corners, StickCOP )
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

NewCorners = zeros(3,4,ForcePlateStructure.NrOfFrames);
force = decimateData(ForcePlateStructure.COP, ForcePlateStructure.SamplingFactor);
for i = 1:ForcePlateStructure.NrOfFrames
    NewCorners(:,:,i) = repmat(StickCOP(:,:,i), 1, 4) - repmat(force,1, 4) + Corners;
end
end

