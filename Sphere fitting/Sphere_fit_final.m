function [ COP, variance, centers ] = Sphere_fit_final(startFrame, endFrame, identifiedMarkers )
%SPHERE_FIT Computes the center of rotation of a set of markers
%   [COP, variance, centers] = SPHERE_FIT(filename, startFrame, endFrame)
%OUTPUT:
%   COP - the calculated location of the COP, a 1x3 matrix
%   variance - the variance of the centers for all markers
%   centers - the calculated centers for each marker
%INPUT:
%   filename - path to the recording.mat file, exported by QTM
%   startFrame - the frame the sphere fitting is to be started
%   endFrame - the frame the spehere fitting should end

markers = zeros(size(identifiedMarkers,1),endFrame - startFrame,3);

%Reshape array Nx4xM array (N = Number of markers, M = Number of
%measurements) to an N x M x 3 matrix
for i = 1 : endFrame - startFrame
    markers(:,i,:) = identifiedMarkers(:,1:3,i + startFrame);
end

centers = zeros(size(identifiedMarkers,1),3);
radius = zeros(size(identifiedMarkers,1),1);

for i = 1 : size(identifiedMarkers,1)
    [centers(i,:), radius(i)] = sphereFit(squeeze(markers(i,:,:)));
end
COP = mean(centers)';
variance = var(centers)';
end