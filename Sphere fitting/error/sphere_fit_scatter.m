function [ totalVariance, totalLength, totalCOP, totalError, totalI ] = sphere_fit_scatter(startFrames, endFrames, frames, R, V, identifiedMarkers, map)
%SPHERE_FIT_SCATTER Create a scatter plot of length of length of necessary
%sphere fit vs COP error
%   Performs a sphere fit to calculate the position of the stick tip, and
%   compares this to a measured location for different sphere fit lengths.
% Inputs:
%    path - path to current recording
%    startFrames - 1 x F vector with start frame numbers
%    endFrames - 1 x F vector with end frame numbers
%    frames - 1 x X with correct measurement frame numbers
%    R - rotation matrix used to transform sphere fit location to force
%       plate reference frame
%    V - translation vector used to transform sphere fit location to force
%       plate reference frame
%    identifiedMarkers - M x 4 x N matrix of positions of markers as
%       measured by camera's 
%   map - locations where force is applied on the force plate


totalVariance = [];
totalLength = [];
totalCOP = [];
totalError = [];
totalI = [];
for i = 1 : length(frames)
    len = endFrames(frames(i)) - startFrames(frames(i));
    lenSphereFit = 100 : 100 : len;
    for j = 1 : length(lenSphereFit)
        [COP2 , variance, ~] = Sphere_fit_final(startFrames(frames(i)), startFrames(frames(i)) + lenSphereFit(j), identifiedMarkers);
        totalVariance = [totalVariance, variance];
        totalLength = [totalLength, lenSphereFit(j)];
        totalCOP = [totalCOP, COP2];
        totalError = [totalError, (map(i,:) - transform(COP2',R, V))'];
        totalI = [totalI, i];
    end
end
figure('Name', 'Sphere_fitting_cop_error'); 
hold on
title('Sphere fitting COP error');
xlabel('Length of measurement [s]');
ylabel('Normalized error [mm]')
axis([0 max(totalLength)/100 0 50])
scatter(totalLength/100, vecnorm(totalError),'x'); 
end

