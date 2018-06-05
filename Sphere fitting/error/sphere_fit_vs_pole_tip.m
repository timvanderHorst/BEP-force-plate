function [ COPSphereFit, varianceSphereFit, meanErrorSphereFit ] = sphere_fit_vs_pole_tip( startFrames, endFrames, stickTip, identifiedMarkers )
%SPHERE_FIT_VS_POLE_TIP Calculates the difference in location calculated by
%sphere fitting and rigid body transformation pole tip calculations
%For each valid measurement frame, calculate the location of the COP as
%determined by sphere fitting, and the difference between this value and
%the mean stick tip location.
%
% Syntax:  [startFrames, endFrames] = getCalibrationFrames(COP)
%
% Inputs:
%    startFrames - a 1 x M vector, containing the starts of each possible
%    measurements
%    endFrames - a 1 x M vector, containing the end frames of each possible
%    measurements
%    frames - 1 x F with correct measurement frame numbers
%    stickTip - 3 x 1 x N vector of position of pole tip in global frame
%    identifiedMarkers - M x 4 x N matrix of positions of markers as
%       measured by camera's 
%
% Outputs:
%   COPSphereFit - 3 x M matrix containing the COP determined by sphere
%   fitting, for every measurement
%   varianceSphereFit - 3 x M matrix containing the variance of the COP
%   meanErrorSphereFit - 3 x M matrix containing mean difference between
%   sphere fitting COP and rigid body transformation COP

COPSphereFit = [];
varianceSphereFit = [];
error = [];
measurementFrames = [];
meanErrorSphereFit = [];
stickTip = squeeze(stickTip);
for i=1:length(startFrames)
    [COP2, variance, centers] = Sphere_fit_final(startFrames(i), endFrames(i), identifiedMarkers);
    measurementFrames = [measurementFrames startFrames(i):endFrames(i)];
    COPSphereFit = [COPSphereFit COP2];
    varianceSphereFit = [varianceSphereFit variance];
    error = [error stickTip(:,startFrames(i):endFrames(i)) - COP2];
    mean_error = mean(stickTip(:,startFrames(i):endFrames(i)) - COP2,2);
    meanErrorSphereFit = [meanErrorSphereFit mean_error];
    
    startmin = (startFrames(i)/6000);
    endmin = (endFrames(i)/6000);
    startsec = floor(60*(startmin - floor(startmin)));
    endsec = floor(60*(endmin - floor(endmin)));
    startmin = floor(startmin);
    endmin = floor(endmin);
    
    fprintf('Frames %d to %d (%2.d:%2.d - %2.d:%2.d) \n',startFrames(i),...
        endFrames(i), startmin, startsec, endmin, endsec);
    fprintf('Average difference between sphere fitting COP and pole COP\n');
    fprintf('X: %fmm\nY: %fmm\nZ: %fmm\n', mean_error);
end

figure('Name', 'comparison_Spherefit_Stick_COP')
subplot(2,1,1)
hold on
plot(measurementFrames, vecnorm(error));
xlabel('Frame number')
ylabel('Distance [mm]')
title('Distance between Sphere Fitting COP and Pose calculation COP')
subplot(2,1,2)
hold on;
plot(measurementFrames, error);
xlabel('Frame number')
ylabel('Distance [mm]')
title('Distance components')
legend('X', 'Y', 'Z')


end

