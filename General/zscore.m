function [ score ] = zscore( forceVector, d, f )
%ZSCORE - Perform Z-test over relation between force > mean force and 
%angular error
% Syntax:  zscore( forceVector, d, f )
%
% Inputs:
%    forceVector - 3 X N vector containing force plate force vector
%    d - 3 x 1 x N vector containing vector along pole's central axis
%    f - 1 x X vector with frame numbers
% Outputs:
%    score - z-score value

angle = vectorAngle(forceVector, d);
force_norm = vecnorm(forceVector);
validFrames = zeros(length(forceVector),1);
validFrames(f) = 1; validFrames = logical(validFrames);
force_norm_greater = force_norm > mean(force_norm(f));
angle_greater = angle(and(force_norm_greater,validFrames));
angle_less = angle(and(not(force_norm_greater),validFrames));
fprintf('Average of angle where force > average force: %f, sd: %f\n', mean(angle_greater), sqrt(var(angle_greater)));
fprintf('Average of angle where force <= average force: %f, sd: %f\n', mean(angle_less), sqrt(var(angle_less)));
score = sqrt(length(f))*norm(mean(angle_greater) - mean(angle_less))/max(sqrt(var(angle_greater)), sqrt(var(angle_less)));
fprintf('Z score: %f\n', score);

end

