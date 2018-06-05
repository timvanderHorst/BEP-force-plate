function [ minIndex ] = height_correction(plate, data, stickTip, R, V, f)
%HEIGHT_CORRECTION Calculates the optimal h-offset value for a given
%dataset.
%Uses values of h between -100 and 0mm to minimize the distance between the
%tip of the pole and the force plate data.
% Inputs:
%    plate - plate id ('a', 'c')
%    data - N x 8 analog data for all frames
%    stickTip - 3 x 1 x N vector of position of pole tip in global frame
%    R - rotation matrix used to transform sphere fit location to force
%       plate reference frame
%    V - translation vector used to transform sphere fit location to force
%       plate reference frame
%    f - 1 x X vector with frame numbers
% Outputs:
%    minIndex - best h-offset value in mm

stickCOP = transform(squeeze(stickTip(:,:,f))',R, V);
min = 10e10;

for i = -100:0
    [analogFpCOP,~] = analogDataAnalysis(plate, data, i);
    q = sum(vecnorm((analogFpCOP(:,f)' - stickCOP)'));
    if (q < min)
        min = q;
        minIndex = i;
    end
end
fprintf('Found h that minimizes rmse: %dmm\n',minIndex); 

end

