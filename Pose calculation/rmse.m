function [ RMSE, E ] = rmse(frames, stickTip, fpCOP, R, V)
    stickTip = squeeze(stickTip);
    stickTip = (R*(stickTip - V));
    E = fpCOP - stickTip;
    Enorm = vecnorm(E);
    RMSE = sum(Enorm(frames))/length(frames);
end

