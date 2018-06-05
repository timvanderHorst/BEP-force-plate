function errorData = total_rmse_calc( filename, f, frames, startFrames, endFrames, stickTip, d, forceCOP, forceVector, R1final, V1final, CornerR, CornerV, Rtest, Vtest )
COP_RMSE = [rmse(f, stickTip, forceCOP, R1final, V1final), ...
    rmse(f, stickTip, forceCOP, mean(CornerR,3), mean(CornerV,3)'),...
    rmse(f, stickTip, forceCOP, Rtest, Vtest)];
fprintf('COP\n')
fprintf('RMSE of rabuffetti using fminsearch: %9.6fmm\n', COP_RMSE(1));
fprintf('RMSE of corner calibration: %9.6fmm\n', COP_RMSE(2));
fprintf('RMSE of challis: %9.6fmm\n', COP_RMSE(3));

angle = vectorAngle(forceVector, R1final*squeeze(-1*d));
angle3 = vectorAngle(forceVector, mean(CornerR,3)*squeeze(-1*d));
angle4 = vectorAngle(forceVector, Rtest*squeeze(-1*d));    
anglesRMSE = zeros(3,length(frames));
forceTotal = zeros(3,length(frames));
COP_error = zeros(3, length(frames));
measurementLength = zeros(length(frames),1);
anglestotalRMSE = [(sum(angle(f).^2)/length(f))^0.5, ...
    (sum(angle3(f).^2)/length(f))^0.5, ...
    (sum(angle4(f).^2)/length(f))^0.5];
for i = 1 : length(frames)
    fr = startFrames(frames(i)) : endFrames(frames(i));
    anglesRMSE(1,i) = (sum(angle(fr).^2)/length(fr))^0.5;
    anglesRMSE(3,i) = (sum(angle3(fr).^2)/length(fr))^0.5;
    anglesRMSE(4,i) = (sum(angle4(fr).^2)/length(fr))^0.5;
    COP_error(:,i) = mean(mean(CornerR,3)*(squeeze(stickTip(:,:,fr)) - mean(CornerV,3)') - forceCOP(:,fr),2);
    forceTotal(:,i) = mean(forceVector(:,fr),2);
    measurementLength(i) = length(fr);
end

errorData = struct('angles',anglesRMSE,'anglesTotal',anglestotalRMSE,'length',measurementLength,'COP_error', COP_error, 'force',forceTotal,'COP_RMSE',COP_RMSE);
save(strcat(filename,'errorData2.mat'),'errorData')
fprintf('Angle\n')
fprintf('RMSE of angle using fminsearch: %9.6f%c\n', anglestotalRMSE(1), char(176));
fprintf('RMSE of angle using corner calibration: %9.6f%c\n', anglestotalRMSE(2), char(176));
fprintf('RMSE of angle using challis: %9.6f%c\n', anglestotalRMSE(3), char(176));
end

