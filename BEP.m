%Force plate analysis - 15/03/2018
close all;

% Force plate dimensions
%xWidth = 500.0;
%yWidth = 600.0;
%yWidth = 298.5;

% Define marker positions in reference frame
m1 = [0; -160; 200];
m2 = [0; 160; 200];
m3 = [0; 160; 400];
m4 = [0; -160; 800];
m5 = [0; 160; 800];

% Define pole mass
mass = 1.4; %kg

% Plate thickness
plateThickness = 2.4; %mm

% Constants
g = 9.81;


% Load recording. If there is a variable called 'recording' in scope, ask
% if the user wants to use this variable.
if ((exist('recording') + exist('file') + exist('path')) ~= 3)
    [file,path] = uigetfile('*.mat', 'Main recording')
    recording = load(strcat(path, file))
    rec = fieldnames(recording);
    rec = rec{1};
    
    % Check if the file is in local or global coordinates.
    worldCoordinates = input('Is this file using World (lab) coordinates? Y/N [Y]: ','s');
    worldCoordinates = ~(~isempty(worldCoordinates) && any(lower(worldCoordinates) == 'n'));
else
    x = input('Would you like to use the current file: ? Y/N [Y]: ','s');
    if(~isempty(x) && any(lower(x) == 'n'))
        [file,path] = uigetfile
        recording = load(strcat(path, file))
        rec = fieldnames(recording);
        rec = rec{1};
        worldCoordinates = input('Is this file using World (lab) coordinates? Y/N [Y]: ','s');
        worldCoordinates = ~(~isempty(worldCoordinates) && any(lower(worldCoordinates) == 'n'));
    end
end


% The coordinates + residual of every marker, for every frame
identifiedMarkers = recording.(rec).Trajectories.Labeled.Data;

%Define marker positions in reference frame. This function finds out what
%ordering of markers is used in the camera data by comparing relative
%distances between markers.
i = 1;
while any(any(isnan(identifiedMarkers(:,1:3,i)))) || ~exist('a')
    i = i + 1;
    try
        [a, markerNumbers] = getMarkerLocations(identifiedMarkers(:,1:3,i)', [m1 m2 m3 m4 m5]);
    catch ME
    fprintf('Finding markers failed, trying next frame\n');
    end
end
clear i
%% Geometric stick tip
topMarkers = identifiedMarkers([find(markerNumbers == 5), find(markerNumbers == 4)],:,:);
bottomMarkers = identifiedMarkers([find(markerNumbers == 1), find(markerNumbers == 2)],:,:);
meanTop = squeeze(mean(topMarkers,1));
meanBottom = squeeze(mean(bottomMarkers,1));
d2 = meanTop(1:3,:) - meanBottom(1:3,:);
d2 = d2./vecnorm(d2)';
vGeometric = meanBottom(1:3,:) - 200*d2;
% vGeometric = vGeometric(:,find(valid == 1));
%
% filterfreq = 25;  samplefreq = 100;  [w,q] = butter(3,filterfreq/(samplefreq/2)); newV = filtfilt(w,q,vGeometric')';
%
% figure
% plot(newV' - squeeze(vChallis(:,find(valid == 1)))')
%% Veldpaus, 1988 and Rose
% %Perform rigid body transformation using veldpaus' method. Tries to
% %calculate a rotation matrix and translation vector that maps the markers
% %in the reference frame defined above (variable a) to the markers in the
% %global coordinate system. This function returns a Rotation matrix and
% %translation vector for every frame.
[RVeldpaus_1, vVeldpaus_1, valid] = spoorVeldpausRotation(a, identifiedMarkers);

%Calculate the RMS, by mapping the markers in the global frame back onto
%the reference frame, and calculating the rms of the total distance between
%the mapped markers and the markers in the reference frame.
[EVeldpaus_1] = rms(a, identifiedMarkers, RVeldpaus_1, vVeldpaus_1, valid);
fprintf('RMSE Veldpaus eigenvector decomposition: %d\n', sum(EVeldpaus_1 .* valid) / sum(valid));

[RVeldpaus_2,vVeldpaus_2] = spoorVeldpausRotationSVD(a, identifiedMarkers);
[EVeldpaus_2] = rms(a, identifiedMarkers, RVeldpaus_2, vVeldpaus_2, valid);
fprintf('RMSE using Veldpaus SVD: %d\n', sum(EVeldpaus_2 .* valid) / sum(valid));

%% Challis, 1995
[RChallis,vChallis, valid] = challis(a, identifiedMarkers);
[EChallis] = rms(a, identifiedMarkers, RChallis, vChallis, valid);
fprintf('RMSE using Challis SVD: %d\n', sum(EChallis .* valid) / sum(valid));
% figure
% plot(EChallis)
% title('RMS of Challis (1995)')
% xlabel('RMS [mm]')
% ylabel('Frame number [-]')
%% Data analysis
%Get the location of the COP and the vector pointing along the pole's
%central axis in the global frame, by mapping the vector along the pole's
%central axis in the reference frame to the global frame. The tip of the
%stick is defined as [0,0,0] in the reference frame and is therefore equal
%to the translation vector.
[stickTip, d] = getStickVector(RChallis, vChallis, [0; 0; 1]);

%Find the frames in which the pole's tip are kept stationary, indicating a
%calibration is taking place
[startFrames,endFrames] = getCalibrationFrames(stickTip);
for i = 1 : length(startFrames)
    startFrames(i) = startFrames(i) + 300;
    endFrames(i) = endFrames(i) - 120;
end
f = [];
for i = 1 : length(frames)
    f = [f startFrames(frames(i)):endFrames(frames(i))];
end

%% Zero force offset
% Calculates the zero force offset of the force plate, and subtracts this
% from all force data. This is done by calculating the mean force during
% a time period specified by the user.

%Resample vectors to match camera sampling rate.
forceVector = decimateData(recording.(rec).Force.Force,recording.(rec).Force.SamplingFactor);
forceCOP = decimateData(recording.(rec).Force.COP,recording.(rec).Force.SamplingFactor);


% % Check if there are any events defined in QTM which could indicate a
% % zero-offset measurement.
% events = recording.(rec).Events;
% strings = {'zero','empty','offset','leeg'};
% clear zero_start;
% for i = 1 : size(events,1)
%     event = events(i);
%     if any(strcmp(strings, lower(event.Label)))
%         frame = int2str(event.Frame);
%         zero_start = input(strcat('Use event:',' ', lower(event.Label),' ','with frame start: ', frame,' ', '[Y] ?: '));
%         if(isempty(zero_start) || strcmp(lower(zero_start), 'y'))
%             zero_start = event.Frame;
%         end
%     end
% end
% 
% %Ask user for input if no variable has been entered.
% if ~exist('zero_start')
%     zero_start = input('Please input start frame of zero offset: ');
% end
% zero_end = input('Please input end frame of zero offset [END]: ');
% if(isempty(zero_end) || any(zero_end == 'END'))
%     zero_end = length(forceVector);
% end
% 
% %Subtract mean force between the zero_start and zero_end frames
% forceZeroffset = forceVector(:,zero_start:zero_end);
% meanOffset = mean(forceZeroffset,2);
% forceVector = forceVector - meanOffset;

%% Save measurement frames variables to memory.
% if ~isfield(recording.(rec),'MeasurementFrames')
    frames = [];
    THRESHOLD_Z_FORCE = 15;
    THRESHOLD_LENGTH = 20; %%5s measurement
    THRESHOLD_Z = 200;
    for i = 1:length(startFrames)
        p = mean(stickTip(:,:,startFrames(i):endFrames(i)),3);
        f = mean(forceVector(:,startFrames(i):endFrames(i)),2);
        if abs(f(3)) > THRESHOLD_Z_FORCE && (endFrames(i) - startFrames(i)) > THRESHOLD_LENGTH && p(3) < THRESHOLD_Z
            frames = [frames i];
        else
            if abs(f(3)) < THRESHOLD_Z_FORCE
                fprintf('Force in Z direction (%3.f) during frames %d - %d below threshold \n', f(3), startFrames(i), endFrames(i));
            elseif p(3) > THRESHOLD_Z
                fprintf('Invalid stick COP <%3.f %3.f %3.f> \n',p);
            else
                fprintf('Measurement from frames %d - %d too short \n',startFrames(i), endFrames(i));
            end
            
        end
    end
    string = sprintf('%d ',frames);
    x = input(strcat('Would you like to use these valid frames: [', string, '] [Y]/N? '),'s');
    if(~isempty(x) && any(lower(x) == 'n'))
        frames = input('Please input frames:\n');
    end
    recording.(rec).MeasurementFrames = frames;
    fprintf('Saving frames to file\n')
    data = recording.(rec);
    save(strcat(path,file), 'data')
    clear data
    fprintf('Done saving \n')
    
% else
%     frames = recording.(rec).MeasurementFrames;
%     string = sprintf('%d ',frames);
%     fprintf(strcat('Using these frames from recording: [', string, ']\n'));
% end
% 
clear f THRESHOLD_Z_FORCE THRESHOLD_LENGTH string

%% Force plate location
%Calculate the positions of the tip of the stick and the COP of the force
%plate during the frames which were calculated above.
[meanStickCOP, varianceStick] = getCOP(stickTip, startFrames, endFrames);
[meanFpCOP, varianceFp] = getCOP(forceCOP,startFrames, endFrames);

% Use challis to calculate a rotation and translation matrix which would
% map the two coordinate systems onto each other. This is done by entering
% the FpCOP values as the reference frame positions and the tip of the
% stick location as the global position. Uses the locations in each of the
% frames specified above to perform this calculation.
[Rtest,Vtest, validtest] = challis(meanFpCOP(:,frames), [meanStickCOP(:,frames)' zeros(size(meanStickCOP(:,frames),2),1)]);

% Calculate the rotation matrix and translation vector of the corners, by
% asking the user to input the files they would like to analyze, and
% performing the procedure described in Rabuffetti, 2012
[CornerMean, CornerVar,positions, CornerR, CornerV] = cornerCalibration(path, 9);

%% Rabuffetti, 2002
% T = transformation vector
% T = zeros(12, length(startFrames));
T = zeros(6, 1);
%Rstart = [-1 0 0; 0 1 0; 0 0 -1];
%Rstart = eul2rotm([1.5 -3 -1], 'ZYZ');

Rstart = Rtest;
Vstart = Vtest';
% T(:) = rabuffetti( d(:,:,f), ...
%     stickTip(:,:,f), ...
%     -1*forceVector(:,f), ...
%     forceCOP(:,f), ...
%     Rstart,Vstart');
% 
% Rfinal = eul2rotm(T(1:3)','ZYZ');
% Vfinal = T(end-2:end);
% fprintf('\n\n                         fmincon\n');
% printOptimisationResults(Rstart,Rfinal,Vstart,Vfinal);
Rfinal = Rtest;
Vfinal = Vstart';
T1 = T;
T1 = minsearch( d(:,:,f), ...
    stickTip(:,:,f), ...
    -1*forceVector(:,f), ...
    forceCOP(:,f), ...
    Rstart,Vstart');
R1final = eul2rotm(T1(1:3)','ZYZ');
V1final = T1(end-2:end);

fprintf('\n\n                         fminsearch\n');
printOptimisationResults(Rstart,R1final,Vstart,V1final);
%% Plot figures

%If using world coordinates, force = opposite direction force
forceReal = forceVector;

if(worldCoordinates)
    forceReal = forceReal*-1;
else
    
end
angle = vectorAngle(forceReal, d);
angleMean = [];
angleVar = [];
for i = 1:length(startFrames)
    angleMean = [angleMean mean(angle(startFrames(i) : endFrames(i)))];
    angleVar = [angleVar var(angle(startFrames(i) : endFrames(i)))];
end

%Projection of stick vector onto XY plane
%proj(V onto plane with normal N) = v - N*(v.N)/||N||^2
projectStick = d - dot(d,repmat([0; 0; 1],1,1,size(d,3))).*repmat([0; 0; 1],1,1,size(d,3));
angle2 = vectorAngle(projectStick, repmat([0;1;0],1,1,size(d,3)));

figure('Name','angle_forcevector_polevector')
hold on;
plot(1:size(identifiedMarkers,3), angle);
plot(1:size(identifiedMarkers,3), angle2);
xlabel('Frame number [-]')
ylabel('Angle [deg]')
title('Angle between Force vector and pole vector')
legend('Pole - force vector','projected pole onto XY plane - y axis')

figure('Name','cop_location_fp')
hold on;
plot(1:length(forceCOP), forceCOP);
legend('COP X','COP Y', 'COP Z')
ylabel('Coordinates [mm]')
xlabel('Frame number')
title('Location of COP in force plate reference frame')

filterfreq = 25;  samplefreq = 1000;  [w,q] = butter(3,filterfreq/(samplefreq/2)); newforceCOP = filtfilt(w,q,forceCOP')';
filterfreq = 25;  samplefreq = 1000;  [w,q] = butter(3,filterfreq/(samplefreq/2)); newForce = filtfilt(w,q,forceReal')';

figure('Name','filtered_cop_fp_stick')
hold on;
plot(1:length(newforceCOP(:,f)), newforceCOP(:,f));

R = mean(CornerR,3);
V = mean(CornerV,3);
Q = squeeze(stickTip); Q = R * (Q - V');
plot(1:length(Q(:,f)), Q(1:2,f));
legend('COP X','COP Y', 'COP Z', 'Stick X', 'Stick Y')
ylabel('Coordinates [mm]')
xlabel('Frame number')
title('Location of FILTERED COP in force plate reference frame')

figure('Name','filtered_force_fp')
hold on;
plot(1:length(newForce(:,f)), newForce(:,f));
legend('Fx','Fy', 'Fz')
title('FILTERED Force plate XYZ force components');
ylabel('Force [N]')
xlabel('Frame number')

figure('Name','unfiltered_stick_tip')
hold on;
plot(1:length(stickTip), squeeze(stickTip));
legend('X','Y', 'Z')
title('Location of tip of stick in global reference frame');
ylabel('Coordinates [mm]')
xlabel('Frame number')

figure('Name','unfiltered_force_fp')
hold on;
plot(1:length(forceReal), forceReal);
legend('Fx','Fy', 'Fz')
title('Force plate XYZ force components');
ylabel('Force [N]')
xlabel('Frame number')

figure('Name','stick_vector')
hold on;
plot(1:length(d), squeeze(d));
legend('Fx','Fy', 'Fz')
title('Stick vector xyz components');
ylabel('-')
xlabel('Frame number')
%% Calculate error between stick tip and 
if(~worldCoordinates)
    COP_RMSE = [rmse(f, stickTip, forceCOP, R1final, V1final), ...
        rmse(f, stickTip, forceCOP, Rfinal, Vfinal),...
        rmse(f, stickTip, forceCOP, mean(CornerR,3), mean(CornerV,3)'),...
        rmse(f, stickTip, forceCOP, Rtest, Vtest)];
    fprintf('COP\n')
    fprintf('RMSE of rabuffetti using fminsearch: %9.6fmm\n', COP_RMSE(1));
    fprintf('RMSE of rabuffetti using fmincon: %9.6fmm\n', COP_RMSE(2));
    fprintf('RMSE of corner calibration: %9.6fmm\n', COP_RMSE(3));
    fprintf('RMSE of challis: %9.6fmm\n', COP_RMSE(4));
   
    angle = vectorAngle(forceVector, R1final*squeeze(-1*d));
    angle2 = vectorAngle(forceVector, Rfinal*squeeze(-1*d));
    angle3 = vectorAngle(forceVector, mean(CornerR,3)*squeeze(-1*d));
    angle4 = vectorAngle(forceVector, Rtest*squeeze(-1*d));    
    anglesRMSE = zeros(4,length(frames));
    forceTotal = zeros(3,length(frames));
    COP_error = zeros(3, length(frames));
    measurementLength = zeros(length(frames),1);
    anglestotalRMSE = [(sum(angle(f).^2)/length(f))^0.5, ...
        (sum(angle2(f).^2)/length(f))^0.5, ...
        (sum(angle3(f).^2)/length(f))^0.5, ...
        (sum(angle4(f).^2)/length(f))^0.5];
    for i = 1 : length(frames)
        fr = startFrames(frames(i)) : endFrames(frames(i));
        anglesRMSE(1,i) = (sum(angle(fr).^2)/length(fr))^0.5;
        anglesRMSE(2,i) = (sum(angle2(fr).^2)/length(fr))^0.5;
        anglesRMSE(3,i) = (sum(angle3(fr).^2)/length(fr))^0.5;
        anglesRMSE(4,i) = (sum(angle4(fr).^2)/length(fr))^0.5;
        COP_error(:,i) = mean(R*(squeeze(stickTip(:,:,fr)) - V') - forceCOP(:,fr),2);
        forceTotal(:,i) = mean(forceVector(:,fr),2);
        measurementLength(i) = length(fr);
    end
    
    errorData = struct('angles',anglesRMSE,'anglesTotal',anglestotalRMSE,'length',measurementLength,'COP_error', COP_error, 'force',forceTotal,'meanFpCOP',meanFpCOP(:,frames),'COP_RMSE',COP_RMSE);
    save(strcat(path,file,'errorData.mat'),'errorData')
    fprintf('Angle\n')
    fprintf('RMSE of angle using fminsearch: %9.6f%c\n', anglestotalRMSE(1), char(176));
    fprintf('RMSE of angle using fmincon: %9.6f%c\n', anglestotalRMSE(2), char(176));
    fprintf('RMSE of angle using corner calibration: %9.6f%c\n', anglestotalRMSE(3), char(176));
    fprintf('RMSE of angle using challis: %9.6f%c\n', anglestotalRMSE(4), char(176));
    clear fr
end


%% Z-score - compare forces
% Performs a z - test, to check if the size of the force has an influence
% on the angle between the force vector and the stick vector.
force_norm = zeros(length(forceVector),1);
validFrames = zeros(length(angle),1);
for i = 1 : length(startFrames)
    validFrames(startFrames(i):endFrames(i)) = 1;
end
validFrames = (validFrames == 1);
for i = 1:length(forceVector)
    force_norm(i) = norm(forceVector(:,i));
end
force_norm_greater = force_norm > mean(force_norm);
angle_greater = angle(and(force_norm_greater,validFrames));
angle_less = angle(and(not(force_norm_greater),validFrames));
fprintf('Average of angle where force > average force: %f, sd: %f\n', mean(angle_greater), sqrt(var(angle_greater)));
fprintf('Average of angle where force <= average force: %f, sd: %f\n', mean(angle_less), sqrt(var(angle_less)));
zscore = sqrt(sum(endFrames - startFrames))*norm(mean(angle_greater) - mean(angle_less))/max(sqrt(var(angle_greater)), sqrt(var(angle_less)));
fprintf('Z score: %f\n', zscore);
%% Sphere fit accuracy
totalVariance = [];
totalLength = [];
totalCOP = [];
totalError = [];
totalI = [];
for i = 1 : length(frames)
    len = endFrames(frames(i)) - startFrames(frames(i));
    lenSphereFit = 100 : 100 : len;
    for j = 1 : length(lenSphereFit)
        [COP2 , variance, ~] = Sphere_fit_final(strcat(path, file), startFrames(frames(i)), startFrames(frames(i)) + lenSphereFit(j), identifiedMarkers);
        totalVariance = [totalVariance, variance];
        totalLength = [totalLength, lenSphereFit(j)];
        totalCOP = [totalCOP, COP2];
        totalError = [totalError, (map(i,:) - transform(COP2',mean(CornerR,3), mean(CornerV,3)))'];
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

%% Compare sphere fitting and rigid body transformation
% For each of the frames where the tip of the stick is held still, output
% the average difference between sphere fitting COP and the tip of the
% stick calculated by the rigid body transformation
COPSphereFit = [];
varianceSphereFit = [];
error = [];
measurementFrames = [];
meanErrorSphereFit = [];
varianceSphereFit2 = [];
for i=1:length(startFrames)
    [COP2, variance, centers] = Sphere_fit_final(strcat(path, file), startFrames(i), endFrames(i), identifiedMarkers);
    measurementFrames = [measurementFrames startFrames(i):endFrames(i)];
    COPSphereFit = [COPSphereFit COP2];
    varianceSphereFit = [varianceSphereFit variance];
    error = [error squeeze(stickTip(:,:,startFrames(i):endFrames(i)) - COP2)];
    mean_error = mean(squeeze(stickTip(:,:,startFrames(i):endFrames(i)) - COP2),2);
    meanErrorSphereFit = [meanErrorSphereFit mean_error];
    varianceSphereFit2 = [varianceSphereFit2 var(squeeze(stickTip(:,:,startFrames(i):endFrames(i)) - COP2),0,2)];
    
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

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
map = errorMap();

%% Calculate force vector components - PILS
if (~worldCoordinates)
    
    %Calculate magnitude of force vector
    forceMagnitude = vecnorm(squeeze(forceVector))';
    
    %Define reference force as stick vector (unit vector) * magnitude
    Fref = squeeze(d*-1).*forceMagnitude;
    
    %Rotate force (which is in global reference frame due to stick vector)
    %to local reference frame of force plate
    Fref =mean(CornerR,3) * Fref;
    
    %Define as stickTip in local coordinate system of force plate.
    COPf = squeeze(stickTip) - mean(CornerV,3)';
    
    %Rotate COP into local coordinate system
    COPf = mean(CornerR,3) * COPf;
    
    %Define moment as cross product of COPf [mm] and force
    Mref = cross(COPf/1000, Fref);
    
    %Define force plate data
    Ffp = decimateData(recording.(rec).Force.Force,recording.(rec).Force.SamplingFactor);
    Mfp = decimateData(recording.(rec).Force.Moment,recording.(rec).Force.SamplingFactor);
    
    %S = signal vector from force plate (in N). This could be wrong, it is
    %not clear what unit this should be in.
    % S = [Fx Fy Fz Mx My Mz] (force plate data)
    S = [Ffp; Mfp];
    
    %Stack into 12 x N matrix (see line 174 of insitucalibration.m)
    % RR = [Fx Fy Fz Mx My Mz] defined by stick vector & force magnitude
    RR = [Fref; Mref];
    
    %Calculate calibration matrix, which is RR / S (pseudoinverse) (line
    %177)
    %C = RR(:,f)'\S(:,f)';
    C = RR(:,f)/S(:,f);
    
    %Apply new calibration matrix to Ffp and Mfp to calculate new
    %forces/moments after calculating calibration matrix. (182)
    F1n = [Ffp; Mfp]'*C(:,1:3);
    M1n = [Ffp; Mfp]'*C(:,4:6);
    
    %COPZ = z component of force plate COP. 2.4 = thickness of plate.
    copz = mean(CornerV,3); copz = copz(3) - 2.4; copz = copz/1000;
    copz = 0;
    
    %Caclulated corrected COP (186)
    COP1nF = [ (copz*F1n(:,1)-M1n(:,2))./F1n(:,3)     (M1n(:,1)+copz*F1n(:,2))./F1n(:,3)     copz*ones(length(F1n),1) ];
    
    COP1nF = COP1nF * 1000;
    RErrorMap = mean(CornerR,3);
    VErrorMap = mean(CornerV,3);
    [meanCorrectedFpCOP, varianceCorrectedFp] = getCOP(COP1nF',startFrames, endFrames);
    transformedStickCOP = (RErrorMap*(meanStickCOP' - VErrorMap)')';
    sphereFit = (RErrorMap*(COPSphereFit' - VErrorMap)')';
    
    variancePlot(map, meanCorrectedFpCOP(1:2,frames)',varianceCorrectedFp(1:2,frames)',...
        transformedStickCOP(frames,1:2),varianceStick(1:2,frames)',...
        sphereFit(frames,1:2),varianceSphereFit(1:2,frames)');
    title('Error Map with applied calibration matrix')
    set(gcf,'name','CalibrationErrorMap')
end
%% Error map plots - using QTMs COP.
% map = [[500/2, -298.5/2, 0];
%     [500/2,  298.5/2, 0];
%     [-500/2, 298.5/2, 0];
%     [-500/2, -298.5/2, 0]];
% map = map - [[88.0, -51.2, 0];
%      [96.6, 70, 0];
%     [-89.1, 41.9,0];
%     [-92.1, -41.0, 0]];


%Transforms the force plate COP to the reference frame using the rotation
%and translation matrices defined by the corner calibration.
transformedFpCOP = meanFpCOP;
RErrorMap = mean(CornerR,3);
VErrorMap = mean(CornerV,3);

if(worldCoordinates)
    transformedFpCOP = transform(meanFpCOP', RErrorMap, VErrorMap);
else
    transformedFpCOP = meanFpCOP';
end
transformedStickCOP = transform(meanStickCOP', RErrorMap, VErrorMap);
sphereFit = transform(COPSphereFit', RErrorMap, VErrorMap);

variancePlot(map, transformedFpCOP(frames,1:2),varianceFp(1:2,frames)',...
    transformedStickCOP(frames,1:2),varianceStick(1:2,frames)',...
    sphereFit(frames,1:2),varianceSphereFit(1:2,frames)');
set(gcf,'name','NormalErrorMap')

%% Make error map with distance from COP as arrow
% arrowPlot(map, transformedFpCOP(frames,1:2))

distances = [vecnorm((transformedStickCOP(frames,:) - map)'),...
    vecnorm((transformedFpCOP(frames,:) - map)'),...
    vecnorm((sphereFit(frames,:) - map)')];

errorMapDistances = struct('NormalizedDistances',distances,...
    'stick',transformedStickCOP(frames,:) - map,...
    'sphereFit',sphereFit(frames,:) - map,...
    'force_plate',transformedFpCOP(frames,:) - map,...
    'var_stick',varianceStick(:,frames)',...
    'var_spherefit', varianceSphereFit(:,frames)',...
    'var_fp', varianceFp(:,frames)');

save(strcat(path,'/',file,'errorMapDistances.mat'),'errorMapDistances')
%% Force plate manual calculation using analog data
data = recording.(rec).Analog.Data;
channels = 17:24;
data = data(channels,:)';
data = data(1:10:length(data),:);
min = 10e10;

for i = -100:0
    [analogFpCOP,analogFpForce] = analogDataAnalysis(data, i);
    q = sum(vecnorm((analogFpCOP(:,f)' - transform(squeeze(stickTip(:,:,f))',RErrorMap, VErrorMap))'));
    if (q < min)
        min = q
        minIndex = i
    end
end
clear q;
[analogFpCOP,analogFpForce] = analogDataAnalysis(data, minIndex);

[meanAnalogFpCOP, varianceAnalogFpCOP] = getCOP(analogFpCOP, startFrames, endFrames);
[meanAnalogFpForce, varianceAnalogFpForce] = getCOP(analogFpForce, startFrames, endFrames);
variancePlot(map,meanAnalogFpCOP(1:2,frames)',varianceAnalogFpCOP(1:2,frames)',...
    transformedStickCOP(frames,1:2),varianceStick(1:2,frames)',...
    sphereFit(frames,1:2),varianceSphereFit(1:2,frames)')
title('Analog data error map')
set(gcf,'name','Analog error map')

filterfreq = 25;  samplefreq = 1000;  [w,q] = butter(3,filterfreq/(samplefreq/2)); newanalogforceCOP = filtfilt(w,q,analogFpCOP')';
filterfreq = 25;  samplefreq = 1000;  [w,q] = butter(3,filterfreq/(samplefreq/2)); newForce = filtfilt(w,q,analogFpForce')';

figure('Name','analog_filtered_cop_fp_stick')
hold on;
plot(1:length(newanalogforceCOP(:,f)), newanalogforceCOP(:,f));

R = mean(CornerR,3);
V = mean(CornerV,3);
Q = squeeze(stickTip); Q = R * (Q - V');
plot(1:length(Q(:,f)), Q(1:2,f));
legend('COP X','COP Y', 'COP Z', 'Stick X', 'Stick Y')
ylabel('Coordinates [mm]')
xlabel('Frame number')
title(strcat('Location of FILTERED COP in force plate reference frame h = %dmm', minIndex));

%% Save all figures
% https://nl.mathworks.com/matlabcentral/answers/182574-save-all-the-plots
FolderName = strcat(path,'figures\',file ,'\');   % Your destination folder
mkdir(FolderName)
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
    FigHandle = FigList(iFig);
    FigName   = get(FigHandle, 'Name');
    if isempty(FigName)
        FigName = int2str(iFig);
    end
    savefig(FigHandle, strcat(FolderName,FigName,'.fig'));
end
%% 
close all; save(strcat(path, 'allData',file));