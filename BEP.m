%Force plate analysis - 15/03/2018
close all;

%% IMPORTANT VARIABLES
plate = 'c';    %plate name
MeasurementNumber = 1;
if (lower(plate) == 'c')
    map = error_map_locations();%positions on plate where forces are exerted
else
    map = error_map_locations(1);
end
NOMAP = true;
%% 
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
%Calculate tip of pole location using pole geometry
vGeometric = geometric_tip(identifiedMarkers, markerNumbers);

%% Veldpaus, 1988 and Rose
%Perform rigid body transformation using method described by Spoor, 
%Veldpaus (1980). Tries to calculate a rotation matrix and translation
%vector that maps the markers in the reference frame defined above (
% variable a) to the markers in the %global coordinate system. 
% This function returns a Rotation matrix and translation vector for every 
% frame.
[RVeldpaus_1, vVeldpaus_1, valid] = spoorVeldpausRotation(a, identifiedMarkers);

%Calculate the RMS, by mapping the markers in the global frame back onto
%the reference frame, and calculating the rms of the total distance between
%the mapped markers and the markers in the reference frame.
[EVeldpaus_1] = rms(a, identifiedMarkers, RVeldpaus_1, vVeldpaus_1, valid);
fprintf('RMSE Spoor, Veldpaus eigenvector decomposition: %d\n', sum(EVeldpaus_1 .* valid) / sum(valid));

%Repeat the same procedure using SVD instead of eigenvalue decomposition
[RVeldpaus_2,vVeldpaus_2] = spoorVeldpausRotationSVD(a, identifiedMarkers);
[EVeldpaus_2] = rms(a, identifiedMarkers, RVeldpaus_2, vVeldpaus_2, valid);
fprintf('RMSE using Spoor, Veldpaus SVD: %d\n', sum(EVeldpaus_2 .* valid) / sum(valid));

%% Challis, 1995
%Perform rigid body transformation using method described by Challis (1995)
[RChallis,vChallis, valid] = challis(a, identifiedMarkers);
[EChallis] = rms(a, identifiedMarkers, RChallis, vChallis, valid);
fprintf('RMSE using Challis SVD: %d\n', sum(EChallis .* valid) / sum(valid));

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

%Remove frames from start and end of measurement
for i = 1 : length(startFrames)
    startFrames(i) = startFrames(i) + 120;
    endFrames(i) = endFrames(i) - 120;
end
%% Force vector
%Resample vectors to match camera sampling rate.
forceVector = decimateData(recording.(rec).Force.Force,recording.(rec).Force.SamplingFactor);
forceCOP = decimateData(recording.(rec).Force.COP,recording.(rec).Force.SamplingFactor);
%% Save measurement frames variables to memory.
    frames = [];
    THRESHOLD_Z_FORCE = 15; %Min. force during measurement
    THRESHOLD_LENGTH = 20;  %0.2s measurement
    THRESHOLD_Z = 200;      %Max. distance of pole tip from floor
    
    %for each measurement (startFrame to endFrame), check if the force is
    %high enough, the measurement is long enough and if the tip of the pole
    %is below 200mm from the floor of the lab.
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
    
    %Save these frames to the file
    recording.(rec).MeasurementFrames = frames;
    fprintf('Saving frames to file\n')
    data = recording.(rec);
    save(strcat(path,file), 'data')
    fprintf('Done saving \n')
    
    clear f p x string THRESHOLD_Z_FORCE THRESHOLD_LENGTH data
    
    %Define variable f as the frames that were actually part of a
    %measurement
    f = [];
    for i = 1 : length(frames)
        f = [f startFrames(frames(i)):endFrames(frames(i))];
    end


%% Force plate location
%Calculate the positions of the tip of the stick and the COP of the force
%plate during the frames which were calculated above.
[meanStickCOP, varianceStick] = getCOP(stickTip, startFrames, endFrames);
[meanFpCOP, varianceFp] = getCOP(forceCOP,startFrames, endFrames);


% Use challis to calculate a rotation and translation matrix which would
% map the coordinate system of the force plate onto the global coordinate
% system. 
% This is done by entering the FpCOP values as the reference frame positions and the tip of the
% pole location as the global position. Uses the locations in each of the
% frames specified above to perform this calculation.
[R_fp_challis,V_fp_challis, validtest] = challis(meanFpCOP(:,frames), [meanStickCOP(:,frames)' zeros(size(meanStickCOP(:,frames),2),1)]);
% Calculate the rotation matrix and translation vector of the corners, by
% asking the user to input the files they would like to analyze, and
% performing the procedure described in Rabuffetti, 2012
[CornerMean, CornerVar,positions, CornerR, CornerV, meanCornerR, meanCornerV] = cornerCalibration(path, 9);



%If no error map is used, define the map as the average pole tip positions.
if(NOMAP)
    map = transform(meanStickCOP', RErrorMap, VErrorMap);
    map = map(frames,:);
end

%% Rabuffetti, 2002

Rstart = R_fp_challis;
Vstart = V_fp_challis';

% T = transformation vector
T = rabuffetti( d(:,:,f), ...
    stickTip(:,:,f), ...
    -1*forceVector(:,f), ...
    forceCOP(:,f), ...
    Rstart,Vstart');
%Transform transformation parameters from euler angles, type ZYZ back to a
%rotation vector.
R1final = eul2rotm(T(1:3)','ZYZ');
V1final = T(end-2:end);

fprintf('\n\n                         fminsearch\n');
printOptimisationResults(Rstart,R1final,Vstart,V1final);
clear Rstart Vstart

%% Plot figures
createPlots(worldCoordinates, forceVector, forceCOP, d, stickTip, identifiedMarkers, startFrames, endFrames, f, meanCornerR, meanCornerV)

%% Calculate error between stick tip and 
if(~worldCoordinates)
    errorData = total_rmse_calc( strcat(path,file,'normal'), f, frames, startFrames, endFrames, stickTip, d, forceCOP, forceVector, R1final, V1final, CornerR, CornerV, R_fp_challis, V_fp_challis);
end

%% Z-score - compare forces
% Performs a z - test, to check if the size of the force has an influence
% on the angle between the force vector and the stick vector.
z_score = zscore( forceVector, d, f );
%% Compare sphere fitting and rigid body transformation
% For each of the frames where the tip of the stick is held still, output
% the average difference between sphere fitting COP and the tip of the
% stick calculated by the rigid body transformation
 [COPSphereFit, varianceSphereFit, meanErrorSphereFit ] = sphere_fit_vs_pole_tip( startFrames, endFrames, stickTip, identifiedMarkers );

%% Sphere fit accuracy
%Create a scatter plot of length of length of necessary sphere fit vs
%distance between stick tip and sphere fit
sphere_fit_scatter(startFrames, endFrames, frames, meanCornerR, meanCornerV, identifiedMarkers, map);
%% Calculate calibration matrix - PILS
if (~worldCoordinates)
    momentVector = decimateData(recording.(rec).Force.Moment,recording.(rec).Force.SamplingFactor);
    [C, COP1nF, F1n, M1n] = calibration_matrix( f, forceVector,momentVector ,d, stickTip, meanCornerR, meanCornerV');
    
    RErrorMap = meanCornerR;
    VErrorMap = meanCornerV;
    [meanCorrectedFpCOP, varianceCorrectedFp] = getCOP(COP1nF',startFrames, endFrames);
    transformedStickCOP = (RErrorMap*(meanStickCOP' - VErrorMap)')';
    sphereFit = (RErrorMap*(COPSphereFit' - VErrorMap)')';
    
    error_map_plot(NOMAP, map, meanCorrectedFpCOP(1:2,frames)',varianceCorrectedFp(1:2,frames)',...
        transformedStickCOP(frames,1:2),varianceStick(1:2,frames)',...
        sphereFit(frames,1:2),varianceSphereFit(1:2,frames)');
    title('Error Map with applied calibration matrix')
    set(gcf,'name','CalibrationErrorMap')
end
%% Error map plots - using QTMs COP.
%Transforms the force plate COP to the reference frame using the rotation
%and translation matrices defined by the corner calibration.
transformedFpCOP = meanFpCOP;
RErrorMap = meanCornerR;
VErrorMap = meanCornerV;

if(worldCoordinates)
    transformedFpCOP = transform(meanFpCOP', RErrorMap, VErrorMap);
else
    transformedFpCOP = meanFpCOP';
end
transformedStickCOP = transform(meanStickCOP', RErrorMap, VErrorMap);
sphereFit = transform(COPSphereFit', RErrorMap, VErrorMap);

error_map_plot(NOMAP, map, transformedFpCOP(frames,1:2),varianceFp(1:2,frames)',...
    transformedStickCOP(frames,1:2),varianceStick(1:2,frames)',...
    sphereFit(frames,1:2),varianceSphereFit(1:2,frames)');
set(gcf,'name','NormalErrorMap')

%% Make error map with distance from COP as arrow
arrowPlot(map, transformedFpCOP(frames,1:2))
if(not(NOMAP))
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
end
%% Force plate manual calculation using analog data
data = recording.(rec).Analog.Data;

if(lower(plate) == 'c')
    channels = 17:24;
elseif(lower(plate) == 'a')
    channels = 1:8;
end
data = data(channels,:)';
data = data(1:10:length(data),:);

%Find best h-offset value by minimizing force plate - stick tip distance
minH = height_correction(plate, data, stickTip, meanCornerR, meanCornerV, f);
[analogFpCOP,analogFpForce] = analogDataAnalysis(plate, data, minH);

%Plot values using optimized h-offset in an error map
[meanAnalogFpCOP, varianceAnalogFpCOP] = getCOP(analogFpCOP, startFrames, endFrames);
[meanAnalogFpForce, varianceAnalogFpForce] = getCOP(analogFpForce, startFrames, endFrames);
error_map_plot(NOMAP, map,meanAnalogFpCOP(1:2,frames)',varianceAnalogFpCOP(1:2,frames)',...
    transformedStickCOP(frames,1:2),varianceStick(1:2,frames)',...
    sphereFit(frames,1:2),varianceSphereFit(1:2,frames)')
title('Analog data error map')
set(gcf,'name','Analog error map')

%Plot values using optimized h-offset in a filtered COP graph
filterfreq = 25;  samplefreq = 1000;  [w,q] = butter(3,filterfreq/(samplefreq/2)); newanalogforceCOP = filtfilt(w,q,analogFpCOP')';
filterfreq = 25;  samplefreq = 1000;  [w,q] = butter(3,filterfreq/(samplefreq/2)); newForce = filtfilt(w,q,analogFpForce')';

figure('Name','analog_filtered_cop_fp_stick')
hold on;
plot(1:length(newanalogforceCOP(:,f)), newanalogforceCOP(:,f));

R = meanCornerR;
V = meanCornerV;
Q = squeeze(stickTip); Q = R * (Q - V');
plot(1:length(Q(:,f)), Q(1:2,f));
legend('COP X','COP Y', 'COP Z', 'Stick X', 'Stick Y')
ylabel('Coordinates [mm]')
xlabel('Frame number')
title(strcat('Location of FILTERED COP in force plate reference frame h = %dmm', minH));
clear R V Q
%% Save all figures
% https://nl.mathworks.com/matlabcentral/answers/182574-save-all-the-plots
FolderName = strcat(path,'figures\',file ,'\');   % Your destination folder
mkdir(FolderName)
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
    FigHandle = FigList(iFig);
    FigName   = get(FigHandle, 'Name');
    if (length(FigHandle.Children) == 1)
        i = 1;
    else 
        i = 2;
    end
    OldTitle = FigHandle.Children(i).Title.String;
    FigHandle.Children(i).Title.String = strcat(sprintf('M%d: ', MeasurementNumber), OldTitle);
    if isempty(FigName)
        FigName = int2str(iFig);
    end
    savefig(FigHandle, strcat(FolderName,FigName,'.fig'));
    saveas(FigHandle, strcat(FolderName,FigName,'.png'));
    FigHandle.Children(i).Title.String = OldTitle;
end
%% 
save(strcat(path, 'allData',file));