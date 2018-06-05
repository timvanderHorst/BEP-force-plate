function [] = createPlots(worldCoordinates,  forceVector, forceCOP, d, stickTip, identifiedMarkers, startFrames, endFrames, f, R, V)
%CREATE_PLOTS - Create plots using previously calculated 
%Segments a recording into shorter measurements by analyzing the movement
%of the tip of the pole. If this movement is below a threshold, it is part
%of a measurement
%
% Syntax:  createPlots(worldCoordinates,  forceVector, forceCOP, d, stickTip, identifiedMarkers, startFrames, endFrames, f)
%
% Inputs:
%    worldCoordinates - logical vector, 1 if using world coordinates
%    forceVector - 3 X N vector containing force plate force vector
%    forceCOP - 3 x N vector containing force plate COP coordinates
%    d - 3 x 1 x N vector containing vector along pole's central axis
%    stickTip - 3 x 1 x N vector of position of pole tip in global frame
%    identifiedMarkers - M x 4 x N matrix of positions of markers as
%       measured by camera's 
%    startFrames - 1 x F vector with start frame numbers
%    endFrames - 1 x F vector with end frame numbers
%    f - 1 x X vector with frame numbers
%    R - rotation matrix used to transform sphere fit location to force
%       plate reference frame
%    V - translation vector used to transform sphere fit location to force
%       plate reference frame

%If using world coordinates, force = opposite direction force

if(worldCoordinates)
    forceVector = forceVector*-1;
else
    
end
angle = vectorAngle(forceVector, d);
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
filterfreq = 25;  samplefreq = 1000;  [w,q] = butter(3,filterfreq/(samplefreq/2)); newForce = filtfilt(w,q,forceVector')';

figure('Name','filtered_cop_fp_stick')
hold on;
plot(1:length(newforceCOP(:,f)), newforceCOP(:,f));

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
xlabel('Frame numbeR')

figure('Name','unfiltered_stick_tip')
hold on;
plot(1:length(stickTip), squeeze(stickTip));
legend('X','Y', 'Z')
title('Location of tip of stick in global reference frame');
ylabel('Coordinates [mm]')
xlabel('Frame number')

figure('Name','unfiltered_force_fp')
hold on;
plot(1:length(forceVector), forceVector);
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


end

