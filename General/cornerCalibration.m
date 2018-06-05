function [Mean,STD, positions, R, v, meanR, meanV] = cornerCalibration(varargin)
%CORNERCALIBRATION - Calculates the positions of the corners of the
%force plate.
%The script will ask the user to select a file where the recording of the
%corners is stored. The markers must be numbered (in a clockwise direction
%around the force plate, starting in the +X, +Y corner - 3, 2, 1, 4)
%
% Syntax:  [~, ~, ~, R, V ] = cornerCalibration('C:/')
% Syntax:  [~, ~, ~, R, V ] = cornerCalibration('C:/', 9)

% Inputs:
%    path - the path to the base directory where the file selection window
%    will be opened
%    markerOffset - the offset in mm from the top of the force plate to the
%    center of the markers. If nothing is entered, 7mm is used by default.
% Outputs: 
%    Mean = 4 x 3 matrix containing the positions of each corner
%    STD = 4 x 3 matrix containing standard deviation of each coordinate,
%       calculated using the data for all markers
%    positions = 4 x 3 x N matrix, where N is the number of recordings,
%       containing the positions of the corners in each recording
%    R = 3 x 3 x N matrix containing the rotation matrix that maps the global
%       reference frame to the force plate coordinate system
%    v = 1 x 3 x N matrix containing the translation vector that contains the
%       position of the force plate coordinate system, expressed in the global
%       coordinate system.
markerHeight = 7;
path = varargin{1};
if(nargin == 2)
    markerHeight = varargin{2}
end


%Load files
if ((exist('recording') + exist('file') + exist('path')) ~= 3)
    [file,path] = uigetfile(strcat(path,'*.mat'),'Select input file(s)','MultiSelect','on');
    recording = {};
    rec = {};
    if(~iscell(file))
        file = {file};
    end
    for i = 1:length(file)
        recording{i} = load(strcat(path,file{i}));
        r = fieldnames(recording{i});
        rec{i} = r{1};
    end
else
    x = input('Would you like to use the current file: ? Y/N [Y]: ','s');
    if(~isempty(x) && lower(x) == 'n')
        [file,path] = uigetfile('*.mat','Select input file(s)','MultiSelect','on');
        recording = {};
        rec = {};
        for i = 1:length(file)
            recording{i} = load(strcat(path,file{i}));
            r = fieldnames(recording{i});
            rec{i} = r{1};
        end
    end
end

meanPositions = zeros(4,3,length(file));
center = zeros(3,length(file));
stdPositions = zeros(4,4,length(file));
for i = 1:length(file)
    identifiedMarkers = recording{i}.(rec{i}).Trajectories.Labeled.Data;
    meanPos = mean(identifiedMarkers,3);
    %Remove residual column
    meanPositions(:,:,i) = meanPos(:,1:3);
    
    %Center = geometric center of all markers
    center(:,i) = mean(meanPositions(:,:,i));
    
    stdPositions(:,:,i) = std(identifiedMarkers,0,3);
end


% For each of the markers on the corners, and for each file compare the
% ordering of the markers in the meanPositions matrix to the ordering in
% the first file. Fill a matrix which indicates the row number of the
% marker in that file.
markers = zeros(4,length(file));
for i = 1 : 4
    markers(i,1) = i;
    for j = 2 : length(file)
        minDistance = 1000;
        minIndex = 20;
        for k = 1 : 4
            distance = meanPositions(k,:,j) - meanPositions(i,:,1);
            distance = norm(distance);
            
            %The marker with the smallest distance to the marker in the
            %first file is assumed to be the same marker. Specify the index
            %of that marker in the markers matrix.
            if distance < minDistance
                markers(i,j) = k;
                minDistance = distance;
            end
        end
    end
end

oldPos = meanPositions;
%Reorder the mean positions matrix so the markers are ordered correctly.
for j = 2 : length(file)
    meanPositions(:,:,j) = meanPositions(markers(:,j),:,j);
end

%Get the rotation and translation matrix which maps the local coordinates
%onto the global coordinates for each file.
[R,v] = getTransformation(recording, rec, meanPositions);
v = v - [0 0 markerHeight];
Mean = mean(meanPositions,3) - [0 0 markerHeight];
STD = std(meanPositions,0,3);
positions = meanPositions - [0 0 markerHeight];
[meanR, meanV]  = getTransformation(recording, rec, mean(meanPositions,3));

    function [x, y, z] = getRotation(x1, x2, y1, y2)
        %vector parallel to x,y axis = positive side - negative side
        xAxis = x1 - x2;
        yAxis = y1 - y2;
        
        %perform Gram Schmidt to orthogonalize x axis and y axis
        yOrthogonal = yAxis - projVector(yAxis,xAxis);
        
        %Find z axis = cross(x,y)
        zAxis = cross(xAxis,yOrthogonal);
        x = xAxis./norm(xAxis);
        y = yOrthogonal./norm(yOrthogonal);
        z = zAxis./norm(zAxis);
    end

    function proj = projVector(v1, v2)
        %PROJVECTOR - projection of v1 onto v2
        proj = dot(v1, v2)./dot(v2,v2) .* v2;
    end

    function [R, v] = getTransformation(recording, rec, positions)
        %GETTRANSFORMATION - Calculates a rotation matrix and translation vector
        %corresponding with the coordinates of the 4 corners of the force plate.
        %Marker 4,1 = edge parallel with x-axis (positive x side)
        %Marker 3,4 = edge parallel with y-axis (positive y side)
        %INPUT:
        %   recording = recording cell, with N recording structures
        %   rec = names of recordings
        %   positions = mean positions of corners, 4 x 3 x N matrix
        %OUTPUT:
        %   R = rotation matrix for each recording (3 x 3 x N)
        %   v = translation vector, 1 x 3 x N
        
        %translation vector = geometric mean of coordinates of corners
        v = mean(positions);
        labels = recording{1}.(rec{1}).Trajectories.Labeled.Labels;
        labels = cell2mat(labels');
        
        %define marker 1,2,3,4 coordinates in each recording
        m1 = positions(find(labels == '1'),:,:);
        m2 = positions(find(labels == '2'),:,:);
        m3 = positions(find(labels == '3'),:,:);
        m4 = positions(find(labels == '4'),:,:);
        
        %center of edge on positive x side, parallel to x axis = average
        %coordinates of markers 4,1
        x1 = mean([m4;m1]);
        
        %center of edge on negative x side, parallel to x axis = average
        %coordinates of markers 3,2
        x2 = mean([m3;m2]);
        
        %center of edge on positive y side, parallel to y axis = average
        %coordinates of markers 3,4
        y1 = mean([m3;m4]);
        
        %center of edge on negative y side, parallel to y axis = average
        %coordinates of markers 2,1
        y2 = mean([m2;m1]);
        
        %vector parallel to x,y axis = positive side - negative side
        xAxis = x1 - x2;
        yAxis = y1 - y2;
        
        %perform Gram Schmidt to orthogonalize x axis and y axis
        yOrthogonal = yAxis - projVector(yAxis,xAxis);
        
        %Find z axis = cross(x,y)
        zAxis = cross(xAxis,yOrthogonal);
        if(size(xAxis,3) == 1)
            x = xAxis./norm(xAxis);
            y = yOrthogonal./norm(yOrthogonal);
            z = zAxis./norm(zAxis);
            R = [x;y;z];
        else
            x = squeeze(xAxis)'./vecnorm(squeeze(xAxis));
            y = squeeze(yOrthogonal)'./vecnorm(squeeze(yOrthogonal));
            z = squeeze(zAxis)'./vecnorm(squeeze(zAxis));
            %Rotation matrix is defined as unit vectors of force plate coordinate
            %system, expressed in global coordinate system R = [ex; ey; ez]
            R = zeros(3,3,length(x));
            for i = 1 : length(x)
                R(:,:,i) = [x(i,:);y(i,:);z(i,:)];
            end
        end
    end
end

