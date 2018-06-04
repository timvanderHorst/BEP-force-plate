close all
clc

maxDistance = 2000;

%Define starting frame and ending frame
startFrame = 6613; %1200; %6613
endFrame = 11000; %5400; %11000

%Define locations of markers
m1 = [0; -160; 200]
m2 = [0; 160; 200]
m3 = [0; 160; 400]
m4 = [0; -160; 800]
m5 = [0; 160; 800]
a = [m1 m2 m3 m4 m5];

%Load file
recording = load('././2018-4-23-BEP-forceplate/Data/recording.mat');
rec = fieldnames(recording);
rec = rec{1};
identifiedMarkers = recording.(rec).Trajectories.Labeled.Data;

distances = zeros(1,5);
distances2 = distances;
a2 = identifiedMarkers(:,1:3,1)';
markerNumbers = zeros(1,5);
for i = 1:5
    for j = 1:5
        distances(i) = distances(i) + norm(a(:,i)-a(:,j));
        distances2(i) = distances2(i) + norm(a2(:,i)-a2(:,j));
    end
end
for i = 1:5
    [c index] = min(abs(distances-distances2(i)));
    markerNumbers(i) = index;
end



markers = zeros(size(identifiedMarkers,1),endFrame - startFrame,3);

%Reshape array Nx4xM array (N = Number of markers, M = Number of
%measurements) to an N x M x 3 matrix
for i = 1 : endFrame - startFrame
    markers(:,i,:) = identifiedMarkers(:,1:3,i + startFrame);
end

%Apply sphere fitting to each marker
centers1 = zeros(size(identifiedMarkers,1),3);
radius1 = zeros(size(identifiedMarkers,1),1);
meanError1 = zeros(size(identifiedMarkers,1),1);

centers2 = zeros(size(identifiedMarkers,1),3);
radius2 = zeros(size(identifiedMarkers,1),1);
meanError2 = zeros(size(identifiedMarkers,1),1);

centers3 = zeros(size(identifiedMarkers,1),3);
radius3 = zeros(size(identifiedMarkers,1),1);
meanError3 = zeros(size(identifiedMarkers,1),1);

%METHOD 1: Matlab point cloud sphere fit
for i = 1 : size(identifiedMarkers,1)
    ptCloud = pointCloud(squeeze(markers(i,:,:)));
    center = pcfitsphere(ptCloud,maxDistance);
    %figure;
    %pcshow(ptCloud)
    centers1(i,:) = center.Center;
    radius1(i) = center.Radius;
    meanError1(i,:) = sphere_fit_error(centers1(i,:), radius1(i), squeeze(markers(i,:,:)))/size(markers,2);

end
fprintf('Method 1 mean:     %.2f, %.2f, %.2f\n', mean(centers1))
fprintf('Method 1 variance: %.2f, %.2f, %.2f\n', var(centers1))

%METHOD 2: 
for i = 1 : size(identifiedMarkers,1)
    [centers2(i,:), radius2(i)] = sphereFit(squeeze(markers(i,:,:)));
    meanError2(i,:) = sphere_fit_error(centers2(i,:), radius2(i), squeeze(markers(i,:,:)))/size(markers,2);
end
fprintf('Method 2 mean:     %.2f, %.2f, %.2f\n', mean(centers2))
fprintf('Method 2 variance: %.2f, %.2f, %.2f\n', var(centers2))


%METHOD 3: 
for i = 1 : size(identifiedMarkers,1)
    [centers3(i,:), radius3(i)] = sphere_fit(squeeze(markers(i,:,:)));
    meanError3(i,:) = sphere_fit_error(centers3(i,:), radius3(i), squeeze(markers(i,:,:)))/size(markers,2);
end
fprintf('Method 3 mean:     %.2f, %.2f, %.2f\n', mean(centers3))
fprintf('Method 2 variance: %.2f, %.2f, %.2f\n', var(centers3))


