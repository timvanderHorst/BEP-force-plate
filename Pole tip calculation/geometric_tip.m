function v = geometric_tip( markers, markerNumbers )
%GEOMETRIC_TIP - Calculates the location of the tip of an instrumented pole
%using marker geometry.
%Calculates the point a certain distance below the average position of the
%bottom two markers, along the central axis of the pole
%
% Syntax:  [v] = function_name(markers, markerNumbers)
%
% Inputs:
%    markers - 5 x 4 x N matrix containing [x y z residual] data for each
%       marker, where N is the number of frames
%    markerNumbers - a vector containing the ordering of the markers in the
%    markers matrix
%
% Outputs:
%    v - 3 x N matrix containing location of the tip of the pole for all N
%    frames
%
% The following marker ordering is used: 
% m1 = [0; -160; 200];
% m2 = [0; 160; 200];
% m3 = [0; 160; 400];
% m4 = [0; -160; 800];
% m5 = [0; 160; 800];


%Split markers matrix into top and bottom markers
topMarkers = markers([find(markerNumbers == 5), find(markerNumbers == 4)],:,:);
bottomMarkers = markers([find(markerNumbers == 1), find(markerNumbers == 2)],:,:);

%Calculate mean position of top and bottom markers
meanTop = squeeze(mean(topMarkers,1));
meanBottom = squeeze(mean(bottomMarkers,1));

%Calculate vector along pole's central axis
d2 = meanTop(1:3,:) - meanBottom(1:3,:);

%Normalize vector
d2 = d2./vecnorm(d2)';

%Define point as [0; 0; 0] in pole reference frame, which is 200mm below
%mean position of bottom two markers.
v = meanBottom(1:3,:) - 200*d2;
end

