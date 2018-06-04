function [ markerLocationMatrix, markerNumbers ] = getMarkerLocations( initialPosition, referenceFrameMarkers )
%MARKERS Finds the ordering of markers in the identifiedMarkers matrix,
%outputs their positions in the reference frame.
%   Based on an ideal ordering, finds the ordering in the measurement data
%   which corresponds to this ordering. For example, in the ideal case we
%   may define the matrix as [M1 M2 M3 M4 M5] while in the measurement,
%   these markers may be defined as [M3 M1 M2 M5 M4].
%   This is done by calculating the sum of distances between each marker to
%   each other marker, and comparing this sum to the sum of distances in
%   the ideal case.
%INPUT:
%   initialPosition: 3xN matrix containing marker positions in measurement
%   referenceFrameMarkers: 3xN matrix containing marker position in reference
%   frame of pole
%OUTPUT:
%   markerLocationMatrix: 3xN matrix containing positions of markers in
%   same order as in initialPosition matrix, but in the reference frame.

N = size(initialPosition,2);
if(N ~= size(referenceFrameMarkers,2))
    msgID = 'getMarkerLocations:invalidMarkerNumber';
    msg = 'Number of markers in reference frame and initial positions not equal';
    throw(MException(msgID,msg));
end

distances = zeros(1,N);
distances2 = zeros(1,N);
markerNumbers = zeros(1,N);
for i = 1:N
    for j = 1:N
        distances(i) = distances(i) + norm(referenceFrameMarkers(:,i)-referenceFrameMarkers(:,j));
        distances2(i) = distances2(i) + norm(initialPosition(:,i)-initialPosition(:,j));
    end
end

for i = 1:N
    [c index] = min(abs(distances-distances2(i)));
    markerNumbers(i) = index;
end
if(length(markerNumbers) ~= length(unique(markerNumbers)))
    msgID = 'getMarkerLocations:repeatedMarkers';
    msg = 'Could not determine ordering of markers';
    throw(MException(msgID,msg));
end
markerLocationMatrix = referenceFrameMarkers(:, markerNumbers);
end

