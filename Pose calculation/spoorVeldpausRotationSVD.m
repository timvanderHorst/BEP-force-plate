function [ R, v, valid ] = spoorVeldpausRotationSVD(a, markerData)
%SPOOR_VELDPAUS_ROTATION Uses the algorithm described in Spoor Veldpaus
%(1980) and Rose (2012) to calculate the transformation parameters
%Calculates transformation parameters that maps vectors defined in a
%reference frame onto positions defined in a different frame. Uses singular
%value decomposition.
%
% Syntax:  [R, v, valid] = spoorVeldpausRotationSVD(a, markerData)
%
% Inputs:
%    a - 3 x M matrix containing 3D positions of M markers in reference
%    frame
%    markerData - 5 x 4 x N matrix containing [x y z residual] data for each
%       marker, where N is the number of frames
%
% Outputs:
%    v - 3 x 1 x N matrix containing translation vector that maps a onto
%    markerData
%    R - 3 x 3 x N matrix containing rotation matrix that maps a onto
%    markerData
%    valid - 1 x N vector containing 1 if frame is valid, 0 if frame is
%    invalid
%
% Rose (2012): https://www1.udel.edu/biology/rosewc/kaap427627/reserve/segment_tracking/Estimating%20Body%20Segment%20Motion%20by%20Tracking%20Markers%20R1.pdf

%Calculate geometric center of all markers in single vector (average of x, y, z columns)
aMean = mean(a')'; 

%p = position matrix of markers in global coordinate system
p = zeros(3, size(a,2), size(markerData,3)); 

% initialize empty mean position matrix and translation vector
pMean = zeros(3, 1, size(markerData,3)); v = pMean;

%Initialize empty matrices, eq. (5) of Rose (2012)
v = pMean;
M = zeros(3, 3, size(markerData,3));
V = M; %eigenvectors of M'M will be stored here
U = M; 
S = M; %eigenvalues of M
R = M; %rotation matrix stored here
valid = ones(size(markerData,3),1);

for i = 1 : size(markerData,3)
    %transpose marker position matrix and remove residual column
    p(:,:,i) = markerData(:,1:end-1,i)'; 
    
    %Calculate geometric mean of markers, store as column vector
    pMean(:,:,i) = mean(p(:,:,i)')';     
    
    %Define matrix M using eq. (5)
    M(:,:,i) = 1/size(a,2) * a * p(:,:,i)' - aMean*pMean(:,:,i)';
    
    %Check if marker values contain NaN values
    if sum(sum(isnan(M(:,:,i)))) == 0
        %SV decomposition of MM'
        [U(:,:,i), S(:,:,i), V(:,:,i)] = svd(M(:,:,i));
        MV = M(:,:,i) * V(:,:,i);
        v1 = V(:,1,i);
        v2 = V(:,2,i);
        v3 = V(:,3,i);
        s1 = S(1,1,i);
        s2 = S(2,2,i);
        %Use eq. (12) or (13) to calculate R, V
        if abs(sum(cross(v2,v3)./v1) - 3) < 1e-10 %v2 x v3 = +v1 
            R(:,:,i) = [MV(:,1)/s1 MV(:,2)/s2 cross(MV(:,1),MV(:,2))/(s1 * s2)]*V(:,:,i)';
            v(:,:,i) = pMean(:,:,i) - R(:,:,i)*aMean;
        elseif abs(sum(cross(v2,v3)./v1) + 3) < 1e-10 %v2 x v3 = -v1 
            R(:,:,i) = [MV(:,1)/s1 MV(:,2)/s2 cross(MV(:,1),MV(:,2))/(s1 * s2)]*V(:,:,i)';
            v(:,:,i) = pMean(:,:,i) - R(:,:,i)*aMean;
        else 
            valid(i) = 0;
        end
    else 
        valid(i) = 0;
    end
end
end

