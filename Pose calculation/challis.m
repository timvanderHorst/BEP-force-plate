function [ R, v, valid ] = challis(a, markerData)
%CHALLIS Uses the algorithm described in Challis, 1995 to find the rotation
%matrix and translation matrix
%
%INPUT:
%   a - [N x 3] matrix containing marker positions in cartesian coordinates in reference
%frame
%   markerData - [N x 4 x T] matrix containing marker positions and residual 
%at each of the T frames
%OUTPUT: 
%   R - [3 x 3 x T] matrix containing 3D rotaiton matrix at each frame
%   v - [3 x 1 x T] matrix containing 3D translation vector at each frame
%   valid - [T x 1] vector which indicates if frame is valid or not (all
% markers have to be defined)

%Define empty matrices
aMean = mean(a,2);                                      %Calculate geometric center of all markers in single vector (average of x, y, z columns)
p = zeros(3, size(a,2), size(markerData,3));            %p = position matrix of markers in global coordinate system
pMean = zeros(3, 1, size(markerData,3));
v = pMean;
C = zeros(3, 3, size(markerData,3));
R = C;
valid = ones(size(markerData,3),1);

for i = 1 : size(markerData,3)
    p(:,:,i) = markerData(:,1:end-1,i)';                %transpose matrix and remove residual column
    pMean(:,:,i) = mean(p(:,:,i)')';                    %Calculate geometric mean of markers, store as column vector
    C(:,:,i) = (1/size(a,2))*(p(:,:,i) - pMean(:,:,i))*(a - aMean)'; %Calculate cross-dispersion matrix (19)
    if sum(sum(isnan(C(:,:,i)))) == 0
        [U,S,V] = svd(C(:,:,i));                        %Do singular value decomposition of C matrix (20)
        R(:,:,i) = U*[1 0 0; 0 1 0; 0 0 det(U*V')]*V';  %Calculate R using equation (24) in paper
        v(:,:,i) = pMean(:,:,i) - R(:,:,i)*aMean;       %Calculate translation vector using (28)
    else 
        valid(i) = 0;
    end
end

end

