function [ E, EMarkers ] = rms(a, markerData, R, v, valid)
%RMS Calculates the RMS of the rotation and translation vectors. 
%   Maps the markers in the global coordinate system back to the local coordinate
%   system using the transposed R and translation vector. Calculates the
%   sum of the difference between these coordinates and the reference
%   frame, returnrs the sqrt(sum)/size
%INPUT:
%   a - positions of markers in reference frame
%   markerData - positions of markers in global frame
%   R - rotation matrices
%   v - translation vectors
%   valid - matrix which specifies if this frame is valid.
%OUTPUT:
%   E - total RMS
%   EMarkers - RMS for each marker.

E = zeros(size(markerData,3),1);
EMarkers = zeros(size(a,2),3,size(markerData,3));
p = zeros(3, size(a,2), size(markerData,3)); %p = position matrix of markers in global coordinate system

for i = 1 : size(markerData,3)
    x = 0;
    p(:,:,i) = markerData(:,1:end-1,i)'; %transpose matrix and remove residual column
    if(valid(i) == 1)
        for j = 1 : size(a,2)
            EMarkers(j,:,i) = (R(:,:,i)'*p(:,j,i) - R(:,:,i)'*v(:,:,i)) - a(:,j);
            x = x + (R(:,:,i)*a(:,j) + v(:,:,i) - p(:,j,i))' * (R(:,:,i)*a(:,j) + v(:,:,i) - p(:,j,i));
        end
    end
    E(i) = sqrt(x/size(a,2));
end

end

