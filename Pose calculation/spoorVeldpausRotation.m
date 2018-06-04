function [ R, v, valid ] = spoorVeldpausRotation(a, markerData)
%SPOOR_VELDPAUS_ROTATION Uses the algorithm described in Spoor Veldpaus
%1980 and Rose to calculate the Rotation matrix R and translation v of a
%body

aMean = mean(a')'; %Calculate geometric center of all markers in single vector (average of x, y, z columns)
p = zeros(3, size(a,2), size(markerData,3)); %p = position matrix of markers in global coordinate system
pMean = zeros(3, 1, size(markerData,3));
v = pMean;
M = zeros(3, 3, size(markerData,3));
V = M; %eigenvectors of M'M will be stored here
D = M; %eigenvalue^2 of M'M stored here
R = M; %rotation matrix stored here
valid = ones(size(markerData,3),1);

for i = 1 : size(markerData,3)
    p(:,:,i) = markerData(:,1:end-1,i)'; %transpose matrix and remove residual column
    pMean(:,:,i) = mean(p(:,:,i)')'; %Calculate geometric mean of markers, store as column vector    
    M(:,:,i) = 1/size(a,2) * a * p(:,:,i)' - aMean*pMean(:,:,i)';
    if sum(sum(isnan(M(:,:,i)))) == 0
        [V(:,:,i), D(:,:,i)] = eig(M(:,:,i)' * M(:,:,i));
        MV = M(:,:,i) * V(:,:,i);
        v1 = V(:,1,i);
        v2 = V(:,2,i);
        v3 = V(:,3,i);
        d2 = sqrt(D(2,2,i));
        d3 = sqrt(D(3,3,i));
        if abs(sum(cross(v2,v3)./v1) - 3) < 1e-10
            R(:,:,i) = [-cross(MV(:,2),MV(:,3))/(d2*d3), MV(:,2)/d2, MV(:,3)/d3]*V(:,:,i)';
            v(:,:,i) = pMean(:,:,i) - R(:,:,i)*aMean;
        elseif abs(sum(cross(v2,v3)./v1) + 3) < 1e-10 %v2 x v3 = -v1 
            R(:,:,i) = [cross(MV(:,2),MV(:,3))/(d2*d3), MV(:,2)/d2, MV(:,3)/d3]*V(:,:,i)';
            v(:,:,i) = pMean(:,:,i) - R(:,:,i)*aMean;
        else 
            valid(i) = 0;
        end
    else 
        valid(i) = 0;
    end
end

end

